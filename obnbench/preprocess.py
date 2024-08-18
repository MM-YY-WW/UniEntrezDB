import time
from functools import wraps
from pprint import pformat
from typing import Any, Dict, Optional

import obnb
import numpy as np
import scipy.sparse as sp
import torch
from obnb.ext.attnwalk import attnwalk_embed
from obnb.ext.grape import grape_embed
from obnb.ext.orbital_features import orbital_feat_extract
from obnb.ext.pecanpy import pecanpy_embed
from obnb.graph import SparseGraph
from obnb.util.logger import display_pbar
from omegaconf import DictConfig, open_dict
from sklearn.decomposition import PCA
from sklearn.preprocessing import KBinsDiscretizer
from sklearn.random_projection import GaussianRandomProjection, SparseRandomProjection
from torch_geometric.data import Dataset

from obnbench.model_layers import feature_encoders

precomp_func_register = {}


class PreCompFeatureWrapper:

    def __init__(self, name: str):
        self.name = name
        self.feat_name = f"rawfeat_{name}"
        self.fe_name = f"{name}FeatureEncoder"
        assert hasattr(feature_encoders, self.fe_name)

    def __call__(self, func):

        @wraps(func)
        def wrapped_func(dataset: Dataset, *args, **kwargs) -> Dataset:
            obnb.logger.info(f"Precomputing raw features for {self.fe_name}")
            feat = func(*args, dataset=dataset, **kwargs)

            if isinstance(feat, np.ndarray):
                feat = torch.from_numpy(feat.astype(np.float32))
            elif feat is None:
                return dataset
            elif not isinstance(feat, torch.Tensor):
                raise TypeError(
                    f"Unknown feature type {type(feat)} returned from the "
                    f"preprocessing function {func}",
                )

            # Handle dataset attr
            dataset._data_list = None
            dataset._data[self.feat_name] = feat
            if dataset.slices is not None:
                dataset.slices[self.feat_name] = torch.LongTensor([0, feat.shape[1]])

            return dataset

        precomp_func_register[self.name] = wrapped_func

        return wrapped_func


@PreCompFeatureWrapper("OneHotLogDeg")
def get_onehot_logdeg(feat_dim: int, adj: np.ndarray, **kwargs) -> np.ndarray:
    log_deg = np.log(adj.sum(axis=1, keepdims=True))
    feat = KBinsDiscretizer(
        n_bins=feat_dim,
        encode="onehot-dense",
        strategy="uniform",
    ).fit_transform(log_deg)
    obnb.logger.info(f"Bins stats:\n{feat.sum(0)}")
    return feat


@PreCompFeatureWrapper("Constant")
def get_const(feat_dim: int, adj: np.ndarray, **kwargs) -> np.ndarray:
    if feat_dim != 1:
        raise ValueError(
            "Constant feature only allows dimension of 1, "
            f"got {feat_dim!r}",
        )
    feat = np.ones((adj.shape[0], 1))
    return feat


@PreCompFeatureWrapper("RandomNormal")
def get_random_normal(
    feat_dim: int,
    adj: np.ndarray,
    random_state: Optional[int] = None,
    **kwargs,
) -> np.ndarray:
    feat = np.random.default_rng(random_state).random((adj.shape[0], feat_dim))
    return feat



@PreCompFeatureWrapper("Orbital")
def get_orbital_counts(
    g: SparseGraph,
    num_workers: int = 1,
    feat_kwargs: Optional[Dict[str, Any]] = None,
    show_progress: bool = True,
    **kwargs,
) -> np.ndarray:
    feat_kwargs = feat_kwargs or {}
    feat_kwargs.setdefault("graphlet_size", 3)
    feat = orbital_feat_extract(
        g,
        n_jobs=num_workers,
        as_array=True,
        verbose=show_progress,
        **feat_kwargs,
    )
    return feat


@PreCompFeatureWrapper("SVD")
def get_svd_emb(feat_dim: int, adj: np.ndarray, **kwargs) -> np.ndarray:
    A = sp.csr_matrix(adj)
    feat, _, _ = sp.linalg.svds(A, k=feat_dim, which="LM")
    return feat


@PreCompFeatureWrapper("LapEigMap")
def get_lap_eig_map(feat_dim: int, adj: np.ndarray, **kwargs) -> np.ndarray:
    L = sp.csr_matrix(np.diag(adj.sum(1)) - adj)
    assert (L != L.T).sum() == 0, "The input network must be undirected."

    # Symmetric normalized graph Laplacian
    D_inv_sqrt = sp.diags(1 / np.sqrt(adj.sum(1)))
    L = D_inv_sqrt @ L @ D_inv_sqrt

    eigvals, eigvecs = sp.linalg.eigsh(L, which="SM", k=feat_dim + 1)
    sorted_idx = eigvals.argsort()
    eigvals = eigvals[sorted_idx]
    eigvecs = eigvecs[:, sorted_idx]

    assert (eigvals[1:] > 1e-8).all(), f"Network appears to be disconnected.\n{eigvals=}"
    feat = eigvecs[:, 1:] / np.sqrt((eigvecs[:, 1:] ** 2).sum(0))  # l2 normalize

    return feat


@PreCompFeatureWrapper("RandomWalkDiag")
def get_rand_walk_diag(feat_dim: int, adj: np.ndarray, **kwargs) -> np.ndarray:
    P = adj / adj.sum(0)
    feat = np.zeros((adj.shape[0], feat_dim))
    vec = np.ones(adj.shape[0])
    for i in range(feat_dim):
        vec = P @ vec
        feat[:, i] = vec
    return feat


@PreCompFeatureWrapper("RandProjGaussian")
def get_rand_proj_gaussian(
    feat_dim: int,
    adj: np.ndarray,
    random_state: Optional[int] = None,
    **kwargs,
) -> np.ndarray:
    grp = GaussianRandomProjection(n_components=feat_dim, random_state=random_state)
    feat = grp.fit_transform(adj)
    return feat


@PreCompFeatureWrapper("RandProjSparse")
def get_rand_proj_sparse(
    feat_dim: int,
    adj: np.ndarray,
    random_state: Optional[int] = None,
    **kwargs,
) -> np.ndarray:
    srp = SparseRandomProjection(
        n_components=feat_dim,
        dense_output=True,
        random_state=random_state,
    )
    feat = srp.fit_transform(adj)
    return feat


@PreCompFeatureWrapper("LINE1")
def get_line1_emb(
    feat_dim: int,
    g: SparseGraph,
    random_state: Optional[int] = None,
    show_progress: bool = True,
    **kwargs,
) -> np.ndarray:
    feat = grape_embed(
        g,
        "FirstOrderLINEEnsmallen",
        dim=feat_dim,
        as_array=True,
        random_state=random_state,
        verbose=show_progress,
    )
    return feat


@PreCompFeatureWrapper("LINE2")
def get_line2_emb(
    feat_dim: int,
    g: SparseGraph,
    random_state: Optional[int] = None,
    show_progress: bool = True,
    **kwargs,
) -> np.ndarray:
    feat = grape_embed(
        g,
        "SecondOrderLINEEnsmallen",
        dim=feat_dim,
        as_array=True,
        random_state=random_state,
        verbose=show_progress
    )
    return feat


@PreCompFeatureWrapper("Node2vec")
# def get_n2v_emb(
#     feat_dim: int,
#     g: SparseGraph,
#     num_workers: int = 1,
#     random_state: Optional[int] = None,
#     show_progress: bool = True,
#     feat_kwargs: Optional[Dict[str, Any]] = None,
#     **kwargs,
# ) -> np.ndarray:
#     feat = pecanpy_embed(
#         g,
#         workers=num_workers,
#         verbose=show_progress,
#         dim=feat_dim,
#         as_array=True,
#         random_state=random_state,
#         **feat_kwargs,
#     )
#     print(feat.shape)
def get_emb(feat_dim: int,
    g: SparseGraph,
    num_workers: int = 1,
    random_state: Optional[int] = None,
    show_progress: bool = True,
    feat_kwargs: Optional[Dict[str, Any]] = None,
    **kwargs,)-> np.ndarray:
    # path = "datasets/ProteomeHD/processed/"
    path = "datasets/BioGRID/processed/"
    # path = "datasets/SIGNOR/processed/"
    # path = path + "goa.pt"
    # path = path + "goa_dnabert.pt"
    # path = path + "goa_gene2vec.pt"
    path = path + "gene2vec.pt"
    # path = path + "ontoprotein.pt"
    # path = path + "dnabert.pt"

    feat = torch.load(path).numpy()
    return feat


@PreCompFeatureWrapper("Walklets")
def get_walklets_emb(
    feat_dim: int,
    g: SparseGraph,
    random_state: Optional[int] = None,
    feat_kwargs: Optional[Dict[str, Any]] = None,
    **kwargs,
) -> np.ndarray:
    feat_kwargs = feat_kwargs or {}
    feat_kwargs.setdefault("epochs", 1)
    feat_kwargs.setdefault("window_size", 4)

    # NOTE: The resulding feat is a concatenation of (window_size x 2) number
    # of embeddings, each has the dimension of feat_dim.
    feat_raw = grape_embed(
        g,
        "WalkletsSkipGramEnsmallen",
        dim=feat_dim * feat_kwargs["window_size"],  # one emb per-window (both-sides)
        as_array=True,
        grape_enable=True,
        random_state=random_state,
        **feat_kwargs,
    )

    # Reduce multscale embedding to feat_dim via PCA following arxiv:1605.02115
    if feat_raw.shape[1] > feat_dim:
        pca = PCA(n_components=feat_dim, random_state=random_state)
        feat = pca.fit_transform(feat_raw)
        evr = pca.explained_variance_ratio_.sum()
        obnb.logger.info(
            "Reduced concatenated walklets embedding dimensions from "
            f"{feat_raw.shape[1]} to {feat.shape[1]} (EVR={evr:.2%})."
        )
    else:
        feat = feat_raw

    return feat


@PreCompFeatureWrapper("AttnWalk")
def get_attnwalk_emb(
    feat_dim: int,
    g: SparseGraph,
    show_progress: bool = True,
    feat_kwargs: Optional[Dict[str, Any]] = None,
    **kwargs,
) -> np.ndarray:
    feat, attn = attnwalk_embed(
        g,
        verbose=show_progress,
        dim=feat_dim,
        as_array=True,
        return_attn=True,
        **feat_kwargs
    )
    attn_str = ", ".join(f"{i:.4f}" for i in attn)
    obnb.logger.info(f"AttnWalk attentions: [{attn_str}]")
    return feat


@PreCompFeatureWrapper("Adj")
def get_adj(adj: np.ndarray, **kwargs) -> np.ndarray:
    return adj.copy()


@PreCompFeatureWrapper("AdjEmbBag")
def get_adj_emb_bag(**kwargs) -> None:
    # No need to compute new features
    ...


@PreCompFeatureWrapper("Embedding")
def get_embedding(**kwargs) -> None:
    # No need to compute new features
    ...


@PreCompFeatureWrapper("LabelReuse")
def get_label_resuse(dataset: Dataset, **kwargs) -> torch.Tensor:
    feat = torch.zeros_like(dataset._data.y, dtype=torch.float)
    train_mask = dataset._data.train_mask[:, 0]
    feat[train_mask] = dataset._data.y[train_mask]
    feat /= feat.sum(0)  # normalize
    return feat


def precompute_features(cfg: DictConfig, dataset: Dataset, g: SparseGraph):
    # Catch invalid node encoders before executing
    invalid_fe = []
    for feat_name in (node_encoders := cfg.dataset.node_encoders.split("+")):
        if feat_name not in precomp_func_register:
            invalid_fe.append(feat_name)
    if invalid_fe:
        raise ValueError(
            f"Invalid node encoders {invalid_fe} in {cfg.dataset.node_encoders}",
        )

    # Prepare shared data arguments
    data_dict = {
        "dataset": dataset,
        "g": g,
        "adj": g.to_dense_graph().mat,
        "num_workers": cfg.num_workers,
        "show_progress": display_pbar(cfg.log_level),
        "random_state": cfg.seed,
    }

    tic = time.perf_counter()
    obnb.logger.info("Start pre-computing features")
    for feat_name in node_encoders:
        fe_cfg = cfg.node_encoder_params.get(feat_name)
        feat_dim = fe_cfg.get("raw_dim", None)
        feat_kwargs = dict(fe_cfg.get("feat_kwargs", None) or {})
        precomp_func_register[feat_name](
            feat_dim=feat_dim,
            log_level=cfg.log_level,
            feat_kwargs=feat_kwargs,
            **data_dict,
        )
    elapsed = time.perf_counter() - tic
    obnb.logger.info(f"Precomputation done! Took {elapsed:.2f} seconds.")


def infer_dimensions(cfg: DictConfig, dataset: Dataset):
    # NOTE: this function is called after precompute_features, so we don't need
    # to check the validity of the node_encoders setting again.
    node_encoders = cfg.dataset.node_encoders.split("+")

    # Infer number of nodes and tasks
    num_nodes, dim_out = dataset._data.y.shape

    # Infer feature encoder dimensions
    hid_dim = cfg.model.hid_dim
    fe_raw_dims, fe_processed_dims = [], []
    for feat_name in node_encoders:
        fe_cfg = cfg.node_encoder_params.get(feat_name)

        # Handel special cases, otherwise get raw dim from precomputed features
        if feat_name in ["AdjEmbBag", "Embedding"]:
            raw_dim = fe_cfg.raw_dim
        else:
            raw_dim = dataset._data[f"rawfeat_{feat_name}"].shape[1]
        encoded_dim = raw_dim if fe_cfg.layers == 0 else hid_dim

        fe_raw_dims.append(raw_dim)
        fe_processed_dims.append(encoded_dim)

    # Infer composed feature dimension and message passing input dimension
    if len(node_encoders) == 1:  # single feature encoder
        composed_fe_dim_in = None
        mp_dim_in = fe_processed_dims[0]
    else:  # composed feature encoder
        composed_fe_dim_in = sum(fe_processed_dims)
        fe_cfg = cfg.node_encoder_params.Composed
        mp_dim_in = composed_fe_dim_in if fe_cfg.layers == 0 else hid_dim

    # Infer prediction head input dimension
    pred_head_dim_in = mp_dim_in if cfg.model.mp_layers == 0 else hid_dim

    inferred_dims_dict = {
        "fe_raw_dims": fe_raw_dims,
        "fe_processed_dims": fe_processed_dims,
        "composed_fe_dim_in": composed_fe_dim_in,
        "mp_dim_in": mp_dim_in,
        "pred_head_dim_in": pred_head_dim_in,
        "dim_out": dim_out,
        "num_nodes": num_nodes,
    }
    obnb.logger.info(f"Node encoders: {node_encoders}")
    obnb.logger.info(f"Inferred module dimensions:\n{pformat(inferred_dims_dict)}")

    with open_dict(cfg):
        cfg._shared = inferred_dims_dict
