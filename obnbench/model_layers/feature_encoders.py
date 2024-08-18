from typing import List

import torch
import torch.nn as nn
import torch_geometric.nn as pygnn

from obnbench.model_layers.misc import RawFeatNorm


class BaseFeatureEncoder(nn.Module):

    def __init__(
        self,
        *,
        dim_feat: int,
        dim_encoder: int = 128,
        layers: int = 1,
        dropout: float = 0.0,
        raw_dropout: float = 0.0,
        raw_bn: bool = True,
        raw_bn_affine: bool = True,
        **kwargs,
    ):
        super().__init__()

        self.dim_feat = dim_feat
        self.dim_encoder = dim_encoder

        self.raw_bn = RawFeatNorm(dim_feat, affine=raw_bn_affine) if raw_bn else None
        self.raw_dropout = nn.Dropout(raw_dropout)

        if layers > 0:
            self.enc = pygnn.MLP(
                in_channels=dim_feat,
                out_channels=dim_encoder,
                hidden_channels=dim_encoder,
                num_layers=layers,
                dropout=dropout,
                norm="batch_norm",
                act="relu",
                plain_last=False,
            )
        else:
            self.enc = None

        fe_name = self.__class__.__name__.replace("FeatureEncoder", "")
        self._raw_feat_name = f"rawfeat_{fe_name}"

    def reset_parameters(self):
        if self.raw_bn is not None:
            self.raw_bn.reset_parameters()
        if self.enc is not None:
            self.enc.reset_parameters()

    def get_raw_features(self, batch) -> torch.Tensor:
        return getattr(batch, self._raw_feat_name)

    def forward(self, batch):
        x = self.get_raw_features(batch)
        if self.raw_bn is not None:
            x = self.raw_bn(x)
        x = self.raw_dropout(x)
        batch.x = x if self.enc is None else self.enc(x)
        return batch


class ComposedFeatureEncoder(BaseFeatureEncoder):

    def __init__(self, *, fe_list: List[BaseFeatureEncoder], **kwargs):
        super().__init__(**kwargs)
        if len(fe_list) < 1:
            raise ValueError("Empty sub-feature-encoder list.")
        self.feature_encoders = nn.ModuleList(fe_list)

    def get_raw_features(self, batch):
        xs = [enc(batch).x for enc in self.feature_encoders]
        return torch.cat(xs, dim=1)


class OneHotLogDegFeatureEncoder(BaseFeatureEncoder):
    ...


class ConstantFeatureEncoder(BaseFeatureEncoder):
    ...


class RandomNormalFeatureEncoder(BaseFeatureEncoder):
    ...


class OrbitalFeatureEncoder(BaseFeatureEncoder):
    ...


class SVDFeatureEncoder(BaseFeatureEncoder):
    ...


class LapEigMapFeatureEncoder(BaseFeatureEncoder):
    ...


class RandomWalkDiagFeatureEncoder(BaseFeatureEncoder):
    ...


class RandProjGaussianFeatureEncoder(BaseFeatureEncoder):
    ...


class RandProjSparseFeatureEncoder(BaseFeatureEncoder):
    ...


class LINE1FeatureEncoder(BaseFeatureEncoder):
    ...


class LINE2FeatureEncoder(BaseFeatureEncoder):
    ...


class Node2vecFeatureEncoder(BaseFeatureEncoder):
    ...


class WalkletsFeatureEncoder(BaseFeatureEncoder):
    ...


class AttnWalkFeatureEncoder(BaseFeatureEncoder):
    ...


class AdjFeatureEncoder(BaseFeatureEncoder):
    ...


class AdjEmbBagFeatureEncoder(BaseFeatureEncoder):

    def __init__(self, *, num_nodes: int, dim_feat: int, bias: bool = True, **kwargs):
        super().__init__(dim_feat=dim_feat, **kwargs)
        self.num_nodes = num_nodes
        self.emb = nn.EmbeddingBag(num_nodes, dim_feat, mode="sum")
        self.bias = nn.Parameter(torch.empty(dim_feat)) if bias else None
        self.reset_parameters()

    def reset_parameters(self):
        super().reset_parameters()
        self.emb.reset_parameters()
        if self.bias is not None:
            nn.init.zeros_(self.bias)

    def get_raw_features(self, batch) -> torch.Tensor:
        row, col = batch.edge_index

        # Offset by degree
        ones = torch.ones_like(row, dtype=torch.long)
        offsets = torch.zeros(self.num_nodes + 1, dtype=torch.long, device=row.device)
        offsets = offsets.scatter_add(0, row + 1, ones).cumsum(0)

        # Target to source aggregation
        x = self.emb(col, offsets[:-1], batch.edge_weight)

        if self.bias is not None:
            x = x + self.bias

        return x


class EmbeddingFeatureEncoder(BaseFeatureEncoder):

    def __init__(self, *, num_nodes: int, dim_feat: int, **kwargs):
        super().__init__(dim_feat=dim_feat, **kwargs)
        self.num_nodes = num_nodes
        self.emb = nn.Embedding(num_nodes, dim_feat)

    def reset_parameters(self):
        super().reset_parameters()
        self.emb.reset_parameters()

    def get_raw_features(self, batch) -> torch.Tensor:
        # XXX: assumes full-batch
        return self.emb.weight


class LabelReuseFeatureEncoder(BaseFeatureEncoder):
    ...


__all__ = [
    "AdjFeatureEncoder",
    "ComposedFeatureEncoder",
    "ConstantFeatureEncoder",
    "LINE1FeatureEncoder",
    "LINE2FeatureEncoder",
    "LabelReuseFeatureEncoder",
    "LapEigMapFeatureEncoder",
    "Node2vecFeatureEncoder",
    "OrbitalFeatureEncoder",
    "RandProjGaussianFeatureEncoder",
    "RandProjSparseFeatureEncoder",
    "RandomNormalFeatureEncoder",
    "RandomWalkDiagFeatureEncoder",
    "SVDFeatureEncoder",
    "WalkletsFeatureEncoder",
]
