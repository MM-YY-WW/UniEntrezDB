import os.path as osp
import time
from collections import OrderedDict
from math import ceil
from typing import Any, Dict, List, Optional

import lightning.pytorch as pl
import obnb
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric.nn as pygnn
from omegaconf import DictConfig

import obnbench.metrics
from obnbench import optimizers, schedulers
from obnbench.model_layers import feature_encoders, mp_layers
from obnbench.model_layers.post_proc import CorrectAndSmooth, FeaturePropagation

act_register = {
    "relu": nn.ReLU,
    "prelu": nn.PReLU,
    "gelu": nn.GELU,
    "selu": nn.SELU,
    "elu": nn.ELU,
}

norm_register = {
    "BatchNorm": pygnn.BatchNorm,
    "LayerNorm": pygnn.LayerNorm,
    "PairNorm": pygnn.PairNorm,
    "DiffGroupNorm": pygnn.DiffGroupNorm,
    "none": nn.Identity,
}


class ModelModule(pl.LightningModule):

    splits = ["train", "val", "test"]

    def __init__(self, cfg: DictConfig, node_ids: List[str], task_ids: List[str]):
        super().__init__()

        # Save IDs info for logging
        self.node_ids = node_ids
        self.task_ids = task_ids

        self.save_hyperparameters(cfg)  # register cfg as self.hparams
        self.setup_metrics()

        # Build model
        self.feature_encoder = build_feature_encoder(cfg)
        self.mp_layers = build_mp_module(cfg)
        self.pred_head = build_pred_head(cfg)
        self.post_prop, self.pred_act, self.post_cands = build_post_proc(cfg)

    def setup_metrics(self):
        self.metrics = nn.ModuleDict()
        self.final_pred, self.final_true, self.final_scores = {}, {}, []
        metric_kwargs = {
            "task": "multilabel",
            "num_labels": self.hparams._shared.dim_out,
            "validate_args": False,
        }
        for split in self.splits:
            for metric_name in self.hparams.metric.options:
                metric_cls = getattr(obnbench.metrics, metric_name)

                # Metrics to be calculated step/epoch-wise
                self.metrics[f"{split}/{metric_name}"] = \
                    metric_cls(average="macro", **metric_kwargs)

                # Metrics for final evaluation (record perf of individual task)
                self.metrics[f"final/{split}_{metric_name}"] = \
                    metric_cls(average="none", **metric_kwargs)

            self.final_pred[split], self.final_true[split] = [], []

    def configure_optimizers(self):
        optimizer_cls = getattr(optimizers, self.hparams.optim.optimizer)
        optimizer_kwargs = dict(self.hparams.optim.optimizer_kwargs or {})
        if (weight_decay := self.hparams.optim.weight_decay) is not None:
            optimizer_kwargs["weight_decay"] = weight_decay
        optimizer = optimizer_cls(
            self.parameters(),
            lr=self.hparams.optim.lr,
            **optimizer_kwargs,
        )

        lr_scheduler_config = {"optimizer": optimizer}

        if self.hparams.optim.scheduler != "none":
            scheduler_cls = getattr(schedulers, self.hparams.optim.scheduler)
            scheduler_kwargs = dict(self.hparams.optim.scheduler_kwargs or {})

            eval_interval = self.hparams.trainer.eval_interval
            if (patience := scheduler_kwargs.get("patience", None)):
                # Rescale the scheduler patience for ReduceLROnPlateau to the
                # factor w.r.t. the evaluation interval
                scheduler_kwargs["patience"] = ceil(patience / eval_interval)

            scheduler = scheduler_cls(optimizer, **scheduler_kwargs)

            lr_scheduler_config["lr_scheduler"] = {
                "scheduler": scheduler,
                "monitor": f"val/{self.hparams.metric.best}",
                "frequency": eval_interval,
            }

        return lr_scheduler_config

    def forward(self, batch):
        batch = self.feature_encoder(batch)
        batch = self.mp_layers(batch)
        batch = self.pred_head(batch)
        pred, true = self._post_processing(batch)  # FIX: move to end of steps?
        return pred, true

    def _post_processing(self, batch):
        pred, true = batch.x, batch.y

        if self.post_prop is not None:
            pred = self.post_prop(
                pred,
                batch.edge_index,
                batch.edge_weight,
            )

        pred = self.pred_act(pred)

        if not self.training and self.post_cands is not None:
            pred = self.post_cands(
                pred,
                true,
                batch.train_mask.squeeze(-1),
                batch.edge_index,
                batch.edge_weight,
            )

        # Apply split masking
        if batch.split is not None:
            mask = batch[f"{batch.split}_mask"].squeeze(-1)
            pred, true = pred[mask], true[mask]

        return pred, true

    def _shared_step(self, batch, split, final: bool = False):
        # XXX: only allow full batch training now for several reasons (1) post
        # propagation and correction needs access to the full graph, (2) metrics
        # should be computed on the full graph.
        # To resolve (1), need to extract post processing (accordingly the
        # metrics computations also) to a shared private end_of_epoch func.
        # This will also automatically resolve (2).
        tic = time.perf_counter()
        batch.split = split

        # NOTE: split masking is done in _post_processing
        pred, true = self(batch)

        logger_opts = {
            "on_step": False,
            "on_epoch": True,
            "logger": True,
            "batch_size": pred.shape[0],
        }

        # Compute classification loss for training
        if split == "train":
            loss = F.binary_cross_entropy(pred, true)
            self.log(f"{split}/loss", loss.item(), **logger_opts)
            self._maybe_log_grad_norm(logger_opts)
        else:
            loss = None

        self._maybe_log_metrics(pred, true, split, logger_opts)
        self._maybe_prepare_final_metrics(pred, true, split, final)

        self.log(f"{split}/time_epoch", time.perf_counter() - tic, **logger_opts)

        return loss

    # TODO: reset_parameters

    def _maybe_log_grad_norm(self, logger_opts):
        if (not self.training) or (not self.hparams.trainer.watch_grad_norm):
            return

        grad_norms = [
            p.grad.detach().norm(2)
            for p in self.parameters() if p.grad is not None
        ]
        if grad_norms:
            grad_norm = torch.stack(grad_norms).norm(2).item()
            self.log("train/grad_norm", grad_norm, **logger_opts)

    @torch.no_grad()
    def _maybe_log_metrics(self, pred, true, split, logger_opts):
        if (
            (self.current_epoch + 1) % self.hparams.trainer.eval_interval != 0
            and split == "train"
        ):
            return

        for metric_name, metric_obj in self.metrics.items():
            if metric_name.startswith(f"{split}/"):
                metric_obj(pred, true)
                self.log(metric_name, metric_obj, **logger_opts)

    def _maybe_prepare_final_metrics(self, pred, true, split, final):
        if final:
            self.final_pred[split].append(pred.detach())
            self.final_true[split].append(true.detach())

    def _compute_final_metrics(self, split):
        pred = torch.vstack(self.final_pred[split])
        true = torch.vstack(self.final_true[split])

        for metric_name, metric_obj in self.metrics.items():
            prefix = f"final/{split}_"
            if not metric_name.startswith(prefix):
                continue

            score_type = metric_name.replace(prefix, "")
            scores = metric_obj(pred, true).tolist()
            for task_id, task_score in zip(self.task_ids, scores):
                self.final_scores.append(
                    {
                        "split": split,
                        "task_id": task_id,
                        "score_type": score_type,
                        "score_value": task_score,
                    }
                )

            metric_obj.reset()

    def log_final_results(self):
        score_df = pd.DataFrame(self.final_scores)

        for logger in self.loggers:
            if isinstance(logger, pl.loggers.CSVLogger):
                out_path = osp.join(logger.log_dir, "final_scores.csv")
                score_df.to_csv(out_path)
                obnb.logger.info(f"Final results saved to {out_path}")

            elif isinstance(logger, pl.loggers.WandbLogger):
                logger.log_table("final_scores", dataframe=score_df)

            else:
                obnb.logger.error(f"Unknown logger type {type(logger)}")
        return score_df

    def training_step(self, batch, *args, **kwargs):
        return self._shared_step(batch, split="train")

    def validation_step(self, batch, *args, **kwargs):
        self._shared_step(batch, split="val")
        # HACK: Enable early testing that was deliberaly disabled by Lightning
        # https://github.com/Lightning-AI/lightning/issues/5245
        # Note that the early access to testing performance is **not** used for
        # model selection and hyperparameter tuning by any means. Instead, it is
        # only used to see if the trend of the testing and validation curves
        # differ significantly, which indicates there is some problem witht the
        # data split.
        self._shared_step(batch, split="test")

    def on_test_epoch_start(self):
        # Reset final outputs
        for obj in [self.final_pred, self.final_true]:
            for split in obj:
                obj[split].clear()
        self.final_scores.clear()

    def test_step(self, batch, *args, **kwargs):
        for split in self.splits:
            self._shared_step(batch, split=split, final=True)

    def on_test_epoch_end(self, *args, **kwargs):
        for split in self.splits:
            self._compute_final_metrics(split)


class MPModule(nn.Module):

    _residual_opts = [
        "none",
        "skipsum",
        "skipsumbnorm",
        "skipsumlnorm",
        "catlast",
        "catall",
        "catcompact",
    ]

    def __init__(
        self,
        mp_cls: nn.Module,
        dim: int,
        num_layers: int,
        dropout: float = 0.0,
        norm_type: str = "BatchNorm",
        act: str = "relu",
        act_first: bool = False,
        residual_type: str = "none",
        use_edge_feature: bool = True,
        mp_kwargs: Optional[Dict[str, Any]] = None,
        norm_kwargs: Optional[Dict[str, Any]] = None,
    ):
        super().__init__()
        self.dim = dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.norm_type = norm_type
        self.act = act
        self.act_first = act_first

        self.norm_type = norm_type
        self.norm_kwargs = norm_kwargs or {}
        if norm_type == "LayerNorm":
            self.norm_kwargs.setdefault("mode", "node")
        elif norm_type == "DiffGroupNorm":
            self.norm_kwargs.setdefault("groups", 6)
        if norm_type != "PairNorm":
            # Need to pass feature dimension except for PairNorm
            self.norm_kwargs["in_channels"] = None

        # Set residual_type last to make sure we have registered all params
        # Setting residual_type will automatically set up the forward function
        # and the dimensions for the hidden layers.
        self.residual_type = residual_type

        # Set up the message passing layers using the dimensions prepared
        self.layers = nn.ModuleList()
        self.res_norms = nn.ModuleList()
        mp_kwargs = mp_kwargs or {}
        for dim_in in self._layer_dims:
            self._build_layer(mp_cls, dim_in, dim, use_edge_feature, mp_kwargs)

    def _build_layer(
        self,
        mp_cls: nn.Module,
        dim_in: int,
        dim_out: int,
        use_edge_feature,
        mp_kwargs: Dict[str, Any],
    ) -> nn.Module:
        conv_layer = mp_cls(dim_in, dim_out, use_edge_feature, **mp_kwargs)
        activation = act_register.get(self.act)()
        dropout = nn.Dropout(self.dropout)

        # Check if 'in_channels' is set to determine whether we need to pass
        # the feature dimension to initialize the normalization layer
        if "in_channels" in self.norm_kwargs:
            self.norm_kwargs["in_channels"] = dim_out
        norm_layer = norm_register.get(self.norm_type)(**self.norm_kwargs)

        # Graph convolution layer
        new_layer = nn.ModuleDict()
        new_layer["conv"] = conv_layer

        # Post convolution layers
        post_conv = []
        if self.act_first:
            post_conv.extend([("act", activation), "norm", norm_layer])
        else:
            post_conv.extend([("norm", norm_layer), ("act", activation)])
        post_conv.append(("dropout", dropout))
        new_layer["post_conv"] = nn.Sequential(OrderedDict(post_conv))
        self.layers.append(new_layer)

        # Residual normalizations
        if self.residual_type == "skipsumbnorm":
            self.res_norms.append(norm_register["BatchNorm"](dim_out))
        elif self.residual_type == "skipsumlnorm":
            self.res_norms.append(norm_register["LayerNorm"](dim_out, mode="node"))

    def extra_repr(self) -> str:
        return f"residual_type: {self.residual_type}"

    @property
    def residual_type(self) -> str:
        return self._residual_type

    @residual_type.setter
    def residual_type(self, val: str):
        if val == "none":
            self._forward = self._stack_forward
            self._layer_dims = [self.dim] * self.num_layers
        elif val in ["skipsum", "skipsumbnorm", "skipsumlnorm"]:
            self._forward = self._skipsum_forward
            self._layer_dims = [self.dim] * self.num_layers
        elif val == "catlast":
            self._forward = self._catlast_forward
            self._layer_dims = (
                [self.dim] * (self.num_layers - 1)
                + [self.dim * self.num_layers]
            )
        elif val == "catall":
            self._forward = self._catall_forward
            self._layer_dims = [self.dim * (i + 1) for i in range(self.num_layers)]
        else:
            raise ValueError(
                f"Unknown residual type {val!r}, available options are:\n"
                f"    {self._residual_opts}",
            )
        self._residual_type = val

    @staticmethod
    def _layer_forward(layer, batch):
        batch = layer["conv"](batch)
        batch.x = layer["post_conv"](batch.x)
        return batch

    def _stack_forward(self, batch):
        for layer in self.layers:
            batch = self._layer_forward(layer, batch)
        return batch

    def _skipsum_forward(self, batch):
        for i, layer in enumerate(self.layers):
            x_prev = batch.x
            batch = layer["conv"](batch)
            # batch = self._layer_forward(layer, batch)
            if self.res_norms:
                batch.x = self.res_norms[i](batch.x)
            batch.x = batch.x + x_prev
            batch.x = layer["post_conv"](batch.x)
        return batch

    def _catlast_forward(self, batch):
        xs = []
        for i, layer in enumerate(self.layers):
            if i == self.num_layers - 1:
                batch.x = torch.cat([batch.x] + xs, dim=1)
            batch = self._layer_forward(layer, batch)
            if i < self.num_layers - 1:
                xs.append(batch.x)
        return batch

    def _catall_forward(self, batch):
        for i, layer in enumerate(self.layers):
            x_prev = batch.x
            batch = self._layer_forward(layer, batch)
            if i < self.num_layers - 1:
                batch.x = torch.cat([batch.x, x_prev], dim=1)
        return batch

    def forward(self, batch):
        return self._forward(batch)


class PredictionHeadModule(nn.Module):

    def __init__(
        self,
        dim_in: int,
        dim_out: int,
        num_layers: int = 1,
        dim_inner: int = 128,
    ):
        super().__init__()
        if num_layers > 0:
            self.layers = pygnn.MLP(
                in_channels=dim_in,
                out_channels=dim_out,
                hidden_channels=dim_inner,
                num_layers=num_layers,
                act="relu",
                norm="batch_norm",
                bias=True,
                plain_last=True,
            )
        else:
            self.layers = nn.Identity()

    def forward(self, batch):
        batch.x = self.layers(batch.x)
        return batch


def build_feature_encoder(cfg: DictConfig):
    feat_names = cfg.dataset.node_encoders.split("+")

    fe_list = []
    for i, feat_name in enumerate(feat_names):
        fe_cfg = cfg.node_encoder_params.get(feat_name)
        fe_cls = getattr(feature_encoders, f"{feat_name}FeatureEncoder")
        fe = fe_cls(
            dim_feat=cfg._shared.fe_raw_dims[i],
            dim_encoder=cfg.model.hid_dim,
            layers=fe_cfg.layers,
            dropout=fe_cfg.dropout,
            raw_dropout=fe_cfg.raw_dropout,
            raw_bn=fe_cfg.raw_bn,
            num_nodes=cfg._shared.num_nodes,
        )
        fe_list.append(fe)

    if len(fe_list) == 1:
        return fe_list[0]
    else:
        fe_cfg = cfg.node_encoder_params.Composed
        return feature_encoders.ComposedFeatureEncoder(
            dim_feat=cfg._shared.composed_fe_dim_in,
            dim_encoder=cfg.model.hid_dim,
            layers=fe_cfg.layers,
            dropout=fe_cfg.dropout,
            raw_dropout=fe_cfg.raw_dropout,
            raw_bn=fe_cfg.raw_bn,
            fe_list=fe_list,
        )


def build_mp_module(cfg: DictConfig):
    mp_cls = getattr(mp_layers, cfg.model.mp_type)
    return MPModule(
        mp_cls,
        dim=cfg.model.hid_dim,
        num_layers=cfg.model.mp_layers,
        dropout=cfg.model.dropout,
        norm_type=cfg.model.norm_type,
        act=cfg.model.act,
        act_first=cfg.model.act_first,
        residual_type=cfg.model.residual_type,
        use_edge_feature=cfg.model.use_edge_feature,
        mp_kwargs=cfg.model.mp_kwargs,
        norm_kwargs=cfg.model.norm_kwargs,
    )


def build_pred_head(cfg: DictConfig):
    return PredictionHeadModule(
        dim_in=cfg._shared.pred_head_dim_in,
        dim_out=cfg._shared.dim_out,
        num_layers=cfg.model.pred_head_layers,
        dim_inner=cfg.model.hid_dim,
    )


def build_post_proc(cfg: DictConfig):
    post_prop = None
    if cfg.model.post_prop.enable:
        post_prop = FeaturePropagation(
            num_layers=cfg.model.post_prop.num_layers,
            alpha=cfg.model.post_prop.alpha,
            norm=cfg.model.post_prop.norm,
            cached=cfg.model.post_prop.cached,
        )

    pred_act = nn.Identity() if cfg.model.skip_pred_act else nn.Sigmoid()

    post_cands = None
    if cfg.model.post_cands.enable:
        post_cands = CorrectAndSmooth(
            num_correction_layers=cfg.model.post_cands.num_correction_layers,
            num_smoothing_layers=cfg.model.post_cands.num_smoothing_layers,
            correction_alpha=cfg.model.post_cands.correction_alpha,
            smoothing_alpha=cfg.model.post_cands.smoothing_alpha,
            correction_norm=cfg.model.post_cands.correction_norm,
            smoothing_norm=cfg.model.post_cands.smoothing_norm,
            autoscale=cfg.model.post_cands.autoscale,
            scale=cfg.model.post_cands.scale,
            cached=cfg.model.post_cands.cached,
        )

    return post_prop, pred_act, post_cands
