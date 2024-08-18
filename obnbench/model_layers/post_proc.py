import warnings
from typing import Optional, Tuple

import torch
import torch.nn.functional as F
from torch import Tensor
from torch_geometric import nn as pygnn
from torch_geometric.typing import SparseTensor
from torch_geometric.utils import scatter, spmm


class FeaturePropagation(pygnn.MessagePassing):
    """Feature propagation similar to APPNP.

    Args:
        num_layers: Number of propagation steps. If set to None, then will run
            till convergence (element wise mean absolute error below tol) or
            max_iter reached.
        norm: Adjacency matrix normalization schemes ("left", "right", or
            "sym").
        dropout: Dropout to propagated features at each propagation step.
        cached: If set to True, then save the computed normalized adjacency
            for future usage. Otherwise, compute the noramlized adjacency
            everytime the forward function is called.
        max_iter: Maximum iteration for the convergence mode.
        tol: Tolerence for the convergence mode.

    """

    def __init__(
        self,
        num_layers: Optional[int],
        alpha: float,
        norm: str = "left",
        dropout: float = 0.0,
        cached: bool = False,
        max_iter: int = 100,
        tol: float = 1e-12,
    ):
        super().__init__(aggr="add", flow="source_to_target")

        if num_layers is not None:
            tol = None  # disable convergence mode
        else:
            num_layers = max_iter  # enable convergence mode

        self.num_layers = num_layers
        self.alpha = alpha
        self.dropout = dropout
        self.norm = norm
        self.cached = cached
        self.tol = tol

    @property
    def cached_norm_edge_weight(self) -> Optional[Tensor]:
        return getattr(self, "_cached_norm_edge_weight", None)

    def forward(
        self,
        x: Tensor,
        edge_index: Tensor,
        edge_weight: Optional[Tensor],
    ):
        # Try to load normalized adjacency edge weights from cache
        # WARNING: Caching assumes full-batch training with no shuffling.
        if (norm_edge_weight := self.cached_norm_edge_weight) is None:
            norm_edge_weight = adj_norm(
                self.norm,
                x.shape[0],
                edge_index,
                edge_weight,
            )

            if self.cached:
                self._cached_norm_edge_weight = norm_edge_weight

        x0 = x
        for i in range(self.num_layers):
            x_prev = F.dropout(x, p=self.dropout, training=self.training)

            # propagate_type: (x: Tensor, edge_ewight: Tensor)
            x = self.propagate(edge_index, x=x_prev, edge_weight=norm_edge_weight)
            x = self.alpha * x + (1 - self.alpha) * x0

            # Check convergence
            if self.tol is not None:
                with torch.no_grad():
                    err = (x - x_prev).abs().mean()
                    if err < self.tol:
                        # print(f"FeaturePropagation converged at {i} step")
                        break
        else:
            # Check if convergence failed
            if self.tol is not None:
                warnings.warn(
                    "FeaturePropagation failed to converge within "
                    f"{self.num_layers} iterations (tol={self.tol:.2e}, "
                    f"err={err:.2e}). Consider increasing the max_iter.",
                    RuntimeWarning,
                    stacklevel=2,
                )

        return x

    def message(self, x_j: Tensor, edge_weight: Tensor) -> Tensor:
        return edge_weight.view(-1, 1) * x_j

    def message_and_aggregate(self, adj_t: SparseTensor, x: Tensor) -> Tensor:
        return spmm(adj_t, x, reduce=self.aggr)


class CorrectAndSmooth(pygnn.CorrectAndSmooth):
    """Correct and smooth with improved forward to simplify calling.

    Normalized adjacency used for propagation is computed on-the-fly upon
    calling the forward function. Furthermore, three different normalization
    schemes ("left", "right", and "sym") can be selected independently for the
    correction and smoothing steps. Optionally, the computed normalized
    adjacencies can be cached for future computation.

    """

    def __init__(
        self,
        num_correction_layers: int,
        correction_alpha: float,
        num_smoothing_layers: int,
        smoothing_alpha: float,
        autoscale: bool = True,
        scale: float = 1.0,
        correction_norm: str = "sym",
        smoothing_norm: str = "left",
        cached: bool = False,
    ):
        super().__init__(
            num_correction_layers=num_correction_layers,
            correction_alpha=correction_alpha,
            num_smoothing_layers=num_smoothing_layers,
            smoothing_alpha=smoothing_alpha,
            autoscale=autoscale,
            scale=scale,
        )
        self.correction_norm = correction_norm
        self.smoothing_norm = smoothing_norm
        self.cached = cached

    @property
    def cached_correct_edge_weight(self) -> Optional[Tensor]:
        return getattr(self, "_cached_correct_edge_weight", None)

    @property
    def cached_smooth_edge_weight(self) -> Optional[Tensor]:
        return getattr(self, "_cached_smooth_edge_weight", None)

    def forward(
        self,
        y_soft: Tensor,
        y_true: Tensor,
        mask: Tensor,
        edge_index: Tensor,
        edge_weight: Optional[Tensor] = None,
    ):
        correct_edge_weight, smooth_edge_weight = self._get_adjs(
            y_soft.shape[0],
            edge_index,
            edge_weight,
        )
        y_true_masked = y_true[mask]
        y_soft = self.correct(y_soft, y_true_masked, mask, edge_index, correct_edge_weight)
        y_soft = self.smooth(y_soft, y_true_masked, mask, edge_index, smooth_edge_weight)
        return y_soft

    def _get_adjs(
        self,
        num_nodes: int,
        edge_index: Tensor,
        edge_weight: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tensor]:
        # Try to load from cache first
        # WARNING: Caching assumes full-batch training with no shuffling.
        correct_edge_weight = self.cached_correct_edge_weight
        smooth_edge_weight = self.cached_smooth_edge_weight

        # Compute if cache not available
        if correct_edge_weight is None or smooth_edge_weight is None:
            if edge_weight is None:
                edge_weight = torch.ones_like(edge_index[0], dtype=torch.float)

            deg_inv = get_degree(num_nodes, edge_index, edge_weight).pow(-1)

            args = (num_nodes, edge_index, edge_weight, deg_inv)
            correct_edge_weight = adj_norm(self.correction_norm, *args)
            if self.smoothing_norm == self.correction_norm:
                smooth_edge_weight = correct_edge_weight
            else:
                smooth_edge_weight = adj_norm(self.smoothing_norm, *args)

            if self.cached:
                self._cached_correct_edge_weight = correct_edge_weight
                self._cached_smooth_edge_weight = smooth_edge_weight

        return correct_edge_weight, smooth_edge_weight


def get_degree(
    num_nodes: int,
    edge_index: Tensor,
    edge_weight: Optional[Tensor] = None,
) -> Tensor:
    if edge_weight is None:
        edge_weight = torch.ones_like(edge_index[0], dtype=torch.float)

    row, col = edge_index
    deg = scatter(edge_weight, col, dim=0, dim_size=num_nodes, reduce="sum")

    return deg


def adj_norm(
    norm_type: str,
    num_nodes: int,
    edge_index: Tensor,
    edge_weight: Optional[Tensor],
    deg_inv: Optional[Tensor] = None,
) -> Tensor:
    if edge_weight is None:
        edge_weight = torch.ones_like(edge_index[0], dtype=torch.float)

    if deg_inv is None:
        deg_inv = get_degree(num_nodes, edge_index, edge_weight).pow(-1)

    row, col = edge_index
    if norm_type == "sym":
        deg_inv_sqrt = deg_inv.pow(0.5)
        edge_weight = deg_inv_sqrt[row] * edge_weight * deg_inv_sqrt[col]
    elif norm_type == "left":
        edge_weight = deg_inv[row] * edge_weight
    elif norm_type == "right":
        edge_weight = edge_weight * deg_inv[col]
    else:
        raise ValueError(f"Unknown norm_type {norm_type!r}")

    return edge_weight
