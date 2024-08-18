import warnings
from typing import Optional

import torch_geometric.nn as pygnn

# TODO: GAT/GINE edge attr mod bining


class BaseConvMixin:

    _edge_usage = "none"

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        use_edge_feature: bool = True,
        **kwargs,
    ):
        # Set up forward function dependending on the edge usage capabilities
        # of the wrapped convolution module
        if self._edge_usage == "none" or use_edge_feature:
            self._forward = self._forward_simple
        elif self._edge_usage == "edge_weight":
            self._forward = self._forward_edgeweight
        elif self._edge_usage == "edge_attr":
            self._forward = self._forward_edgeattr
        else:
            raise ValueError(
                f"Unknown edge usage mode {self._edge_usage!r}, "
                "available options are: 'none', 'edge_weight', 'edge_attr'",
            )

        super().__init__(in_channels, out_channels, **kwargs)

    def _forward_simple(self, batch):
        return super().forward(batch.x, batch.edge_index)

    def _forward_edgeweight(self, batch):
        return super().forward(batch.x, batch.edge_index, edge_weight=batch.edge_weight)

    def _forward_edgeattr(self, batch):
        if (edge_attr := batch.edge_attr) is None:
            warnings.warn(
                "Implicitly use edge_attr in place of edge_weight because "
                "edge_attr is unavailable.",
                stacklevel=2,
            )
            # Try to use edge weight as attr if edge attr is unavailable
            edge_attr = batch.edge_weight
        return super().forward(batch.x, batch.edge_index, edge_attr=edge_attr)

    def forward(self, batch):
        batch.x = self._forward(batch)
        return batch


class GATConv(BaseConvMixin, pygnn.GATConv):

    _edge_usage = "edge_attr"


class GATv2Conv(BaseConvMixin, pygnn.GATv2Conv):

    _edge_usage = "edge_attr"


class GCNConv(BaseConvMixin, pygnn.GCNConv):

    _edge_usage = "edge_weight"


class GENConv(BaseConvMixin, pygnn.GENConv):

    _edge_usage = "edge_weight"


class PatchedGINConv(pygnn.GINConv):

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        *,
        num_layers: int = 2,
        hidden_channels: Optional[int] = None,
        eps: float = 0.0,
        train_eps: bool = False,
        **kwargs,
    ):
        hidden_channels = hidden_channels or out_channels
        mlp = pygnn.MLP(
            in_channels=in_channels,
            out_channels=out_channels,
            hidden_channels=hidden_channels,
            num_layers=num_layers,
            act="relu",
            norm="batch_norm",
        )
        super().__init__(mlp, eps=eps, train_eps=train_eps, **kwargs)


class GINConv(BaseConvMixin, PatchedGINConv):

    _edge_usage = "none"


class PatchedGINEConv(pygnn.GINEConv):

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        use_edge_feature: bool = True,
        *,
        num_layers: int = 2,
        hidden_channels: Optional[int] = None,
        eps: float = 0.0,
        train_eps: bool = False,
        edge_dim: Optional[int] = None,
        **kwargs,
    ):
        hidden_channels = hidden_channels or out_channels
        mlp = pygnn.MLP(
            in_channels=in_channels,
            out_channels=out_channels,
            hidden_channels=hidden_channels,
            num_layers=num_layers,
            act="relu",
            norm="batch_norm",
        )
        super().__init__(mlp, eps=eps, train_eps=train_eps, edge_dim=edge_dim, **kwargs)


class GINEConv(BaseConvMixin, PatchedGINEConv):

    _edge_usage = "edge_attr"


class SAGEConv(BaseConvMixin, pygnn.SAGEConv):

    _edge_usage = "none"


class GatedGCNConv(BaseConvMixin, pygnn.ResGatedGraphConv):

    _edge_usage = "none"


__all__ = [
    "GATConv",
    "GATv2Conv",
    "GCNConv",
    "GENConv",
    "GINConv",
    "GINEConv",
    "SAGEConv",
    "GatedGCNConv",
]


if __name__ == "__main__":
    import torch
    from torch_geometric.data import Data

    m = GATConv(2, 5)
    data = Data(
        x=torch.ones(4, 2),
        edge_index=torch.LongTensor([[1, 2, 3], [0, 0, 0]]),
    )

    print(f"{m=}")
    print(f"{data=}")
    print(f"{m.forward(data)=}")
