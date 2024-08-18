import torch
from torch import nn, Tensor


class RawFeatNorm(nn.Module):

    def __init__(
        self,
        num_features: int,
        eps: float = 1e-5,
        affine: bool = True,
    ):
        super().__init__()
        self.num_features = num_features
        self.eps = eps
        self.affine = affine

        if self.affine:
            self.weight = nn.Parameter(torch.empty(num_features))
            self.bias = nn.Parameter(torch.empty(num_features))

        self.reset_parameters()

    def __repr__(self) -> str:
        paramstr = f"{self.num_features}, eps={self.eps}, affine={self.affine}"
        return f"{self.__class__.__name__}({paramstr})"

    def reset_parameters(self):
        if self.affine:
            nn.init.ones_(self.weight)
            nn.init.zeros_(self.bias)

    def forward(self, x: Tensor) -> Tensor:
        x = (x - x.mean(0)) / (x.std(0) + self.eps)

        if self.affine:
            x = x * self.weight + self.bias

        return x
