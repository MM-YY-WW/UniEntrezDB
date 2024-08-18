import torch
import torchmetrics
from torchmetrics.utilities.data import dim_zero_cat


class AUROC(torchmetrics.classification.MultilabelAUROC):
    ...


class AP(torchmetrics.classification.MultilabelAveragePrecision):
    ...


class APOP(AP):

    def __init__(self, task="multilabel", *args, average="macro", **kwargs):
        if task != "multilabel":
            raise NotImplementedError(
                "AveragePrecisionOverPrior is ony implemented for "
                "multilabel task for now.",
            )
        self.__average = average
        super().__init__(task=task, *args, average="none", **kwargs)

    def compute(self) -> torch.Tensor:
        scores = super().compute()

        # XXX: Does not consider negative selection currently.
        target = dim_zero_cat(self.target)
        prior = target.sum(0).clamp(1) / target.shape[0]
        scores = torch.log2(scores / prior)

        if self.__average == "macro":
            return scores.mean()
        elif self.__average == "none":
            return scores
        else:
            raise ValueError(f"Unknown averaging option {self.__average!r}")


__all__ = [
    "AP",
    "APOP",
    "AUROC",
]
