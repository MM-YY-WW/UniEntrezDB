import lightning.pytorch as pl
from torch_geometric.data import Dataset
from torch_geometric.loader import DataLoader


class DataModule(pl.LightningModule):
    def __init__(
        self,
        dataset: Dataset,
        sampler: str = "full_batch",
        num_workers: int = 1,
        pin_memory: bool = True,
    ):
        self.node_ids = dataset.node_ids
        self.task_ids = dataset.task_ids
        self._setup_loader(dataset, sampler, num_workers, pin_memory)
        super().__init__()

    def _setup_loader(self, dataset, sampler, num_workers, pin_memory):
        pw = num_workers > 1
        if sampler == "full_batch":
            # NOTE: Technically we only need one dataloader, whcih is the same
            # regardless of train/val/test split. What matters is the ``split``
            # attribute within the batch that signals the usage of mask.
            self.dataloaders = [
                DataLoader(
                    dataset,
                    batch_size=1,
                    shuffle=False,
                    num_workers=num_workers,
                    pin_memory=pin_memory,
                    persistent_workers=pw,
                )
                for _ in range(3)
            ]
        else:
            raise ValueError(f"Unknown sampler {sampler!r}")

    def train_dataloader(self) -> DataLoader:
        return self.dataloaders[0]

    def val_dataloader(self) -> DataLoader:
        return self.dataloaders[1]

    def test_dataloader(self) -> DataLoader:
        return self.dataloaders[2]
