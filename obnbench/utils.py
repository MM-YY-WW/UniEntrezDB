from functools import partial
from typing import Optional

import numba
import numpy as np
import torch
from obnb import logger


def get_num_workers(num_workers: int) -> int:
    default_num_workers = numba.config.NUMBA_DEFAULT_NUM_THREADS
    if 0 < num_workers <= default_num_workers:
        return num_workers
    elif num_workers > default_num_workers:
        logger.warning(
            f"The specified {num_workers=} exceeds the maximum number of "
            f"available threads {default_num_workers}. "
            f"Setting to {default_num_workers} instead.",
        )
    elif num_workers < 0:
        raise ValueError(f"num_workers should be a non-negative, got {num_workers}")
    logger.info(f"Setting number of workers to: {default_num_workers}")
    return default_num_workers


def replace_random_split(
    dataset,
    train_ratio: float = 0.6,
    val_ratio: float = 0.2,
    min_num_pos: int = 1,
    max_tries: int = 100,
    random_state: Optional[int] = 42,
):
    """Replace the original study-bias-holdout splits with random splits.

    Args:
        dataset: OBNB dataset object.
        train_ratio: Training ratio.
        val_ratio: Validation ratio.
        min_num_pos: Minimum number of positive examples across tasks and
            splits requirement.
        max_tries: Number of tries for setting up random split so that the
            minimum number of positive requirement is met.
        random_state: Random state for controlling the randomness.

    """
    logger.info(f"Preparing random splits ({train_ratio=}, {val_ratio=})")

    labeled_mask = dataset.train_mask | dataset.val_mask | dataset.test_mask
    labeled_node_idx = torch.where(labeled_mask)[0].numpy()

    size = labeled_node_idx.size
    train_size = int(size * train_ratio)
    train_val_size = int(size * (train_ratio + val_ratio))

    rng = np.random.default_rng(random_state)
    for _ in range(max_tries):
        perm = rng.permutation(size)

        train_idx = labeled_node_idx[perm[:train_size]]
        val_idx = labeled_node_idx[perm[train_size:train_val_size]]
        test_idx = labeled_node_idx[perm[train_val_size:]]

        # Check each task has sufficient number of positives in all splits
        if min(
            dataset.y[train_idx].sum(0).min().item(),
            dataset.y[val_idx].sum(0).min().item(),
            dataset.y[test_idx].sum(0).min().item(),
        ) > min_num_pos:
            break

    else:
        raise ValueError(
            f"Failed to prepare split within {max_tries} tries "
            f"({min_num_pos=}, {random_state=})"
        )

    dataset.train_mask[:], dataset.val_mask[:], dataset.test_mask[:] = map(
        partial(idx_ary_to_mask_tensor, size=dataset.num_nodes),
        (train_idx, val_idx, test_idx),
    )


def idx_ary_to_mask_tensor(idx_ary: np.ndarray, size: int) -> torch.Tensor:
    mask = torch.zeros(size, 1, dtype=bool)
    mask[idx_ary] = True
    return mask
