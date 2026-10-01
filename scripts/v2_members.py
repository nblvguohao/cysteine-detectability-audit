"""Level-1 members of the v2 stack: LightGBM rankers and a torch fusion net."""
from __future__ import annotations

import numpy as np

from v2_stack import group_sorted, objective_on

import lightgbm as lgb


LGB_BASE = {
    "verbosity": -1,
    "num_threads": 0,
    "max_bin": 127,
    "force_row_wise": True,
    "deterministic": True,
    "seed": 20260913,
    "feature_fraction_seed": 20260913,
    "bagging_seed": 20260913,
}

LGB_RANK_GRID = (
    {"num_leaves": 31, "min_data_in_leaf": 40, "learning_rate": 0.05, "feature_fraction": 0.30, "lambda_l2": 1.0},
    {"num_leaves": 63, "min_data_in_leaf": 20, "learning_rate": 0.05, "feature_fraction": 0.20, "lambda_l2": 5.0},
    {"num_leaves": 15, "min_data_in_leaf": 80, "learning_rate": 0.08, "feature_fraction": 0.40, "lambda_l2": 1.0},
)
LGB_BIN_GRID = (
    {"num_leaves": 31, "min_data_in_leaf": 40, "learning_rate": 0.05, "feature_fraction": 0.30, "lambda_l2": 1.0},
    {"num_leaves": 63, "min_data_in_leaf": 20, "learning_rate": 0.05, "feature_fraction": 0.20, "lambda_l2": 5.0},
    {"num_leaves": 15, "min_data_in_leaf": 80, "learning_rate": 0.08, "feature_fraction": 0.40, "lambda_l2": 1.0},
)

MAX_ROUNDS = 1200
BLOCK = 25
PATIENCE_BLOCKS = 8
MIN_DELTA = 1e-4


def _params(kind, config):
    params = dict(LGB_BASE)
    params.update(config)
    if kind == "rank":
        params.update({
            "objective": "lambdarank",
            "metric": "None",
            "lambdarank_truncation_level": 12,
            "label_gain": [0, 1],
        })
    else:
        params.update({"objective": "binary", "metric": "None"})
    return params


def _dataset(X, y, index, proteins, kind):
    if kind == "rank":
        ordered, sizes = group_sorted(index, proteins)
        data = lgb.Dataset(
            X[ordered], label=y[ordered], group=sizes, free_raw_data=False
        )
        return data, ordered
    data = lgb.Dataset(X[index], label=y[index], free_raw_data=False)
    return data, index


def lgb_select(kind, X, y, proteins, components, subtrain, validation, grid=None):
    """Stage 1: grid search plus boosting-round selection on the inner fold."""
    grid = grid or (LGB_RANK_GRID if kind == "rank" else LGB_BIN_GRID)
    X_valid = X[validation]
    best = None
    trace = []
    for config_index, config in enumerate(grid):
        params = _params(kind, config)
        data, _ = _dataset(X, y, subtrain, proteins, kind)
        booster = None
        accumulated = np.zeros(len(validation), dtype=np.float64)
        rounds_done = 0
        best_rounds, best_objective, stale = BLOCK, -np.inf, 0
        best_accumulated = accumulated.copy()
        while rounds_done < MAX_ROUNDS:
            booster = lgb.train(
                params,
                data,
                num_boost_round=BLOCK,
                init_model=booster,
                keep_training_booster=True,
            )
            accumulated += booster.predict(
                X_valid,
                start_iteration=rounds_done,
                num_iteration=BLOCK,
                raw_score=True,
            )
            rounds_done += BLOCK
            objective = objective_on(
                y, accumulated, proteins, components, validation
            )
            trace.append({
                "config_index": config_index,
                "rounds": rounds_done,
                "validation_objective": float(objective),
            })
            if objective > best_objective + MIN_DELTA:
                best_objective, best_rounds, stale = objective, rounds_done, 0
                best_accumulated = accumulated.copy()
            else:
                stale += 1
            if stale >= PATIENCE_BLOCKS:
                break
        record = {
            "config": config,
            "config_index": config_index,
            "rounds": int(best_rounds),
            "validation_objective": float(best_objective),
            # exact partial sum of the first ``best_rounds`` trees on the inner fold
            "validation_raw_scores": best_accumulated,
        }
        if best is None or record["validation_objective"] > best["validation_objective"]:
            best = record
    return {
        "kind": f"lgb_{kind}",
        "config": best["config"],
        "config_index": best["config_index"],
        "rounds": best["rounds"],
        "validation_objective": best["validation_objective"],
        "trace": trace,
    }, np.asarray(best["validation_raw_scores"], dtype=np.float64)


def lgb_refit(kind, X, y, proteins, outer_train, test, selection):
    data, _ = _dataset(X, y, outer_train, proteins, kind)
    booster = lgb.train(
        _params(kind, selection["config"]), data, num_boost_round=selection["rounds"]
    )
    return booster.predict(X[test], raw_score=True), booster


# --------------------------------------------------------------------------- #
# torch fusion network


FUSION_MAX_EPOCHS = 400
FUSION_EVAL_EVERY = 5
FUSION_PATIENCE = 10
FUSION_HIDDEN = 256
FUSION_DROPOUT = 0.2
FUSION_LR = 1e-3
FUSION_WEIGHT_DECAY = 1e-4
FUSION_LISTWISE_WEIGHT = 0.5


def _torch():
    import torch

    device = "mps" if torch.backends.mps.is_available() else "cpu"
    return torch, device


def _standardize(train_block, blocks):
    median = np.nanmedian(train_block, axis=0)
    median = np.where(np.isfinite(median), median, 0.0)
    filled = np.where(np.isfinite(train_block), train_block, median)
    mean = filled.mean(axis=0)
    scale = filled.std(axis=0)
    scale[scale < 1e-6] = 1.0
    out = []
    for block in blocks:
        block = np.where(np.isfinite(block), block, median)
        out.append(((block - mean) / scale).astype(np.float32))
    return out


class _Fusion:
    def __init__(self, widths, seed):
        torch, device = _torch()
        self.torch = torch
        self.device = device
        torch.manual_seed(seed)
        self.branches = torch.nn.ModuleList([
            torch.nn.Sequential(
                torch.nn.Linear(width, FUSION_HIDDEN),
                torch.nn.GELU(),
                torch.nn.Dropout(FUSION_DROPOUT),
            )
            for width in widths
        ])
        self.head = torch.nn.Sequential(
            torch.nn.Linear(FUSION_HIDDEN * len(widths), FUSION_HIDDEN),
            torch.nn.GELU(),
            torch.nn.Dropout(FUSION_DROPOUT),
            torch.nn.Linear(FUSION_HIDDEN, 1),
        )
        self.model = torch.nn.ModuleList([self.branches, self.head]).to(device)

    def forward(self, blocks):
        torch = self.torch
        hidden = [branch(block) for branch, block in zip(self.branches, blocks)]
        return self.head(torch.cat(hidden, dim=1)).squeeze(1)


def _padded_groups(group_index_np, n_groups):
    """Dense (n_groups, max_size) gather layout for MPS-friendly listwise ops."""
    order = np.argsort(group_index_np, kind="stable")
    sorted_groups = group_index_np[order]
    counts = np.bincount(sorted_groups, minlength=n_groups)
    max_size = int(counts.max())
    layout = np.zeros((n_groups, max_size), dtype=np.int64)
    mask = np.zeros((n_groups, max_size), dtype=bool)
    cursor = 0
    for group, count in enumerate(counts):
        layout[group, :count] = order[cursor:cursor + count]
        mask[group, :count] = True
        cursor += count
    assert cursor == len(group_index_np)
    return layout, mask


def _listwise_loss(torch, logits, labels, layout, mask, usable):
    """Softmax cross-entropy against the in-group annotation distribution."""
    grouped = logits[layout]
    grouped = torch.where(mask, grouped, torch.full_like(grouped, -1e30))
    log_probability = grouped - torch.logsumexp(grouped, dim=1, keepdim=True)
    grouped_labels = torch.where(labels[layout], mask, torch.zeros_like(mask)).to(
        logits.dtype
    )
    positives = grouped_labels.sum(dim=1, keepdim=True)
    target = grouped_labels / torch.clamp(positives, min=1.0)
    per_group = -(target * torch.where(mask, log_probability, torch.zeros_like(log_probability))).sum(dim=1)
    return per_group[usable].mean()


def fusion_run(
    X_tab,
    embeddings,
    centred,
    y,
    proteins,
    components,
    train_index,
    eval_index,
    seed,
    epochs=None,
):
    """Train the fusion net; select epochs when ``epochs`` is None."""
    torch, device = _torch()
    blocks_train = []
    blocks_eval = []
    for source in (X_tab, embeddings, centred):
        train_block, eval_block = _standardize(
            source[train_index].astype(np.float64),
            [source[train_index].astype(np.float64), source[eval_index].astype(np.float64)],
        )
        blocks_train.append(torch.from_numpy(train_block).to(device))
        blocks_eval.append(torch.from_numpy(eval_block).to(device))

    labels = torch.from_numpy(y[train_index].astype(np.float32)).to(device)
    unique, group_index_np = np.unique(proteins[train_index], return_inverse=True)
    n_groups = len(unique)
    layout_np, mask_np = _padded_groups(group_index_np.astype(np.int64), n_groups)
    layout = torch.from_numpy(layout_np).to(device)
    mask = torch.from_numpy(mask_np).to(device)
    label_bool = torch.from_numpy((y[train_index] > 0)).to(device)
    positives_per_group = torch.from_numpy(
        np.bincount(group_index_np, weights=y[train_index], minlength=n_groups).astype(
            np.float32
        )
    ).to(device)
    usable = positives_per_group > 0

    net = _Fusion([block.shape[1] for block in blocks_train], seed)
    optimizer = torch.optim.AdamW(
        net.model.parameters(), lr=FUSION_LR, weight_decay=FUSION_WEIGHT_DECAY
    )
    bce = torch.nn.BCEWithLogitsLoss()

    def predict():
        net.model.eval()
        with torch.no_grad():
            return torch.sigmoid(net.forward(blocks_eval)).float().cpu().numpy()

    trace = []
    best_epoch, best_objective, stale = FUSION_EVAL_EVERY, -np.inf, 0
    total = epochs if epochs is not None else FUSION_MAX_EPOCHS
    for epoch in range(1, total + 1):
        net.model.train()
        optimizer.zero_grad(set_to_none=True)
        logits = net.forward(blocks_train)
        loss = FUSION_LISTWISE_WEIGHT * _listwise_loss(
            torch, logits, label_bool, layout, mask, usable
        ) + (1.0 - FUSION_LISTWISE_WEIGHT) * bce(logits, labels)
        loss.backward()
        optimizer.step()
        if epochs is not None or epoch % FUSION_EVAL_EVERY != 0:
            continue
        scores = predict()
        objective = objective_on(y, scores, proteins, components, eval_index)
        trace.append({"epoch": epoch, "loss": float(loss.detach().cpu()), "validation_objective": float(objective)})
        if objective > best_objective + MIN_DELTA:
            best_objective, best_epoch, stale = objective, epoch, 0
        else:
            stale += 1
        if stale >= FUSION_PATIENCE:
            break
    scores = predict()
    return {
        "kind": "fusion_mlp",
        "seed": seed,
        "epochs": int(epochs if epochs is not None else best_epoch),
        "validation_objective": float(best_objective) if epochs is None else None,
        "trace": trace,
        "device": device,
        "hyperparameters": {
            "hidden": FUSION_HIDDEN,
            "dropout": FUSION_DROPOUT,
            "learning_rate": FUSION_LR,
            "weight_decay": FUSION_WEIGHT_DECAY,
            "listwise_weight": FUSION_LISTWISE_WEIGHT,
            "max_epochs": FUSION_MAX_EPOCHS,
            "evaluate_every": FUSION_EVAL_EVERY,
            "patience_evaluations": FUSION_PATIENCE,
        },
    }, scores


def fusion_select_and_refit(
    X_tab, embeddings, centred, y, proteins, components,
    subtrain, validation, outer_train, test, seeds,
):
    """Stage 1 epoch selection per seed, then stage 2 refit on outer-train."""
    validation_scores = []
    test_scores = []
    records = []
    for seed in seeds:
        selection, valid_scores = fusion_run(
            X_tab, embeddings, centred, y, proteins, components,
            subtrain, validation, seed,
        )
        refit, fold_scores = fusion_run(
            X_tab, embeddings, centred, y, proteins, components,
            outer_train, test, seed, epochs=selection["epochs"],
        )
        validation_scores.append(valid_scores)
        test_scores.append(fold_scores)
        records.append(selection)
    return (
        records,
        np.mean(np.stack(validation_scores), axis=0),
        np.mean(np.stack(test_scores), axis=0),
    )
