# RF_benchmark.py
# Random-forest benchmark implementing Section 2.6.2.
# Python 3.11; numpy/pandas/scipy/scikit-learn.
#
# Only the 92 samples common to all three enzymes are used.
# A separate 69/23 Kennard-Stone partition is generated from VIS-NIR.
# RF tuning is strictly calibration-only.

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple, List, Any

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

SEED = 125
N_TREES = 500

MIR_RANGE = {
    "BG": (853.6, 3954.7),
    "ACP": (853.6, 3928.6),
    "ARS": (853.6, 3954.7),
}

RF_COMMON = {
    "BG":  None,
    "ACP": None,
    "ARS": None,
}

def _numeric_axis(columns: List[Any]) -> np.ndarray:
    vals = []
    for c in columns:
        try:
            vals.append(float(str(c)))
        except ValueError:
            vals.append(np.nan)
    return np.asarray(vals, dtype=float)

def _spectral_columns(df: pd.DataFrame, lo: float, hi: float) -> Tuple[List[str], np.ndarray]:
    axis = _numeric_axis(list(df.columns))
    mask = np.isfinite(axis) & (axis >= min(lo, hi)) & (axis <= max(lo, hi))
    if not mask.any():
        raise ValueError(f"No spectral columns found in range {lo}-{hi}.")
    idx = np.where(mask)[0]
    order = np.argsort(axis[idx])
    idx = idx[order]
    return [str(df.columns[i]) for i in idx], axis[idx][order]

def reflectance_to_absorbance(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float).copy()
    X[X <= 0] = 1e-6
    return -np.log10(X)

def savgol_rows(X: np.ndarray, window: int, poly: int, deriv: int = 0) -> np.ndarray:
    if window % 2 == 0 or window <= poly:
        raise ValueError("Savitzky-Golay window must be odd and > polynomial order.")
    return savgol_filter(X, window_length=window, polyorder=poly, deriv=deriv, axis=1)

def rf_preprocess_vis(X: np.ndarray) -> np.ndarray:
    X = reflectance_to_absorbance(X)
    X = savgol_rows(X, 5, 2, 0)
    X = savgol_rows(X, 11, 2, 1)
    return X

def rf_preprocess_mir(X: np.ndarray) -> np.ndarray:
    # MIR is acquired in absorbance mode. The common RF pipeline is
    # smoothing followed by first derivative.
    X = savgol_rows(X, 5, 2, 0)
    X = savgol_rows(X, 11, 2, 1)
    return X

def _mahalanobis_matrix(Z: np.ndarray) -> np.ndarray:
    cov = np.cov(Z, rowvar=False)
    cov = np.atleast_2d(cov)
    inv = np.linalg.pinv(cov + np.eye(cov.shape[0]) * 1e-10)
    n = Z.shape[0]
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        d = Z - Z[i]
        D[i] = np.sqrt(np.einsum("ij,jk,ik->i", d, inv, d))
    return D

def kennard_stone(X: np.ndarray, calibration_fraction: float = 0.75,
                  variance_explained: float = 0.95) -> Tuple[np.ndarray, np.ndarray, int]:
    """Kennard-Stone using the minimum number of PCs explaining >=95% variance."""
    pca = PCA()
    Z0 = pca.fit_transform(X)
    cum = np.cumsum(pca.explained_variance_ratio_)
    n_pc = int(np.searchsorted(cum, variance_explained) + 1)
    Z = Z0[:, :n_pc]
    D = _mahalanobis_matrix(Z)
    i, j = np.unravel_index(np.argmax(D), D.shape)
    selected = [int(i)]
    if j != i:
        selected.append(int(j))
    k = int(round(X.shape[0] * calibration_fraction))
    k = min(max(k, 2), X.shape[0] - 1)
    while len(selected) < k:
        remain = [q for q in range(X.shape[0]) if q not in selected]
        nearest = np.min(D[np.ix_(remain, selected)], axis=1)
        selected.append(remain[int(np.argmax(nearest))])
    cal = np.array(sorted(selected), dtype=int)
    val = np.array(sorted(set(range(X.shape[0])) - set(cal)), dtype=int)
    return cal, val, n_pc

def _max_features_values(p: int):
    return [
        max(1, int(round(math.sqrt(p)))),
        max(1, int(round(0.10 * p))),
        max(1, int(round(0.33 * p))),
    ]

def tune_and_evaluate(X: np.ndarray, y: np.ndarray,
                      cal_idx: np.ndarray, val_idx: np.ndarray):
    p = X.shape[1]
    max_features = _max_features_values(p)

    base = RandomForestRegressor(
        n_estimators=N_TREES,
        random_state=SEED,
        n_jobs=-1,
    )

    grid = {
        "max_features": max_features,
        "max_depth": [None, 10],
        "min_samples_leaf": [1, 3, 5],
    }

    cv = KFold(n_splits=5, shuffle=True, random_state=SEED)
    search = GridSearchCV(
        estimator=base,
        param_grid=grid,
        scoring="neg_root_mean_squared_error",
        cv=cv,
        refit=True,
        n_jobs=-1,
        return_train_score=True,
    )
    search.fit(X[cal_idx], y[cal_idx])

    pred_cal = search.predict(X[cal_idx])
    pred_val = search.predict(X[val_idx])

    return {
        "search": search,
        "best_params": search.best_params_,
        "best_cv_rmse": -float(search.best_score_),
        "calibration_prediction": pred_cal,
        "validation_prediction": pred_val,
    }

def metrics(obs: np.ndarray, pred: np.ndarray) -> Dict[str, float]:
    rmse = float(np.sqrt(mean_squared_error(obs, pred)))
    mae = float(mean_absolute_error(obs, pred))
    r2 = float(r2_score(obs, pred))
    rpd = float(np.std(obs, ddof=1) / rmse)
    q = float(np.percentile(obs, 75) - np.percentile(obs, 25))
    rpiq = q / rmse
    return {"R2": r2, "RMSE": rmse, "MAE": mae, "RPD": rpd, "RPIQ": rpiq}

def run_enzyme(vis: pd.DataFrame, mir: pd.DataFrame, enzyme: str) -> Dict[str, Any]:
    if enzyme not in MIR_RANGE:
        raise ValueError(f"enzyme must be one of {list(MIR_RANGE)}")

    if len(vis) != 92 or len(mir) != 92:
        raise ValueError(
            "The RF benchmark requires exactly the 92 samples common to BG, ACP and ARS."
        )
    if not np.array_equal(vis["sampleID"].astype(str).to_numpy(),
                          mir["sampleID"].astype(str).to_numpy()):
        raise ValueError("VIS and MIR sampleID order must be identical.")

    vis_cols, vis_axis = _spectral_columns(vis, 394, 1674)
    mir_cols, mir_axis = _spectral_columns(mir, *MIR_RANGE[enzyme])

    Xvis_raw = vis[vis_cols].to_numpy(dtype=float)
    Xmir_raw = mir[mir_cols].to_numpy(dtype=float)

    Xvis = rf_preprocess_vis(Xvis_raw)
    Xmir = rf_preprocess_mir(Xmir_raw)
    y = vis[enzyme].to_numpy(dtype=float)

    cal_idx, val_idx, n_pc = kennard_stone(
        Xvis, calibration_fraction=0.75, variance_explained=0.95
    )

    results = {}
    for name, X in {
        "RF_VIS_NIR": Xvis,
        "RF_MIR": Xmir,
        "RF_Direct_Fusion": np.column_stack([Xvis, Xmir]),
    }.items():
        tuned = tune_and_evaluate(X, y, cal_idx, val_idx)
        results[name] = {
            "best_params": tuned["best_params"],
            "best_cv_rmse": tuned["best_cv_rmse"],
            "calibration_metrics": metrics(y[cal_idx], tuned["calibration_prediction"]),
            "validation_metrics": metrics(y[val_idx], tuned["validation_prediction"]),
            "validation_prediction": tuned["validation_prediction"],
            "calibration_prediction": tuned["calibration_prediction"],
        }

    performance = pd.DataFrame({
        name: value["validation_metrics"]
        for name, value in results.items()
    }).T

    return {
        "enzyme": enzyme,
        "calibration_index": cal_idx,
        "validation_index": val_idx,
        "n_pc_for_95pct": n_pc,
        "performance": performance,
        "best_params": {k: v["best_params"] for k, v in results.items()},
        "models": results,
        "vis_axis": vis_axis,
        "mir_axis": mir_axis,
    }

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--vis", required=True)
    parser.add_argument("--mir", required=True)
    parser.add_argument("--enzyme", required=True, choices=["BG", "ACP", "ARS"])
    args = parser.parse_args()

    vis = pd.read_csv(args.vis)
    mir = pd.read_csv(args.mir)
    out = run_enzyme(vis, mir, args.enzyme)
    print(out["performance"].round(4))
    print("\nSelected hyperparameters:")
    for name, params in out["best_params"].items():
        print(name, params)
