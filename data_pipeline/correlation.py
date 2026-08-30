"""Correlation / agreement analysis between two pathway scoring methods."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def correlate_methods(
    scores_a: pd.DataFrame,
    scores_b: pd.DataFrame,
    method: str = "both",
) -> pd.DataFrame:
    """Compute per-pathway correlation between two scoring methods.

    Both DataFrames must have same index (samples) and overlapping columns (pathways).
    Typically scores_a = ssGSEA, scores_b = z-score.

    Args:
        scores_a, scores_b: DataFrames samples x pathways
        method: \"pearson\" or \"spearman\" or \"both\"

    Returns:
        DataFrame indexed by pathway with columns:
            pearson_r, pearson_p, spearman_r, spearman_p, n_samples
        If method != both, still returns all but caller can filter.
    """
    # Align
    common_samples = scores_a.index.intersection(scores_b.index)
    if len(common_samples) == 0:
        raise ValueError("No overlapping samples")
    common_pathways = scores_a.columns.intersection(scores_b.columns)
    if len(common_pathways) == 0:
        raise ValueError("No overlapping pathways")

    rows = []
    for pw in common_pathways:
        a = scores_a.loc[common_samples, pw].values
        b = scores_b.loc[common_samples, pw].values
        # drop nan pairs
        mask = ~(np.isnan(a) | np.isnan(b))
        a_clean = a[mask]
        b_clean = b[mask]
        n = len(a_clean)
        if n < 3:
            rows.append((pw, float("nan"), float("nan"), float("nan"), float("nan"), n))
            continue
        # Handle constant input (pearson undefined)
        if np.std(a_clean) == 0 or np.std(b_clean) == 0:
            rows.append((pw, float("nan"), float("nan"), float("nan"), float("nan"), n))
            continue
        try:
            pr, pp = pearsonr(a_clean, b_clean)
        except Exception:
            pr, pp = float("nan"), float("nan")
        try:
            sr, sp = spearmanr(a_clean, b_clean)
        except Exception:
            sr, sp = float("nan"), float("nan")
        rows.append((pw, float(pr), float(pp), float(sr), float(sp), n))

    df = pd.DataFrame(
        rows, columns=["pathway", "pearson_r", "pearson_p", "spearman_r", "spearman_p", "n_samples"]
    ).set_index("pathway")

    if method == "pearson":
        return df[["pearson_r", "pearson_p", "n_samples"]]
    if method == "spearman":
        return df[["spearman_r", "spearman_p", "n_samples"]]
    return df
