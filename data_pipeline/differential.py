"""Differential pathway activity: Wilcoxon + Benjamini-Hochberg FDR.

- Wilcoxon rank-sum (Mann-Whitney U, two-sided) per pathway comparing
  pathway scores between two sample groups.
- BH-FDR correction across pathways (independent or positively dependent tests).

References:
- Mann & Whitney 1947; Wilcoxon 1945 for rank-sum.
- Benjamini & Hochberg 1995 JRSS-B for FDR control.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def benjamini_hochberg(pvals: np.ndarray | list) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Args:
        pvals: array of p-values (0..1)

    Returns:
        array of q-values (FDR-adjusted), monotonic non-decreasing with sorted p order
        but returned in original order, clipped to 1.0.
    """
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    # argsort
    order = np.argsort(p)
    sorted_p = p[order]
    q = np.empty(n, dtype=float)
    # BH: q_i = p_i * n / rank, then cumulative min from largest to smallest
    # rank is 1-indexed
    prev = 1.0
    # iterate reversed to enforce monotonicity (step-up)
    for i in reversed(range(n)):
        rank = i + 1
        val = sorted_p[i] * n / rank
        if val > 1.0:
            val = 1.0
        # ensure monotonic non-increasing when going backwards (i.e., q sorted ascending monotonic)
        if val > prev:
            val = prev
        prev = val
        q[i] = val
    # reorder to original order
    inv_order = np.empty(n, dtype=int)
    inv_order[order] = np.arange(n)
    return q[inv_order]


def differential_analysis(
    scores_df: pd.DataFrame,
    labels: pd.Series | dict | list,
    group_a: str | None = None,
    group_b: str | None = None,
) -> pd.DataFrame:
    """Differential pathway activity between two groups.

    Args:
        scores_df: DataFrame samples x pathways (index samples)
        labels: mapping sample -> group label. Can be Series indexed by sample,
                dict, or list aligned to scores_df.index order.
        group_a, group_b: explicit group names. If None, inferred as the two unique labels
                          (sorted). If more than 2 groups present, must specify.

    Returns:
        DataFrame indexed by pathway with columns:
            mean_A, mean_B, delta (B-A), log2FC_approx, U_stat, p_value, q_value
        Sorted by q_value ascending then p_value.

    Raises:
        ValueError if groups invalid or sample mismatch.
    """
    # Normalize labels to Series
    if isinstance(labels, pd.Series):
        label_series = labels
    elif isinstance(labels, dict):
        label_series = pd.Series(labels)
    elif isinstance(labels, list):
        if len(labels) != len(scores_df.index):
            raise ValueError("labels list length must match scores_df rows")
        label_series = pd.Series(labels, index=scores_df.index)
    else:
        raise TypeError("labels must be Series, dict, or list")

    # Align
    # Ensure all samples in scores_df have label
    missing = set(scores_df.index) - set(label_series.index)
    if missing:
        raise ValueError(f"Missing labels for samples: {missing}")
    label_series = label_series.reindex(scores_df.index)

    unique = sorted(label_series.dropna().unique().tolist())
    if group_a is None or group_b is None:
        if len(unique) != 2:
            raise ValueError(f"Expected 2 groups, found {unique}; specify group_a/group_b")
        group_a, group_b = unique[0], unique[1]
    else:
        if group_a not in unique or group_b not in unique:
            raise ValueError(f"Specified groups not found in labels: {unique}")

    idx_a = label_series[label_series == group_a].index
    idx_b = label_series[label_series == group_b].index

    if len(idx_a) < 2 or len(idx_b) < 2:
        # mannwhitney requires at least 1 but warn if too small; allow but note
        pass
    if len(idx_a) == 0 or len(idx_b) == 0:
        raise ValueError("One group has zero samples")

    rows = []
    for pathway in scores_df.columns:
        vals_a = scores_df.loc[idx_a, pathway].dropna().values
        vals_b = scores_df.loc[idx_b, pathway].dropna().values
        mean_a = float(np.mean(vals_a)) if len(vals_a) else float("nan")
        mean_b = float(np.mean(vals_b)) if len(vals_b) else float("nan")
        delta = mean_b - mean_a
        # Mann-Whitney U two-sided
        try:
            # Use asymptotic if large, exact if small; scipy handles
            res = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
            u = float(res.statistic)
            p = float(res.pvalue)
        except Exception:
            u = float("nan")
            p = 1.0
        rows.append((pathway, mean_a, mean_b, delta, u, p))

    df = pd.DataFrame(rows, columns=["pathway", "mean_A", "mean_B", "delta", "U_stat", "p_value"]).set_index("pathway")
    # BH correction
    q = benjamini_hochberg(df["p_value"].values)
    df["q_value"] = q
    df["group_A"] = group_a
    df["group_B"] = group_b
    df = df.sort_values(["q_value", "p_value"])
    return df
