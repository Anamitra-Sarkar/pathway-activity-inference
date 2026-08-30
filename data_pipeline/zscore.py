"""Combined z-score pathway scoring (Lee et al. 2008).

Reference:
    Lee E et al. The activity of pathway-based inferred networks.
    PLoS Comput Biol. 2008;4(11):e1000218. doi:10.1371/journal.pcbi.1000218
    - Introduced "combined z-score" for pathway activity:
      pathway score = average of gene-wise z-scores across genes in pathway,
      optionally scaled by sqrt(k).

    Also related:
    Levine et al. Pathway-based compressed sensing, similar z-score aggregates.

Method (this implementation):
  Given expression matrix E (samples x genes):
    1. Compute per-gene mean mu_g and std sigma_g across samples.
       (ddof=1; genes with zero variance set sigma=1 to avoid div/0).
    2. Standardize: z_{s,g} = (E_{s,g} - mu_g) / sigma_g
    3. For each pathway G with k_effective genes present:
         score_{s,G} = mean_{g in G} z_{s,g}
       (equivalent to (sum z)/k; Lee's scaled version (sum z)/sqrt(k) = mean*z * sqrt(k)
        is available via `scale_by_sqrt_k=True` for direct comparability.)

    This is a within-dataset (across-sample) standardization; it captures
    context-specific relative activity. For single-sample without cohort, use
    within-sample z-score variant (not implemented here, documented).

Returns: DataFrame samples x pathways.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def zscore_scores(
    expr_df: pd.DataFrame,
    gene_sets: dict,
    scale_by_sqrt_k: bool = False,
) -> pd.DataFrame:
    """Compute combined z-score pathway activity.

    Args:
        expr_df: DataFrame samples x genes.
        gene_sets: dict pathway_name -> {\"genes\": [...]} or list.
        scale_by_sqrt_k: if True, score = sum(z)/sqrt(k) (Lee's original scaling).
                         if False, score = mean(z).

    Returns:
        DataFrame samples x pathways.
    """
    if expr_df.empty:
        raise ValueError("Expression matrix is empty")

    normalized: dict[str, list] = {}
    for name, entry in gene_sets.items():
        if isinstance(entry, dict):
            genes = list(entry.get("genes", []))
        else:
            genes = list(entry)
        normalized[name] = genes

    # Compute gene-wise z
    means = expr_df.mean(axis=0)
    stds = expr_df.std(axis=0, ddof=1)
    # avoid zero sd
    stds = stds.replace(0, 1.0).fillna(1.0)
    z = (expr_df - means) / stds

    result = pd.DataFrame(index=expr_df.index, columns=normalized.keys(), dtype=float)

    for pathway, genes in normalized.items():
        # intersect with available columns
        present = [g for g in genes if g in z.columns]
        k = len(present)
        if k == 0:
            result[pathway] = 0.0
            continue
        sub = z[present]
        mean_z = sub.mean(axis=1)
        if scale_by_sqrt_k:
            mean_z = mean_z * np.sqrt(k)  # sum/sqrt(k) = mean * sqrt(k)
        result[pathway] = mean_z

    return result
