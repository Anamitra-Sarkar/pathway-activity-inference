"""ssGSEA implementation (Barbie et al. 2009).

Reference:
    Barbie DA et al. Systematic RNA interference reveals that oncogenic KRAS-driven
    cancers require TBK1. Nature. 2009;462(7269):108-112. doi:10.1038/nature08460
    - Introduced single-sample GSEA (ssGSEA).

    GSEA original:
    Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach
    for interpreting genome-wide expression profiles. PNAS 2005.

    Clarification / weighting details:
    Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for
    microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.

This module implements ssGSEA as a rank-weighted Kolmogorov-Smirnov running-sum
statistic per sample per gene set, with the integrated ES (area) variant.

Algorithm (per sample, per gene set G, size k, total N genes):
  1. Rank genes in sample by expression descending (highest first).
  2. Define alpha = 0.25 (Barbie default, as in GenePattern ssGSEA).
  3. Compute weight for gene at sorted position i (0-indexed):
        w_i = (N - i) ** alpha   if gene in G, else 0
     N_R = sum_{g in G} w_{pos(g)}   (normalizer)
  4. Walk the ranked list i=0..N-1 iteratively:
        P_hit(i)  cumulates w_i/N_R for genes in G
        P_miss(i) cumulates 1/(N-k) for genes not in G
        ECDF difference at i: D(i) = P_hit(i) - P_miss(i)
  5. ssGSEA enrichment score = sum_{i=0}^{N-1} D(i)
     (integrated difference; some formulations return max(D) but we use sum
      following GSVA's description of ssGSEA).

Edge handling:
  - Genes in gene set but not in expression matrix are ignored (k_effective).
  - If k_effective == 0 or k_effective == N: score = 0 (undefined enrichment).
  - Expression NaNs are dropped per sample.

Returns: DataFrame samples x pathways (samples as index).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def ssgsea_score_single(
    expr: pd.Series,
    gene_set: set | list,
    alpha: float = 0.25,
) -> float:
    """Compute ssGSEA score for a single sample.

    Args:
        expr: Series indexed by gene symbol, values = expression.
        gene_set: iterable of gene symbols.
        alpha: weighting exponent (0.25 Barbie default).

    Returns:
        float enrichment score (integrated ES).
    """
    # Drop NA and ensure numeric
    expr = expr.dropna()
    if expr.empty:
        return 0.0
    gene_set = set(gene_set)
    # intersect
    present = set(expr.index) & gene_set
    k = len(present)
    N = len(expr)
    if k == 0 or k == N:
        return 0.0

    # Sort descending by expression
    sorted_genes = expr.sort_values(ascending=False).index.tolist()
    # Map position to weight: (N - i) ** alpha
    # Precompute NR
    # Build dict gene->pos for fast lookup? But we loop in order anyway
    NR = 0.0
    for i, g in enumerate(sorted_genes):
        if g in present:
            NR += (N - i) ** alpha
    if NR == 0:
        return 0.0

    cum_hit = 0.0
    cum_miss = 0.0
    es_sum = 0.0
    miss_increment = 1.0 / (N - k)

    for i, g in enumerate(sorted_genes):
        if g in present:
            cum_hit += (N - i) ** alpha / NR
        else:
            cum_miss += miss_increment
        es_sum += cum_hit - cum_miss

    return float(es_sum)


def ssgsea_scores(
    expr_df: pd.DataFrame,
    gene_sets: dict,
    alpha: float = 0.25,
) -> pd.DataFrame:
    """Compute ssGSEA scores for all samples x gene sets.

    Args:
        expr_df: DataFrame samples x genes (index=samples, columns=genes) OR genes x samples?
                 This function expects samples as rows, genes as columns.
                 If you have genes as rows, transpose before calling.
        gene_sets: dict pathway_name -> {\"genes\": [...]} OR pathway_name -> list of genes.
        alpha: weighting exponent, must be >0 and <=5.

    Returns:
        DataFrame samples x pathways (index same as expr_df.index, columns pathway names).
    """
    if expr_df.empty:
        raise ValueError("Expression matrix is empty")
    if not isinstance(alpha, (int, float)) or not 0 < alpha <= 5:
        raise ValueError(f"alpha must be in (0, 5], got {alpha}")
    if not isinstance(gene_sets, dict) or len(gene_sets) == 0:
        raise ValueError("gene_sets must be non-empty dict")
    if expr_df.isna().all().all():
        raise ValueError("Expression matrix contains only NaN")

    # Normalize gene_sets to dict name -> set
    normalized: dict[str, set] = {}
    for name, entry in gene_sets.items():
        if isinstance(entry, dict):
            genes = entry.get("genes", [])
        else:
            genes = entry
        normalized[name] = set(genes)

    # Prepare output
    result = pd.DataFrame(index=expr_df.index, columns=normalized.keys(), dtype=float)

    for sample in expr_df.index:
        row = expr_df.loc[sample]
        # row is Series genes -> expression
        row.index = row.index.astype(str)
        for pathway, gene_set in normalized.items():
            result.at[sample, pathway] = ssgsea_score_single(row, gene_set, alpha=alpha)

    return result
