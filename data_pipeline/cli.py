"""CLI for pathway activity pipeline.

Example real-run (requires downloaded GMT and expression matrix):
  python -m data_pipeline.cli \
    --gmt-path data/h.all.v2023.2.Hs.symbols.gmt \
    --expr data/expression.csv \
    --method both \
    --outdir results/

Expression CSV expected: rows=samples, columns=genes, first column sample ID or index.
Or provide --expr-transpose if genes are rows.

For testing with synthetic fixtures, use tests.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .gmt_parser import parse_gmt
from .ssgsea import ssgsea_scores
from .zscore import zscore_scores
from .differential import differential_analysis
from .correlation import correlate_methods


def main():
    parser = argparse.ArgumentParser(description="Pathway activity inference")
    parser.add_argument("--gmt-path", required=True, help="Path to GMT file (MSigDB Hallmark or Reactome)")
    parser.add_argument("--expr", required=True, help="Path to expression CSV (samples x genes)")
    parser.add_argument("--method", choices=["ssgsea", "zscore", "both"], default="both")
    parser.add_argument("--alpha", type=float, default=0.25, help="ssGSEA alpha weighting (0,5]")
    parser.add_argument("--groups", help="Optional CSV with sample,group columns for differential analysis")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--expr-transpose", action="store_true", help="Transpose expression (genes x samples input)")
    args = parser.parse_args()

    if not 0 < args.alpha <= 5:
        parser.error("--alpha must be in (0, 5]")
    if not Path(args.gmt_path).exists():
        parser.error(f"GMT file not found: {args.gmt_path}")
    if not Path(args.expr).exists():
        parser.error(f"Expression CSV not found: {args.expr}")
    if args.groups and not Path(args.groups).exists():
        parser.error(f"Groups CSV not found: {args.groups}")

    gmt = parse_gmt(args.gmt_path)
    print(f"Loaded {len(gmt)} pathways from {args.gmt_path}")

    try:
        expr = pd.read_csv(args.expr, index_col=0)
    except Exception as e:
        parser.error(f"Failed to read expression CSV: {e}")
    if args.expr_transpose:
        expr = expr.T
    if expr.empty:
        parser.error("Expression matrix is empty")
    if expr.shape[0] == 0 or expr.shape[1] == 0:
        parser.error(f"Expression matrix has invalid shape {expr.shape}")
    # Check all columns numeric
    non_numeric = expr.select_dtypes(exclude="number").columns.tolist()
    # Try to coerce object columns that are numeric strings
    if non_numeric:
        for col in non_numeric:
            try:
                expr[col] = pd.to_numeric(expr[col], errors="raise")
            except Exception:
                parser.error(f"Expression column '{col}' contains non-numeric values")
    if expr.isna().all().all():
        parser.error("Expression matrix contains only NaN")
    print(f"Expression matrix: {expr.shape[0]} samples x {expr.shape[1]} genes")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.method in ("ssgsea", "both"):
        ss = ssgsea_scores(expr, gmt, alpha=args.alpha)
        ss.to_csv(outdir / "ssgsea_scores.csv")
        print(f"Wrote {outdir / 'ssgsea_scores.csv'}")

    if args.method in ("zscore", "both"):
        zs = zscore_scores(expr, gmt)
        zs.to_csv(outdir / "zscore_scores.csv")
        print(f"Wrote {outdir / 'zscore_scores.csv'}")

    if args.method == "both":
        ss = pd.read_csv(outdir / "ssgsea_scores.csv", index_col=0)
        zs = pd.read_csv(outdir / "zscore_scores.csv", index_col=0)
        corr = correlate_methods(ss, zs)
        corr.to_csv(outdir / "correlation.csv")
        print(f"Wrote {outdir / 'correlation.csv'}")
        print(corr.head())

    if args.groups:
        try:
            groups_df = pd.read_csv(args.groups, index_col=0)
        except Exception as e:
            parser.error(f"Failed to read groups CSV: {e}")
        if groups_df.empty:
            parser.error("Groups CSV is empty")
        # assume column 'group' or first column; handle whitespace in column names
        groups_df.columns = [c.strip() for c in groups_df.columns.astype(str)]
        if "group" in groups_df.columns:
            labels = groups_df["group"]
        else:
            labels = groups_df.iloc[:, 0]
        labels.index = labels.index.astype(str).str.strip()
        labels = labels.astype(str).str.strip()
        # drop empty labels
        labels = labels[labels != ""]
        if labels.nunique() < 2:
            parser.error(f"Groups must contain at least 2 distinct labels, found {labels.unique().tolist()}")
        # use ssGSEA or zscore whichever available
        scores_path = outdir / "ssgsea_scores.csv"
        if not scores_path.exists():
            scores_path = outdir / "zscore_scores.csv"
        scores = pd.read_csv(scores_path, index_col=0)
        diff = differential_analysis(scores, labels)
        diff.to_csv(outdir / "differential.csv")
        print(f"Wrote {outdir / 'differential.csv'}")
        print(diff.head())


if __name__ == "__main__":
    main()
