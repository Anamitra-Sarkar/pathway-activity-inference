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
    parser.add_argument("--alpha", type=float, default=0.25, help="ssGSEA alpha weighting")
    parser.add_argument("--groups", help="Optional CSV with sample,group columns for differential analysis")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--expr-transpose", action="store_true", help="Transpose expression (genes x samples input)")
    args = parser.parse_args()

    gmt = parse_gmt(args.gmt_path)
    print(f"Loaded {len(gmt)} pathways from {args.gmt_path}")

    expr = pd.read_csv(args.expr, index_col=0)
    if args.expr_transpose:
        expr = expr.T
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
        groups_df = pd.read_csv(args.groups, index_col=0)
        # assume column 'group' or first column
        if "group" in groups_df.columns:
            labels = groups_df["group"]
        else:
            labels = groups_df.iloc[:, 0]
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
