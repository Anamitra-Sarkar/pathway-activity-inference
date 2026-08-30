import pandas as pd
import numpy as np
from data_pipeline.zscore import zscore_scores
from tests.conftest import make_injected_fixture


def test_zscore_basic_shape_and_values():
    expr_df, pathways, _ = make_injected_fixture()
    scores = zscore_scores(expr_df, pathways)
    assert scores.shape == (4, 3)
    # Pathway IFN_RESPONSE: group A should have positive z (upregulated), B negative
    assert scores.loc["A1", "IFN_RESPONSE"] > 0
    assert scores.loc["B1", "IFN_RESPONSE"] < 0
    # Mean of z-scores per gene across samples should be ~0, so pathway mean approx...


def test_zscore_missing_genes_handled():
    expr_df = pd.DataFrame({"G1": [1, 2], "G2": [3, 4]}, index=["S1", "S2"])
    pathways = {"PW_MISSING": {"genes": ["NOT_A_GENE"]}}
    scores = zscore_scores(expr_df, pathways)
    assert scores.loc["S1", "PW_MISSING"] == 0.0


def test_zscore_scale_by_sqrt_k():
    # Use fixture, check scale factor
    expr_df, pathways, _ = make_injected_fixture()
    s1 = zscore_scores(expr_df, pathways, scale_by_sqrt_k=False)
    s2 = zscore_scores(expr_df, pathways, scale_by_sqrt_k=True)
    # s2 = s1 * sqrt(k)
    # IFN pathway k=3
    import math
    assert abs(s2.loc["A1", "IFN_RESPONSE"] - s1.loc["A1", "IFN_RESPONSE"] * math.sqrt(3)) < 1e-9


def test_zscore_recovers_injected_signal():
    expr_df, pathways, labels = make_injected_fixture()
    scores = zscore_scores(expr_df, pathways)
    mean_A = scores.loc[["A1", "A2"], "IFN_RESPONSE"].mean()
    mean_B = scores.loc[["B1", "B2"], "IFN_RESPONSE"].mean()
    assert mean_A > mean_B + 1.0, f"zscore should recover injected signal: {mean_A} vs {mean_B}"
