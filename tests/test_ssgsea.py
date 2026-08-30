import pandas as pd
import numpy as np
import pytest
from data_pipeline.ssgsea import ssgsea_score_single, ssgsea_scores
from tests.conftest import make_injected_fixture


def test_ssgsea_high_rank_scores_high():
    # Single sample: 10 genes, pathway top-ranked should score high
    expr = pd.Series({
        "G_top1": 10.0, "G_top2": 9.0, "G_top3": 8.0,
        "G_mid1": 5.0, "G_mid2": 4.5, "G_mid3": 4.0,
        "G_low1": 1.0, "G_low2": 0.5, "G_low3": 0.1, "G_low4": 0.0
    })
    # gene set with all top genes
    score_high = ssgsea_score_single(expr, {"G_top1", "G_top2", "G_top3"})
    # gene set with random mid
    score_random = ssgsea_score_single(expr, {"G_mid1", "G_mid2"})
    # gene set with all bottom genes
    score_low = ssgsea_score_single(expr, {"G_low1", "G_low2", "G_low3"})

    assert score_high > score_random, f"high {score_high} should > random {score_random}"
    assert score_random > score_low or abs(score_random) < abs(score_high), "random near zero, high vs low separation"
    # high should be positive, low negative
    assert score_high > 0, "top-ranked gene set should give positive ES"
    assert score_low < 0, "bottom-ranked gene set should give negative ES"
    # random near zero (within 2)
    assert abs(score_random) < abs(score_high), "random score magnitude smaller than enriched"


def test_ssgsea_empty_and_full_set_returns_zero():
    expr = pd.Series({"G1": 1.0, "G2": 2.0})
    assert ssgsea_score_single(expr, set()) == 0.0
    assert ssgsea_score_single(expr, {"G1", "G2"}) == 0.0  # k==N
    assert ssgsea_score_single(expr, {"NOT_PRESENT"}) == 0.0


def test_ssgsea_scores_dataframe_shape():
    expr_df, pathways, _ = make_injected_fixture()
    scores = ssgsea_scores(expr_df, pathways)
    assert scores.shape == (4, 3)
    assert list(scores.index) == ["A1", "A2", "B1", "B2"]
    assert set(scores.columns) == {"IFN_RESPONSE", "RANDOM_SET", "P53_PATHWAY"}


def test_ssgsea_recovers_injected_signal():
    """Honest synthetic verification: injected IFN genes up in group A should be recovered as high ssGSEA in A vs B.
    This is not a biological finding; it validates correctness."""
    expr_df, pathways, labels = make_injected_fixture()
    scores = ssgsea_scores(expr_df, pathways)
    # Pathway IFN_RESPONSE should have higher scores in A than B
    mean_A = scores.loc[["A1", "A2"], "IFN_RESPONSE"].mean()
    mean_B = scores.loc[["B1", "B2"], "IFN_RESPONSE"].mean()
    assert mean_A > mean_B + 1.0, f"IFN_RESPONSE mean_A {mean_A} should >> mean_B {mean_B} (injected signal)"
    # Random set should not show strong differential
    mean_A_rand = scores.loc[["A1", "A2"], "RANDOM_SET"].mean()
    mean_B_rand = scores.loc[["B1", "B2"], "RANDOM_SET"].mean()
    assert abs(mean_A_rand - mean_B_rand) < abs(mean_A - mean_B), "random set differential smaller than injected"


def test_ssgsea_alpha_sensitivity():
    expr = pd.Series({"G1": 10, "G2": 9, "G3": 8, "G4": 5, "G5": 1})
    s1 = ssgsea_score_single(expr, {"G1", "G2"}, alpha=0.25)
    s2 = ssgsea_score_single(expr, {"G1", "G2"}, alpha=1.0)
    # Different alpha should give different scores but same sign (k=2 so weighting matters)
    assert s1 > 0 and s2 > 0
    assert abs(s1 - s2) > 1e-6
