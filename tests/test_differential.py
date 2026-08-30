import pandas as pd
import numpy as np
from data_pipeline.differential import benjamini_hochberg, differential_analysis
from data_pipeline.ssgsea import ssgsea_scores
from data_pipeline.zscore import zscore_scores
from tests.conftest import make_injected_fixture


def test_bh_monotonicity_and_bounds():
    pvals = [0.01, 0.04, 0.03, 0.2]
    q = benjamini_hochberg(pvals)
    # q should be monotonic with sorted p
    order = np.argsort(pvals)
    sorted_q = q[order]
    for i in range(1, len(sorted_q)):
        assert sorted_q[i] >= sorted_q[i-1] - 1e-12, "q should be non-decreasing with p"
    assert all(0 <= x <= 1 for x in q)
    # p=0.2 with n=4 -> q <=1
    assert q[3] <= 1.0


def test_bh_empty():
    assert len(benjamini_hochberg([])) == 0


def test_differential_recovers_injected_signal_ssgsea():
    expr_df, pathways, labels = make_injected_fixture()
    scores = ssgsea_scores(expr_df, pathways)
    diff = differential_analysis(scores, labels)
    # IFN_RESPONSE should be most significant (smallest q) among 3; with n=2+2 power is limited, q may be 0.33-1.0
    top = diff.index[0]
    assert top == "IFN_RESPONSE", f"Expected IFN_RESPONSE top, got {top}: {diff}"
    # q may be up to 1 with small n, just ensure it is the smallest
    assert diff.loc["IFN_RESPONSE", "q_value"] <= diff.loc["RANDOM_SET", "q_value"]
    # delta should be B - A negative (A high)
    assert diff.loc["IFN_RESPONSE", "delta"] < -1.0


def test_differential_recovers_injected_signal_zscore():
    expr_df, pathways, labels = make_injected_fixture()
    scores = zscore_scores(expr_df, pathways)
    diff = differential_analysis(scores, labels)
    assert diff.index[0] == "IFN_RESPONSE"


def test_differential_explicit_groups():
    scores = pd.DataFrame({"PW1": [1, 2, 10, 11], "PW2": [5,5,5,5]}, index=["S1","S2","S3","S4"])
    labels = pd.Series(["A","A","B","B"], index=["S1","S2","S3","S4"])
    diff = differential_analysis(scores, labels, group_a="A", group_b="B")
    # PW1 should be significant, PW2 not
    assert diff.loc["PW1", "p_value"] < diff.loc["PW2", "p_value"]
    # Also test group swapping gives opposite delta sign
    diff_swapped = differential_analysis(scores, labels, group_a="B", group_b="A")
    assert diff_swapped.loc["PW1", "delta"] == - diff.loc["PW1", "delta"]
