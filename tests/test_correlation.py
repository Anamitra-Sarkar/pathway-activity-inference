import pandas as pd
import numpy as np
from data_pipeline.ssgsea import ssgsea_scores
from data_pipeline.zscore import zscore_scores
from data_pipeline.correlation import correlate_methods
from tests.conftest import make_injected_fixture


def test_correlation_agreement_on_injected_data():
    """Honest quantitative comparison: ssGSEA and z-score should positively correlate
    on the injected signal pathway (both capture same upregulation), but not perfect."""
    expr_df, pathways, _ = make_injected_fixture()
    # Use more samples for correlation stability: replicate to 8 samples
    # Extend to 8 by duplicating with noise? Simpler: generate 12 samples
    import pandas as pd
    rng = np.random.default_rng(123)
    genes = expr_df.columns.tolist()
    big_data = {}
    labels = {}
    for i in range(6):
        s = f"A{i}"
        big_data[s] = {g: (9.0 if g in ("IFIT1","MX1","ISG15") else 5.0) + rng.normal(0,0.4) for g in genes}
        labels[s]="A"
    for i in range(6):
        s = f"B{i}"
        big_data[s] = {g: (3.0 if g in ("IFIT1","MX1","ISG15") else 5.0) + rng.normal(0,0.4) for g in genes}
        labels[s]="B"
    big_df = pd.DataFrame.from_dict(big_data, orient="index")
    ss = ssgsea_scores(big_df, pathways)
    zs = zscore_scores(big_df, pathways)
    corr = correlate_methods(ss, zs)
    # IFN_RESPONSE should have positive correlation (both methods agree on up in A)
    r_ifn = corr.loc["IFN_RESPONSE", "spearman_r"]
    assert r_ifn > 0.5, f"Expected positive correlation for injected pathway, got {r_ifn}"
    assert r_ifn <= 1.0
    # Check n_samples correct
    assert corr.loc["IFN_RESPONSE", "n_samples"] == 12
    # Pearson also positive
    assert corr.loc["IFN_RESPONSE", "pearson_r"] > 0.5


def test_correlation_handles_constant_and_missing():
    ss = pd.DataFrame({"PW1": [1,1,1], "PW2": [1,2,3]}, index=["S1","S2","S3"])
    zs = pd.DataFrame({"PW1": [1,1,1], "PW2": [1,2,3]}, index=["S1","S2","S3"])
    corr = correlate_methods(ss, zs)
    # PW1 constant => nan
    assert np.isnan(corr.loc["PW1", "pearson_r"])
    # PW2 perfect correlation
    assert abs(corr.loc["PW2", "pearson_r"] - 1.0) < 1e-9
    assert abs(corr.loc["PW2", "spearman_r"] - 1.0) < 1e-9
