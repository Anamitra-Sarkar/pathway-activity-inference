"""Shared synthetic fixtures for pathway activity tests.

Synthetic expression matrix with INJECTED signal (honest: not a real biological finding):
- 4 samples: A1, A2 (group A), B1, B2 (group B)
- Genes: 10 background + 3 IFN pathway genes (IFIT1, MX1, ISG15)
- Group A: IFN genes upregulated (mean 9) vs Group B (mean 3)
- Background genes: ~5 in both groups (+ noise)

This verifies pipeline recovers known injected signal.
"""

import pandas as pd
import numpy as np


def make_injected_fixture():
    genes = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "IFIT1", "MX1", "ISG15"]
    pathways = {
        "IFN_RESPONSE": {"description": "Injected interferon response", "genes": ["IFIT1", "MX1", "ISG15"]},
        "RANDOM_SET": {"description": "Random background", "genes": ["G1", "G2"]},
        "P53_PATHWAY": {"description": "Unperturbed", "genes": ["G3", "G4", "G5"]},
    }
    # Fixed seed for determinism
    rng = np.random.default_rng(42)
    # Build expression matrix samples x genes
    samples = ["A1", "A2", "B1", "B2"]
    data = {}
    for s in samples:
        row = {}
        is_A = s.startswith("A")
        for g in genes:
            if g in ("IFIT1", "MX1", "ISG15"):
                base = 9.0 if is_A else 3.0
            else:
                base = 5.0
            # add small noise
            row[g] = float(base + rng.normal(0, 0.3))
        data[s] = row
    expr_df = pd.DataFrame.from_dict(data, orient="index")
    labels = pd.Series({"A1": "A", "A2": "A", "B1": "B", "B2": "B"})
    return expr_df, pathways, labels
