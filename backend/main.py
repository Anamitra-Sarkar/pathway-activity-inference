"""FastAPI backend for pathway activity inference.

- Health/readiness honestly reflect gate state, never fabricate.
- Scoring endpoints (ssgsea, zscore) are ungated (deterministic stats).
- Curated artifact endpoint is gated fail-closed.
"""

from __future__ import annotations

import os
from typing import Dict, List

import pandas as pd
from fastapi import Depends, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from .auth import require_auth, verify_firebase_token, optional_auth
from .release_gate import gate_status, is_release_approved, require_release_approved

# Reuse pipeline implementations
import sys
from pathlib import Path

# Ensure data_pipeline is importable when backend is run as module
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from data_pipeline.ssgsea import ssgsea_scores  # noqa: E402
from data_pipeline.zscore import zscore_scores  # noqa: E402
from data_pipeline.differential import differential_analysis  # noqa: E402
from data_pipeline.correlation import correlate_methods  # noqa: E402


app = FastAPI(
    title="Pathway Activity Inference API",
    version="0.1.0",
    description="Context-specific biological pathway activity inference (ssGSEA Barbie 2009, combined z-score Lee 2008).",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    """Liveness: always returns ok if process alive, plus honest gate status."""
    return {"status": "ok", "gate": gate_status(), "version": app.version}


@app.get("/ready")
def ready():
    """Readiness: honest. Scoring ready always; curated artifact ready only if gate satisfied."""
    gate = gate_status()
    # Scoring subsystem is stateless, always ready if process up
    scoring_ready = True
    curated_ready = gate["gate_satisfied"]
    # Overall ready if scoring ready; we report curated separately (never fabricate)
    return {
        "status": "ready" if scoring_ready else "not_ready",
        "scoring_ready": scoring_ready,
        "curated_artifact_ready": curated_ready,
        "gate": gate,
    }


@app.get("/api/v1/pathway-db/status")
def pathway_db_status():
    """Curated pathway DB status: gated info still honest even when closed, but full bundle gated."""
    gate = gate_status()
    if not is_release_approved():
        # Return honest status without leaking artifact content
        return {"gate": gate, "artifact": None, "detail": "Curated artifact not approved for serving"}
    # Gate satisfied: return artifact info
    rev = gate["revision"]
    return {
        "gate": gate,
        "artifact": {
            "revision": rev,
            "description": "Approved curated pathway database (MSigDB Hallmark / Reactome GMT) validated for serving",
            "pathways": ["HALLMARK_* (50)", "Reactome pathways"],
            "approved": True,
        },
    }


# --- Scoring endpoints (UNgated, deterministic) ---

class ScoreRequest(BaseModel):
    expression: Dict[str, Dict[str, float]]  # sample -> gene -> value
    gene_sets: Dict[str, List[str]]  # pathway -> genes
    alpha: float = 0.25


class ScoreResponse(BaseModel):
    samples: List[str]
    pathways: List[str]
    scores: Dict[str, Dict[str, float]]  # sample -> pathway -> score


def _to_df(expr: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    # expr: sample -> gene -> value
    df = pd.DataFrame.from_dict(expr, orient="index")
    return df


@app.post("/api/v1/score/ssgsea", response_model=ScoreResponse)
def score_ssgsea(req: ScoreRequest, user=Depends(optional_auth)):
    try:
        df = _to_df(req.expression)
        # Normalize gene_sets to pipeline format
        gs = {k: {"genes": v} for k, v in req.gene_sets.items()}
        scores = ssgsea_scores(df, gs, alpha=req.alpha)
        # scores: samples x pathways
        # Convert to sample->pathway dict
        out = {sample: scores.loc[sample].to_dict() for sample in scores.index}
        return ScoreResponse(samples=list(scores.index), pathways=list(scores.columns), scores=out)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/api/v1/score/zscore", response_model=ScoreResponse)
def score_zscore(req: ScoreRequest, user=Depends(optional_auth)):
    try:
        df = _to_df(req.expression)
        gs = {k: {"genes": v} for k, v in req.gene_sets.items()}
        scores = zscore_scores(df, gs)
        out = {sample: scores.loc[sample].to_dict() for sample in scores.index}
        return ScoreResponse(samples=list(scores.index), pathways=list(scores.columns), scores=out)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/api/v1/score/both")
def score_both(req: ScoreRequest, user=Depends(optional_auth)):
    try:
        df = _to_df(req.expression)
        gs = {k: {"genes": v} for k, v in req.gene_sets.items()}
        ss = ssgsea_scores(df, gs, alpha=req.alpha)
        zs = zscore_scores(df, gs)
        corr = correlate_methods(ss, zs)
        ss_out = {sample: ss.loc[sample].to_dict() for sample in ss.index}
        zs_out = {sample: zs.loc[sample].to_dict() for sample in zs.index}
        # Sanitise NaN to None for JSON compliance
        import numpy as np
        corr_out = corr.reset_index().to_dict(orient="records")
        for rec in corr_out:
            for k, v in list(rec.items()):
                if isinstance(v, float) and np.isnan(v):
                    rec[k] = None
        return {
            "samples": list(ss.index),
            "pathways": list(ss.columns),
            "ssgsea": ss_out,
            "zscore": zs_out,
            "correlation": corr_out,
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


class DifferentialRequest(BaseModel):
    scores: Dict[str, Dict[str, float]]  # sample -> pathway -> score
    groups: Dict[str, str]  # sample -> group
    group_a: str | None = None
    group_b: str | None = None


@app.post("/api/v1/differential")
def differential(req: DifferentialRequest, user=Depends(optional_auth)):
    try:
        scores_df = pd.DataFrame.from_dict(req.scores, orient="index")
        labels = pd.Series(req.groups)
        result = differential_analysis(scores_df, labels, group_a=req.group_a, group_b=req.group_b)
        return result.reset_index().to_dict(orient="records")
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


# --- Gated curated artifact endpoint (fail-closed) ---

@app.get("/api/v1/curated/analysis")
def curated_analysis(authorized: bool = Depends(require_release_approved), user: dict = Depends(verify_firebase_token)):
    """Gated endpoint: only serves when release approved AND authenticated.

    Depends on both release gate and Firebase auth (fail-closed on both).
    """
    rev = gate_status()["revision"]
    return {
        "artifact_revision": rev,
        "data": {
            "note": "This would serve precomputed validated pathway activity for approved revision",
            "pathways": ["HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_RESPONSE"],
            "sample_scores_preview": {},
        },
        "user": {"uid": user.get("uid")},
    }


# Also gated without auth but with gate only (alternative)
@app.get("/api/v1/curated/status-gated")
def curated_status_gated(authorized: bool = Depends(require_release_approved)):
    """Gated but not auth-protected: demonstrates gate pattern for artifact-dependent feature."""
    return {"status": "approved", "revision": gate_status()["revision"]}
