"""Fail-closed release gate for approved-artifact-dependent features.

Design choice documented:
- Core statistical scoring (ssGSEA, z-score) are deterministic, not trained ML models,
  so they can run WITHOUT the gate. They are exposed under /api/v1/score/* ungated.
- Features that depend on a specific curated pathway-database revision / validated
  analysis artifact (e.g., precomputed Hallmark Reactome reference, promoted dashboard
  bundle) are GATED. They require:
      MODEL_RELEASE_APPROVED=true
      APPROVED_ARTIFACT_REVISION=<non-empty revision string, e.g. git SHA or version>
  Without these, gated endpoints return 403 with honest error and health/readiness
  reflect the blocked state. Never fabricate readiness.

Environment variables:
  MODEL_RELEASE_APPROVED : "true"/"1"/"yes" (case-insensitive) means approved
  APPROVED_ARTIFACT_REVISION : opaque string, must be non-empty when approved
"""

from __future__ import annotations

import os

from fastapi import HTTPException


def _is_truthy(val: str | None) -> bool:
    if val is None:
        return False
    return val.strip().lower() in ("true", "1", "yes", "y", "on")


def is_release_approved() -> bool:
    """Check if release gate is satisfied."""
    approved = _is_truthy(os.getenv("MODEL_RELEASE_APPROVED"))
    revision = os.getenv("APPROVED_ARTIFACT_REVISION", "")
    # Revision must be non-empty, non-whitespace
    has_revision = bool(revision and revision.strip())
    return approved and has_revision


def get_approved_revision() -> str | None:
    rev = os.getenv("APPROVED_ARTIFACT_REVISION")
    if rev and rev.strip():
        return rev.strip()
    return None


def require_release_approved():
    """FastAPI dependency: raises 403 if gate not satisfied (fail-closed)."""
    if not is_release_approved():
        rev = get_approved_revision()
        approved_flag = os.getenv("MODEL_RELEASE_APPROVED")
        raise HTTPException(
            status_code=403,
            detail={
                "error": "Release gate closed: artifact not approved for serving",
                "hint": "Set MODEL_RELEASE_APPROVED=true and APPROVED_ARTIFACT_REVISION=<rev> to enable curated artifact endpoints. Scoring endpoints (/api/v1/score/*) remain available without gate.",
                "approved_flag": approved_flag,
                "revision_present": bool(rev),
            },
        )
    return True


def gate_status() -> dict:
    """Honest status dict for health/readiness endpoints."""
    approved_flag = os.getenv("MODEL_RELEASE_APPROVED")
    revision = get_approved_revision()
    approved = is_release_approved()
    return {
        "approved": approved,
        "approved_flag_raw": approved_flag,
        "revision": revision,
        "gate_satisfied": approved,
        "message": "Gate satisfied: approved artifact serving enabled"
        if approved
        else "Gate closed: curated artifact endpoints disabled (scoring endpoints still available)",
    }
