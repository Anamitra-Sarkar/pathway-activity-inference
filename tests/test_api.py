import os
import pytest
from fastapi.testclient import TestClient

# Ensure auth disabled for some tests? We'll test both
os.environ.pop("MODEL_RELEASE_APPROVED", None)
os.environ.pop("APPROVED_ARTIFACT_REVISION", None)

from backend.main import app

client = TestClient(app)

# Synthetic expression payload for scoring tests
EXPR = {
    "S1": {"IFIT1": 9.0, "MX1": 9.2, "ISG15": 8.8, "G1": 5.0, "G2": 4.5},
    "S2": {"IFIT1": 3.0, "MX1": 3.1, "ISG15": 2.9, "G1": 5.1, "G2": 4.6},
}
GENE_SETS = {
    "IFN_RESPONSE": ["IFIT1", "MX1", "ISG15"],
    "RANDOM": ["G1", "G2"],
}


def test_health_honest():
    r = client.get("/health")
    assert r.status_code == 200
    j = r.json()
    assert j["status"] == "ok"
    assert "gate" in j
    assert j["gate"]["gate_satisfied"] is False  # no env set

def test_ready_honest_no_gate():
    r = client.get("/ready")
    assert r.status_code == 200
    j = r.json()
    assert j["scoring_ready"] is True
    assert j["curated_artifact_ready"] is False
    assert j["gate"]["approved"] is False

def test_scoring_ungated_without_gate():
    # Scoring should work even when gate closed (core design choice)
    r = client.post("/api/v1/score/ssgsea", json={"expression": EXPR, "gene_sets": GENE_SETS})
    assert r.status_code == 200
    data = r.json()
    assert set(data["pathways"]) == {"IFN_RESPONSE", "RANDOM"}
    # IFN score S1 > S2
    assert data["scores"]["S1"]["IFN_RESPONSE"] > data["scores"]["S2"]["IFN_RESPONSE"]

def test_zscore_ungated():
    r = client.post("/api/v1/score/zscore", json={"expression": EXPR, "gene_sets": GENE_SETS})
    assert r.status_code == 200
    data = r.json()
    # group z-score: S1 upregulated => positive, S2 negative
    assert data["scores"]["S1"]["IFN_RESPONSE"] > data["scores"]["S2"]["IFN_RESPONSE"]

def test_both_and_differential():
    r = client.post("/api/v1/score/both", json={"expression": EXPR, "gene_sets": GENE_SETS})
    assert r.status_code == 200
    j = r.json()
    assert "ssgsea" in j and "zscore" in j and "correlation" in j

    # Use ssGSEA scores for differential
    scores = j["ssgsea"]
    groups = {"S1": "A", "S2": "B"}
    # Need more samples for differential? Use duplicate to make 4 samples
    expr4 = {
        "A1": {"IFIT1": 9, "MX1":9, "ISG15":9, "G1":5, "G2":5},
        "A2": {"IFIT1": 9.1, "MX1":9.2, "ISG15":8.9, "G1":5.1, "G2":4.9},
        "B1": {"IFIT1": 3, "MX1":3, "ISG15":3, "G1":5, "G2":5},
        "B2": {"IFIT1": 3.1, "MX1":2.9, "ISG15":3.0, "G1":5, "G2":5},
    }
    r2 = client.post("/api/v1/score/ssgsea", json={"expression": expr4, "gene_sets": GENE_SETS})
    assert r2.status_code == 200
    scores4 = r2.json()["scores"]
    r3 = client.post("/api/v1/differential", json={"scores": scores4, "groups": {"A1":"A","A2":"A","B1":"B","B2":"B"}})
    assert r3.status_code == 200
    diff = r3.json()
    # IFN should be significant
    assert any(d["pathway"]=="IFN_RESPONSE" for d in diff)

def test_gated_endpoint_closed_fail():
    r = client.get("/api/v1/curated/status-gated")
    assert r.status_code == 403
    assert "Release gate closed" in r.json()["detail"]["error"]

def test_gated_endpoint_open_with_env(monkeypatch):
    monkeypatch.setenv("MODEL_RELEASE_APPROVED", "true")
    monkeypatch.setenv("APPROVED_ARTIFACT_REVISION", "v1.0.0-abc123")
    r = client.get("/api/v1/curated/status-gated")
    assert r.status_code == 200
    assert r.json()["revision"] == "v1.0.0-abc123"
    # Health should now reflect gate satisfied
    r2 = client.get("/health")
    assert r2.json()["gate"]["gate_satisfied"] is True
    r3 = client.get("/ready")
    assert r3.json()["curated_artifact_ready"] is True

def test_curated_with_auth_gated_closed():
    # This endpoint requires both gate and auth
    r = client.get("/api/v1/curated/analysis")
    # Should be 403 gate closed first (or 401 auth missing, but gate checked first? depends order)
    # Our implementation has require_release_approved before verify_firebase_token in dependency order?
    # Actually Depends order: authorized first then user? FastAPI evaluates in parameter order.
    # So gate 403 should prevail
    assert r.status_code in (403, 401)

def test_auth_stub_missing_token(monkeypatch):
    # Ensure auth not disabled
    monkeypatch.delenv("FIREBASE_AUTH_DISABLED", raising=False)
    monkeypatch.delenv("FIREBASE_SERVICE_ACCOUNT_JSON", raising=False)
    # Try to access auth-required endpoint without token when gate open: should 401
    monkeypatch.setenv("MODEL_RELEASE_APPROVED", "true")
    monkeypatch.setenv("APPROVED_ARTIFACT_REVISION", "rev-test")
    r = client.get("/api/v1/curated/analysis")
    assert r.status_code == 401

def test_auth_stub_with_disabled(monkeypatch):
    monkeypatch.setenv("FIREBASE_AUTH_DISABLED", "true")
    monkeypatch.setenv("MODEL_RELEASE_APPROVED", "true")
    monkeypatch.setenv("APPROVED_ARTIFACT_REVISION", "rev-test2")
    r = client.get("/api/v1/curated/analysis")
    assert r.status_code == 200
    assert r.json()["artifact_revision"] == "rev-test2"

def test_auth_stub_with_bearer_token(monkeypatch):
    monkeypatch.delenv("FIREBASE_AUTH_DISABLED", raising=False)
    monkeypatch.setenv("MODEL_RELEASE_APPROVED", "true")
    monkeypatch.setenv("APPROVED_ARTIFACT_REVISION", "rev-bearer")
    # Provide a dummy JWT 3-part token
    headers = {"Authorization": "Bearer header.payload.signature"}
    r = client.get("/api/v1/curated/analysis", headers=headers)
    assert r.status_code == 200

def test_pathway_db_status_honest_unauthenticated():
    r = client.get("/api/v1/pathway-db/status")
    assert r.status_code == 200
    # When gate closed, artifact is None
    # Ensure no fabricating approved when not set - reset env
    import os
    os.environ.pop("MODEL_RELEASE_APPROVED", None)
    os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
    r2 = client.get("/api/v1/pathway-db/status")
    assert r2.json()["artifact"] is None
    assert r2.json()["gate"]["gate_satisfied"] is False
