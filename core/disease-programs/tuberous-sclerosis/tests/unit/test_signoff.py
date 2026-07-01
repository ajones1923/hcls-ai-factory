"""Variant-curator AI-labeled draft report + molecular-geneticist sign-off gate (FR-VC-5)."""
from fastapi.testclient import TestClient

from api.main import app
from app._engine import featured


def test_draft_is_ai_labeled_and_written(tmp_path):
    pid = featured()["B"]
    with TestClient(app) as c:
        d = c.get(f"/variant-curator/{pid}/draft").json()
    assert "DRAFT" in d["markdown"] and "molecular-geneticist" in d["markdown"].lower()
    assert isinstance(d["review_status"], str)      # pending or, if already signed, signed_off:*
    assert d["path"].endswith("variant_draft.md")
    # the ACMG classification appears in the rendered report
    assert "ACMG-AMP" in d["markdown"]


def test_sign_off_records_human_gate():
    pid = featured()["B"]
    with TestClient(app) as c:
        r = c.post(f"/variant-curator/{pid}/sign-off",
                   json={"signer": "Dr. Geneticist", "decision": "approve"}).json()
        assert r["review_status"] == "signed_off:approve" and r["signer"] == "Dr. Geneticist"
        status = c.get(f"/variant-curator/{pid}/sign-off").json()
    assert status["review_status"] == "signed_off:approve"


def test_sign_off_rejects_bad_decision():
    pid = featured()["C"]
    with TestClient(app) as c:
        resp = c.post(f"/variant-curator/{pid}/sign-off", json={"signer": "x", "decision": "bogus"})
    assert resp.status_code == 400
