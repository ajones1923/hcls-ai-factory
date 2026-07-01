"""Cohort/population endpoint (api/routes/cohort) — aggregates over all 50 projections."""
from fastapi.testclient import TestClient

from api.main import app


def test_cohort_population_summary():
    with TestClient(app) as c:
        d = c.get("/cohort").json()
    assert d["n_patients"] == 50
    assert len(d["patients"]) == 50
    # mosaic recoveries are the headline scale story — the synthetic cohort has 7
    assert d["highlights"]["mosaic_recovered"] == 7
    # classification distribution covers the whole panel
    assert sum(d["distributions"]["classification"].values()) == 50
    # featured A/B/C sort to the top and carry their tags
    top = {r["featured"] for r in d["patients"][:3]}
    assert top == {"A", "B", "C"}
    # at least one patient shows a lesion flag (B crossing or C at threshold)
    assert any(r["lesion_flags"] for r in d["patients"])
