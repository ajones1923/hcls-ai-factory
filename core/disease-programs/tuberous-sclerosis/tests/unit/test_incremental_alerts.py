"""Incremental note dispatch (FR-OR-3) and the per-clinician alert sliding window (FR-SF-3)."""
from src.orchestrator.events import EventType
from src.orchestrator.graph import Orchestrator
from src.orchestrator.state import SqliteEventStore


def _enrolled():
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    cohort = {"phenotypes": [{"hpo_id": "HP:0009717", "label": "Cortical tubers", "onset": "infancy"}],
              "genomics": {"gene": "TSC2", "variant": "c.3037C>T", "vaf": 0.5},
              "notes": []}
    orch.enroll("P", cohort)
    return orch, cohort


def test_ingest_note_reruns_only_phenome_and_tand():
    orch, cohort = _enrolled()
    before = orch.store.events_for("P")
    n_variant_before = sum(e.event_type == EventType.variant_curated for e in before)
    n_traj_before = sum(e.event_type == EventType.trajectory_forecast for e in before)

    text = "[SYNTHETIC] new shagreen patch noted."
    s = text.index("shagreen patch")
    note = {"specialty": "neurology", "date": "2025-12-01", "text": text,
            "spans": [{"start": s, "end": s + len("shagreen patch"), "quote": "shagreen patch",
                       "kind": "hpo", "hpo_id": "HP:0009721", "label": "Shagreen patch",
                       "polarity": "present", "temporality": "current"}]}
    ran = orch.ingest_note("P", note, cohort)
    assert ran == ["phenome_mapper", "tand_surveillance"]

    after = orch.store.events_for("P")
    # exactly one new phenome + one new tand event; variant/trajectory NOT re-run
    assert sum(e.event_type == EventType.phenome_mapped for e in after) == \
        sum(e.event_type == EventType.phenome_mapped for e in before) + 1
    assert sum(e.event_type == EventType.variant_curated for e in after) == n_variant_before
    assert sum(e.event_type == EventType.trajectory_forecast for e in after) == n_traj_before
    # the new phenotype from the note is now in the projection
    ids = {t["hpo_id"] for t in orch.store.projection("P")["hpo_profile"]["hpo_terms"]}
    assert "HP:0009721" in ids


def test_alert_sliding_window_triggers_recalibration():
    orch, _ = _enrolled()
    budget = None
    for _ in range(3):
        r = orch.deliver_alert("dr-a")
        budget = r["weekly_budget"]
    assert orch.deliver_alert("dr-a")["over_budget"] is False or budget == 3
    # one past budget -> recalibrate
    over = None
    for _ in range(5):
        over = orch.deliver_alert("dr-a")
    assert over["over_budget"] is True and over["recalibrate"] is True
    # a different clinician is independent
    assert orch.deliver_alert("dr-b")["week_count"] == 1


def test_alert_surface_reports_window_readonly():
    orch, _ = _enrolled()
    surf = orch.assemble_surface("P", "alerts")
    assert "clinician_week_count" in surf and "weekly_budget" in surf
    # reading the surface does not itself increment the window
    c1 = orch.assemble_surface("P", "alerts")["clinician_week_count"]
    c2 = orch.assemble_surface("P", "alerts")["clinician_week_count"]
    assert c1 == c2
