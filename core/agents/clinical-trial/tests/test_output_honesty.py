"""Tests that machine-readable output does not claim work it did not do.

Two surfaces:

* the SSE stream (audit F15), which narrated "Screening eligibility criteria",
  "Scoring inclusion criteria", "Checking exclusion criteria", "Ranking
  matches" -- and for competitive_intel, "Scanning trial registries" -- over
  ``asyncio.sleep(0.8)``, ending ``"status": "completed"``, ``"<id> workflow
  completed successfully"``. This stream was the only place in the product that
  CLAIMED criteria were checked, and the one place guaranteed not to have
  checked them;
* the report routes (audit F23), which printed a bold Score column with no
  criteria columns and no research-use statement, while the good ``_DISCLAIMER``
  sat in a ``src/export.py`` the HTTP routes do not use.
"""

import ast
import asyncio
import inspect
import textwrap

import pytest

from api.routes import events, reports


def _collect(gen):
    async def run():
        return [frame async for frame in gen]

    return asyncio.run(run())


@pytest.fixture
def demo_frames():
    return _collect(events._demo_progress_generator("patient_matching", {}, None))


# ═══════════════════════════════════════════════════════════════════════════
# F15 — the SSE stream
# ═══════════════════════════════════════════════════════════════════════════


def test_every_demo_frame_is_marked_simulated(demo_frames):
    """No frame on the wire can be read as evidence of real work."""
    assert demo_frames
    for frame in demo_frames:
        assert '"simulated": true' in frame.lower(), frame


def test_no_frame_claims_the_workflow_completed(demo_frames):
    """The terminal frame must not say "completed" or "successfully".

    Clinical consequence of the old payload: a coordinator opens the matcher,
    watches "Screening eligibility criteria -> Scoring inclusion criteria ->
    Checking exclusion criteria -> Ranking matches -> 100% -- patient_matching
    workflow completed successfully", and reasonably concludes that this
    patient's inclusion and exclusion criteria were checked against the trials
    shown. Nothing was.
    """
    joined = " ".join(demo_frames)
    assert '"status": "completed"' not in joined
    assert "completed successfully" not in joined
    assert "NOTHING WAS COMPUTED" in joined


def test_eligibility_step_labels_are_marked_simulated(demo_frames):
    """The steps that assert criterion checking say they are simulated."""
    for needle in ("eligibility criteria", "exclusion criteria"):
        frames = [f for f in demo_frames if needle in f.lower()]
        assert frames, needle
        for frame in frames:
            assert "[SIMULATED]" in frame


def test_registry_scan_step_is_marked_simulated():
    """"Scanning trial registries" asserted a CTG query that never happened."""
    frames = _collect(
        events._demo_progress_generator("competitive_intel", {}, None)
    )
    scan = [f for f in frames if "trial registries" in f.lower()]
    assert scan
    for frame in scan:
        assert "[SIMULATED]" in frame


def test_the_generator_only_ever_awaits_a_sleep():
    """The generator does no work.

    If real work is ever wired in here, this fails and the honesty markers
    above must be revisited in the same change.
    """
    tree = ast.parse(
        textwrap.dedent(inspect.getsource(events._demo_progress_generator))
    )
    awaited = [
        ast.unparse(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Await)
    ]
    assert awaited
    for expr in awaited:
        assert expr.startswith("asyncio.sleep("), expr


def test_stream_endpoint_requires_an_explicit_demo_flag():
    """``workflow_id`` was unvalidated: any string produced a success stream."""
    src = inspect.getsource(events.event_stream)
    assert "demo" in src
    assert "status_code=400" in src


def test_public_endpoint_docstring_does_not_claim_real_time_progress():
    """The docstring said "real-time workflow progress"."""
    doc = events.event_stream.__doc__ or ""
    assert "real-time workflow progress" not in doc
    assert "DEMO" in doc


# ═══════════════════════════════════════════════════════════════════════════
# F23 — report routes
# ═══════════════════════════════════════════════════════════════════════════

_MATCH_DATA = {
    "total_screened": 20,
    "criteria_evaluated": False,
    "matches": [{
        "trial_id": "NCT02511106",
        "trial_title": "ADAURA: Osimertinib as Adjuvant Therapy in EGFR-Mutant NSCLC",
        "phase": None,
        "status": None,
        "retrieval_similarity": 0.91,
        "criteria_evaluated": False,
        "inclusion_met": None,
        "inclusion_total": None,
        "exclusion_clear": None,
        "exclusion_total": None,
    }],
}


def test_markdown_header_carries_the_research_use_statement():
    """A grep for "disclaim|research use|not for clinical" across the trial
    agent's UI, api/main.py, trial_clinical.py and reports.py returned NO hits.
    """
    header = reports._generate_markdown_header("Patient Matching Report")
    assert "RESEARCH USE ONLY" in header
    assert "NOT an eligibility determination" in header


def test_matching_markdown_has_criteria_columns():
    """A Score column with no criteria columns is the whole defect.

    The table printed Trial ID, Title, Phase, **Score** and Confidence. The one
    number a reader's eye lands on was the one with no eligibility content
    behind it.
    """
    md = reports._matching_report_markdown(_MATCH_DATA)
    assert "Inclusion met" in md
    assert "Exclusion clear" in md
    assert "Criteria evaluated" in md


def test_matching_markdown_does_not_bold_the_score():
    """The score is not the emphasised cell any more."""
    md = reports._matching_report_markdown(_MATCH_DATA)
    assert "**0.91**" not in md


def test_matching_markdown_states_that_no_criteria_were_evaluated():
    md = reports._matching_report_markdown(_MATCH_DATA)
    assert "No eligibility criteria were evaluated" in md
    assert "not evaluated" in md


def test_matching_markdown_does_not_invent_phase_or_status():
    """Absent phase and status print as unknown, not as real values."""
    md = reports._matching_report_markdown(_MATCH_DATA)
    assert "unknown" in md
    assert "phase_ii" not in md
    assert "recruiting" not in md


def test_matching_markdown_names_a_missing_registry_id():
    """A row with no NCT ID says so rather than carrying a minted one."""
    data = dict(_MATCH_DATA)
    data["matches"] = [dict(_MATCH_DATA["matches"][0], trial_id=None)]
    md = reports._matching_report_markdown(data)
    assert "no registry ID" in md
    assert "NCT-RESULT" not in md


def test_matching_markdown_does_not_label_retrieval_as_total_screened():
    """"Total Screened: N" was the retrieval hit count, not a screening count."""
    md = reports._matching_report_markdown(_MATCH_DATA)
    assert "Total Screened" not in md
    assert "Candidate records retrieved" in md
