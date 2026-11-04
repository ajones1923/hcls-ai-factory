"""Tests that machine-readable output does not claim work or finality it lacks.

Two surfaces:

* the SSE stream, which narrated "Checking ClinVar records", "Evaluating
  population frequency", "Running in-silico predictors", "Applying ACMG
  criteria" over ``asyncio.sleep(0.8)`` and ended ``"status": "completed"``,
  ``"<id> workflow completed successfully"`` (audit §14);
* the report routes, which emitted a FHIR R4 DiagnosticReport with
  ``status: "final"`` and an empty ``conclusion``, a Phenopacket asserting an
  HPO ontology version for an agent that loads no ontology, unvalidated HPO
  IDs, a "Classification: VUS" default, and no research-use statement in any
  format (audit §23).

Primary source for the FHIR status change: HL7 FHIR R4
``DiagnosticReportStatus`` value set, retrieved 2026-11-03 from
https://hl7.org/fhir/R4/valueset-diagnostic-report-status.html --
``final`` = "The report is complete and verified by an authorized person";
``partial`` = "This is a partial (e.g. initial, interim or preliminary)
report: data in the report may be incomplete or unverified";
``registered`` = "The existence of the report is registered, but there is
nothing yet available."
"""

import inspect

import pytest

from api.routes import events, reports


# ═══════════════════════════════════════════════════════════════════════════
# §14 — the SSE stream must not narrate work it did not do
# ═══════════════════════════════════════════════════════════════════════════


def _collect(gen):
    import asyncio

    async def run():
        return [frame async for frame in gen]

    return asyncio.run(run())


@pytest.fixture
def demo_frames():
    return _collect(
        events._demo_progress_generator("variant_interpretation", {}, None)
    )


def test_every_demo_frame_is_marked_simulated(demo_frames):
    """No frame on the wire can be read as evidence of real work.

    The word "Simulate" appeared only in the internal docstring and an inline
    comment; nothing in the emitted payload said so, and the public endpoint
    docstring called it "real-time workflow progress".
    """
    assert demo_frames
    for frame in demo_frames:
        assert '"simulated": true' in frame.lower(), frame


def test_no_frame_claims_the_workflow_completed(demo_frames):
    """The terminal frame must not say "completed" or "successfully".

    Clinical consequence of the old payload: a clinician opens the variant
    review page for a BRCA2 variant, connects
    ``EventSource('/v1/events/stream?workflow_id=variant_interpretation')``
    and watches six steps tick to 100%, ending "variant_interpretation
    workflow completed successfully". Nothing was computed about the variant.
    They now have positive, specific, step-by-step assurance that ClinVar was
    checked, population frequency evaluated, in-silico predictors run and ACMG
    criteria applied -- the four things they would otherwise verify by hand.
    """
    joined = " ".join(demo_frames)
    assert '"status": "completed"' not in joined
    assert "completed successfully" not in joined
    assert "NOTHING WAS COMPUTED" in joined


def test_step_labels_are_marked_simulated(demo_frames):
    """A step that names a real evidence source says it is simulated."""
    acmg = [f for f in demo_frames if "ACMG" in f]
    assert acmg, "the ACMG step should still be present, but labelled"
    for frame in acmg:
        assert "[SIMULATED]" in frame


def test_the_generator_calls_no_engine():
    """The generator body must not reach into any engine or collection.

    It is 44 lines of sleep and yield. This pins that: if someone later wires
    real work into it, the honesty markers above have to be revisited at the
    same time.
    """
    import ast
    import textwrap

    tree = ast.parse(textwrap.dedent(inspect.getsource(events._demo_progress_generator)))
    awaited = [
        ast.unparse(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Await)
    ]
    # The ONLY thing this generator awaits is a sleep. If real work is ever
    # wired in here, this fails and the honesty markers above must be
    # revisited in the same change.
    assert awaited, "expected at least one await"
    for expr in awaited:
        assert expr.startswith("asyncio.sleep("), expr


def test_stream_endpoint_requires_an_explicit_demo_flag():
    """A workflow_id alone no longer yields a success stream.

    ``workflow_id`` was unvalidated, so ANY string produced a step-by-step
    success stream for work that never ran.
    """
    src = inspect.getsource(events.event_stream)
    assert "demo" in src
    assert "HTTPException" in src
    assert "status_code=400" in src


def test_public_endpoint_docstring_does_not_claim_real_time_progress():
    """The docstring said "real-time workflow progress"."""
    doc = events.event_stream.__doc__ or ""
    assert "real-time workflow progress" not in doc
    assert "DEMO" in doc


# ═══════════════════════════════════════════════════════════════════════════
# §23 — FHIR, Phenopacket and markdown honesty
# ═══════════════════════════════════════════════════════════════════════════


def test_fhir_report_is_never_status_final():
    """``final`` asserts verification by an authorized person. None happens."""
    r = reports._generate_fhir_diagnostic_report(
        {"classification": "Likely pathogenic"}, "Variant report", "p1",
    )
    assert r["status"] != "final"
    assert r["status"] == "partial"


def test_fhir_report_with_no_conclusion_is_registered_and_omits_it():
    """An empty conclusion is not a conclusion of "nothing found".

    Clinical consequence of the old output: ``POST /v1/reports/generate
    {"format":"fhir"}`` after a ``/diagnose`` call that returned an empty
    differential emitted a FHIR R4 DiagnosticReport with ``status: "final"``,
    LOINC 81247-9 ("Master HL7 genetic variant reporting panel") and
    ``conclusion: ""``. An EHR ingesting it files a FINAL genetic diagnostic
    report with an empty conclusion into the patient's chart, where it is read
    as "genetic testing performed, nothing found" and suppresses re-testing.
    """
    r = reports._generate_fhir_diagnostic_report({}, "Variant report", "p1")
    assert r["status"] == "registered"
    assert "conclusion" not in r


def test_fhir_report_carries_the_research_use_statement():
    """The narrative is where a conforming renderer shows it."""
    r = reports._generate_fhir_diagnostic_report(
        {"classification": "VUS"}, "Variant report", "p1",
    )
    assert "RESEARCH USE ONLY" in r["text"]["div"]
    assert "not verified by a qualified person" in r["text"]["div"]


def test_phenopacket_does_not_assert_an_hpo_version():
    """The agent loads no ontology, so it declares no ontology release.

    ``"version": "2024-04-26"`` was asserted for an agent whose entire HPO
    frequency knowledge is a 43-entry hand-written dict. Declaring a resource
    version declares that the terms in the packet were resolved against that
    release.
    """
    pp = reports._generate_phenopacket({}, "Phenotype report", "p1")
    assert pp["metaData"]["resources"] == []


def test_phenopacket_carries_the_research_use_statement():
    pp = reports._generate_phenopacket({}, "Phenotype report", "p1")
    refs = " ".join(
        r.get("description", "") for r in pp["metaData"]["externalReferences"]
    )
    assert "RESEARCH USE ONLY" in refs


INVALID_HPO_IDS = ["HP:9999", "banana", "", "HP:00012500", "hp:0001250", "0001250"]


@pytest.mark.parametrize("term_id", INVALID_HPO_IDS)
def test_phenopacket_rejects_malformed_hpo_ids(term_id):
    """A malformed ID must not reach phenotypicFeatures[].type.id.

    ``p.get("id", "")`` was written through with no validation, so
    ``HP:9999999`` and the literal string ``banana`` -- both accepted
    unvalidated by ``/phenotype/match`` -- landed in an interchange format. A
    mis-coded phenotype then propagates and cannot be caught by reading the
    label, because the label is right and the ID is wrong.
    """
    pp = reports._generate_phenopacket(
        {"phenotypes": [{"id": term_id, "label": "Seizure"}]},
        "Phenotype report", "p1",
    )
    ids = [f["type"]["id"] for f in pp["phenotypicFeatures"]]
    assert term_id not in ids


def test_phenopacket_reports_what_it_rejected():
    """Dropped terms are named, not silently discarded."""
    pp = reports._generate_phenopacket(
        {"phenotypes": [{"id": "banana", "label": "Seizure"}]},
        "Phenotype report", "p1",
    )
    refs = " ".join(
        r.get("description", "") for r in pp["metaData"]["externalReferences"]
    )
    assert "banana" in refs


def test_phenopacket_accepts_a_well_formed_hpo_id():
    """HP:0001250 (Seizure) still goes through."""
    pp = reports._generate_phenopacket(
        {"phenotypes": [{"id": "HP:0001250", "label": "Seizure"}]},
        "Phenotype report", "p1",
    )
    assert pp["phenotypicFeatures"][0]["type"]["id"] == "HP:0001250"


def test_markdown_header_carries_the_research_use_statement():
    """Every markdown report leads with the disclaimer."""
    header = reports._generate_markdown_header("Variant Interpretation Report")
    assert "RESEARCH USE ONLY" in header


def test_variant_markdown_does_not_default_to_vus():
    """A report with no classification must not print "VUS".

    ``data.get('classification', 'VUS')`` asserted a call nobody made.
    """
    md = reports._variant_interpretation_markdown({})
    assert "**VUS**" not in md
    assert "not classified" in md
