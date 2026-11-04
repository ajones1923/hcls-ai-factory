"""A differential diagnosis must say how much evidence stood behind it.

Found by the claim-class sweep, 2026-11-04 (class B), after the same defect was found by
hand in precision-oncology. Seven subjects share it; this is one of three byte-identical
copies.

`RareDiseaseAgent.evaluate_evidence` grades retrieval `sufficient` / `partial` /
`insufficient`. The verdict decided whether to widen the search with sub-questions and was
then dropped on the floor -- `run()` returned the `DiagnosticResult` with no record of how
grounded it was, and unlike precision-oncology there was not even a log line. A differential
built on almost nothing looked exactly like one built on good evidence.

Two things are fixed here, and the second is not cosmetic:

  * the verdict is carried on `DiagnosticResult.evidence_sufficiency` and printed in the
    report header, with an explicit warning block when it is `insufficient`;
  * the grade is recomputed AFTER the sub-question expansion. The original call graded the
    evidence that *triggered* widening, so even had it been reported it would have described
    a state the answer was no longer based on.

`evidence_sufficiency` is deliberately separate from `DiagnosticResult.confidence`:
confidence scores the candidate ranking, this scores the retrieval the ranking came from. A
high confidence over thin evidence is exactly the combination worth being able to see.
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models import DiagnosticResult  # noqa: E402


class TestTheFieldExists:
    def test_diagnostic_result_carries_the_verdict(self):
        assert DiagnosticResult().evidence_sufficiency == "unknown"

    def test_it_is_not_the_same_thing_as_confidence(self):
        """A confident ranking over thin evidence must be visible as exactly that."""
        r = DiagnosticResult(confidence=0.91)
        r.evidence_sufficiency = "insufficient"
        assert r.confidence == 0.91
        assert r.evidence_sufficiency == "insufficient"


class TestItReachesTheReport:
    @pytest.fixture
    def agent(self):
        from src.agent import RareDiseaseAgent
        with patch.object(RareDiseaseAgent, "__init__", lambda self: None):
            a = RareDiseaseAgent()
        a.rag = MagicMock()
        return a

    def test_the_header_states_the_sufficiency(self, agent):
        r = DiagnosticResult()
        r.evidence_sufficiency = "partial"
        report = agent.generate_report(r)
        assert "Evidence sufficiency" in report
        assert "partial" in report

    def test_insufficient_gets_an_explicit_warning_block(self, agent):
        r = DiagnosticResult()
        r.evidence_sufficiency = "insufficient"
        report = agent.generate_report(r)
        assert "Evidence was insufficient" in report
        assert "provisional" in report

    def test_a_well_grounded_result_gets_no_warning(self, agent):
        """The warning has to mean something, so it must not appear on every report."""
        r = DiagnosticResult()
        r.evidence_sufficiency = "sufficient"
        report = agent.generate_report(r)
        assert "Evidence sufficiency" in report
        assert "Evidence was insufficient" not in report


class TestRunWiresItUp:
    def test_run_attaches_the_verdict(self):
        src = (ROOT / "src" / "agent.py").read_text()
        assert "response.evidence_sufficiency = quality" in src, (
            "run() computes `quality` and must attach it; without this line the grade is "
            "dropped again")

    def test_the_grade_is_recomputed_after_widening(self):
        """Grading only before expansion would report a state the answer is not based on."""
        src = (ROOT / "src" / "agent.py").read_text()
        head, _, tail = src.partition("response.results.extend(sub_response)")
        assert tail, "the sub-question expansion block moved; re-read this"
        before_attach = tail.split("response.evidence_sufficiency")[0]
        assert "self.evaluate_evidence(response.results)" in before_attach, (
            "the verdict must be recomputed after the expansion, not carried from before it")
