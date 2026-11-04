"""Clinical trial intelligence: an answer must say how much evidence stood behind it.

Found by the claim-class sweep, 2026-11-04 (class B), after the same defect was found by
hand in precision-oncology. Seven subjects share it; this is one of three byte-identical
copies of the block, alongside rare-disease-diagnostic and neurology.

`evaluate_evidence` grades retrieval `sufficient` / `partial` / `insufficient`. The verdict
decided whether to widen the search with sub-questions and was then dropped -- `run()`
returned `TrialResponse` with no record of how grounded it was, so an answer built on almost
nothing was indistinguishable from one built on good evidence.

Two fixes, and the second is not cosmetic:

  * the verdict is carried on `TrialResponse.evidence_sufficiency` and printed in the report
    header, with an explicit warning block when it is `insufficient`;
  * the grade is recomputed AFTER the sub-question expansion. The original call graded the
    evidence that *triggered* widening, so even had it been reported it would have
    described a state the answer was no longer based on.
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models import TrialResponse  # noqa: E402


class TestTheFieldExists:
    def test_the_response_model_carries_the_verdict(self):
        assert "evidence_sufficiency" in TrialResponse.model_fields
        assert TrialResponse.model_fields["evidence_sufficiency"].default == "unknown"


class TestRunWiresItUp:
    def test_run_attaches_the_verdict(self):
        src = (ROOT / "src" / "agent.py").read_text()
        assert "response.evidence_sufficiency = quality" in src, (
            "run() computes `quality` and must attach it; without this line the grade is "
            "dropped again")

    def test_the_grade_is_recomputed_after_widening(self):
        """Grading only before expansion reports a state the answer is not based on."""
        src = (ROOT / "src" / "agent.py").read_text()
        _, _, tail = src.partition("response.results.extend(sub_response)")
        assert tail, "the sub-question expansion block moved; re-read this"
        before_attach = tail.split("response.evidence_sufficiency")[0]
        assert "self.evaluate_evidence(response.results)" in before_attach, (
            "the verdict must be recomputed after the expansion, not carried from before it")


class TestItReachesTheReport:
    @pytest.fixture
    def agent(self):
        from src.agent import TrialIntelligenceAgent
        with patch.object(TrialIntelligenceAgent, "__init__", lambda self: None):
            a = TrialIntelligenceAgent()
        a.rag = MagicMock()
        return a

    def test_the_header_states_the_sufficiency(self, agent):
        src = (ROOT / "src" / "agent.py").read_text()
        assert '**Evidence sufficiency:** {_sufficiency}' in src
        assert '_sufficiency = getattr(results, "evidence_sufficiency", "unknown")' in src

    def test_insufficient_gets_an_explicit_warning_block(self, agent):
        src = (ROOT / "src" / "agent.py").read_text()
        assert "Evidence was insufficient" in src
        assert "provisional" in src
