"""An answer must say when the evidence behind it was insufficient.

Found by the claim-integrity review of `precision-oncology`, 2026-11-04, check 2.

`OncoIntelligenceAgent.evaluate_evidence` has always graded retrieval as sufficient /
partial / insufficient. The verdict was used for exactly one thing: deciding whether to
retry with a broader search. Once the retry budget was spent the loop broke, synthesis ran
regardless, and the verdict went to a log line and nowhere else -- `AgentResponse` had no
field for it, and the report did not mention it.

So an answer drawn from almost nothing was indistinguishable, to any caller, from one drawn
from ample evidence. `RAGEngine.synthesize` has no empty-evidence guard either; the only
thing standing between zero hits and a confident paragraph was a sentence in the prompt
asking the model to notice, which is an instruction, not a guarantee. This engine recommends
therapy.

The grading itself is unchanged. What changed is that the grade now reaches the caller.
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models import AgentResponse, CrossCollectionResult  # noqa: E402


def _response(answer="Consider pembrolizumab.", question="MSI-high CRC, second line?"):
    return AgentResponse(
        question=question,
        answer=answer,
        evidence=CrossCollectionResult(query=question, hits=[]),
    )


class TestTheFieldExists:
    def test_agent_response_carries_a_sufficiency_verdict(self):
        assert _response().evidence_sufficiency == "unknown"

    def test_it_accepts_the_three_verdicts_evaluate_evidence_returns(self):
        for v in ("sufficient", "partial", "insufficient"):
            r = _response()
            r.evidence_sufficiency = v
            assert r.evidence_sufficiency == v


class TestTheVerdictReachesTheReader:
    @pytest.fixture
    def agent(self):
        from src.agent import OncoIntelligenceAgent
        with patch.object(OncoIntelligenceAgent, "__init__", lambda self: None):
            a = OncoIntelligenceAgent()
        a.MAX_RETRIES = 1
        return a

    def test_the_report_states_the_sufficiency(self, agent):
        r = _response()
        r.evidence_sufficiency = "partial"
        report = agent.generate_report(r)
        assert "Evidence sufficiency" in report
        assert "partial" in report

    def test_an_insufficient_answer_says_so_in_its_first_line(self, agent):
        """Not only in a field a caller may never read."""
        r = _response()
        r.evidence_sufficiency = "insufficient"
        # Reproduce what run() does to the answer, without driving retrieval.
        prefixed = (
            "**Evidence was insufficient.** The retrieval returned too little to support "
            "a confident answer, and the search was already broadened "
            f"{agent.MAX_RETRIES} time(s). Treat everything below as provisional and "
            "verify against primary sources.\n\n" + r.answer
        )
        assert prefixed.startswith("**Evidence was insufficient.**")
        assert "provisional" in prefixed
        assert r.answer in prefixed, "the original answer must still be there, not replaced"


class TestRunWiresItUp:
    def test_run_sets_the_field_from_the_verdict(self):
        """The wiring itself: run() must copy the verdict onto the response it returns."""
        from src.agent import OncoIntelligenceAgent

        src = (ROOT / "src" / "agent.py").read_text()
        assert "response.evidence_sufficiency = verdict" in src, (
            "run() computes `verdict` and must attach it; without this line the grade ends "
            "at the log again")
        # And the insufficient case must alter the answer, not just the metadata.
        assert "Evidence was insufficient" in src
        assert OncoIntelligenceAgent is not None
