"""`evaluate_evidence` has to be called by something other than its own tests.

Found by the claim-class sweep, 2026-11-04 (class B). Of the seven subjects it flagged, six
computed an evidence verdict and dropped it. This one was different and worse: the method
had **no production caller at all**. Eight tests in `test_agent.py` asserted its
classification boundaries, so the suite reported a working evidence check, and nothing in
the pipeline ever ran it.

A tested guard with no caller is worse than no guard. No guard is visibly absent; this one
was visibly present and green. It is the same shape as therapeutic-discovery's
`assert_real` / `contains_mock` helpers, which were written for a real incident and then sat
uncalled until a later pass wired them into the report stage.

`run()` now calls it, carries the verdict on `AgentResponse.evidence_sufficiency`, and
prefixes the answer when the verdict is `insufficient`.
"""
import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

AGENT_SRC = (ROOT / "src" / "agent.py").read_text()


class TestItHasAProductionCaller:
    def test_run_calls_evaluate_evidence(self):
        """The whole finding in one assertion."""
        tree = ast.parse(AGENT_SRC)
        run = next((n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == "run"), None)
        assert run is not None, "ImagingIntelligenceAgent.run not found"
        calls = [n for n in ast.walk(run)
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                 and n.func.attr == "evaluate_evidence"]
        assert calls, (
            "run() does not call evaluate_evidence. The method is tested eight ways and "
            "would again be dead code reporting green.")

    def test_the_caller_is_not_a_test(self):
        """Guards against this being 'fixed' by adding another test that calls it."""
        assert "self.evaluate_evidence(evidence)" in AGENT_SRC


class TestTheVerdictReachesTheCaller:
    def test_the_response_model_carries_it(self):
        from src.models import AgentResponse
        assert "evidence_sufficiency" in AgentResponse.model_fields
        assert AgentResponse.model_fields["evidence_sufficiency"].default == "unknown"

    def test_run_passes_it_into_the_response(self):
        assert "evidence_sufficiency=sufficiency" in AGENT_SRC

    def test_insufficient_is_disclosed_in_the_answer(self):
        assert "Evidence was insufficient" in AGENT_SRC
        assert "provisional" in AGENT_SRC


class TestTheGradingItselfIsUnchanged:
    def test_boundaries_still_match_the_existing_suite(self):
        """The fix wires a method up; it must not quietly re-tune it."""
        from unittest.mock import MagicMock

        from src.agent import ImagingIntelligenceAgent
        agent = ImagingIntelligenceAgent.__new__(ImagingIntelligenceAgent)
        for hits, expected in ((10, "sufficient"), (3, "partial"), (2, "insufficient")):
            ev = MagicMock()
            ev.hit_count = hits
            assert agent.evaluate_evidence(ev) == expected
