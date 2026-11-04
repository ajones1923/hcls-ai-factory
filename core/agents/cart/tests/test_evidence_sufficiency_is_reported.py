"""CAR-T: an answer must say how much evidence stood behind it, on every path.

Found by HAND on 2026-11-04, after `scripts/sweep_claim_classes.py` had reported cart clean.
The scanner's class B asks whether the verdict vocabulary ever reaches the response surface.
It does, in `app/cart_ui.py`:

    if quality == "insufficient" and plan.sub_questions:

which is the discarding pattern itself, copy-pasted into the UI. The word was present, so
the subject was not flagged. A scanner that reads a second copy of the defect as evidence of
the fix is worse than no scanner; the detector now strips comparison lines out of the
surface text before matching, and catches this file against the pre-fix code.

Two real problems, not one:

  * the API path (`CARTAgent.run`) graded the evidence and dropped the grade, exactly as the
    other seven subjects did;
  * the Streamlit UI *did* display it -- before the sub-question expansion ran. So a user
    saw a grade describing the evidence that triggered widening, not the evidence the answer
    was built from. A stale disclosure reads as a current one.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models import AgentResponse  # noqa: E402

AGENT_SRC = (ROOT / "src" / "agent.py").read_text()
UI_SRC = (ROOT / "app" / "cart_ui.py").read_text()


class TestTheApiPath:
    def test_the_response_model_carries_the_verdict(self):
        assert "evidence_sufficiency" in AgentResponse.model_fields
        assert AgentResponse.model_fields["evidence_sufficiency"].default == "unknown"

    def test_run_attaches_it(self):
        assert "evidence_sufficiency=quality" in AGENT_SRC

    def test_the_grade_is_recomputed_after_widening(self):
        _, _, tail = AGENT_SRC.partition("evidence.hits.extend(sub_evidence.hits)")
        assert tail, "the sub-question expansion block moved; re-read this"
        head = tail.split("evidence_sufficiency")[0]
        assert "self.evaluate_evidence(evidence)" in head, (
            "the verdict must be recomputed after the expansion, not carried from before it")

    def test_insufficient_is_disclosed_in_the_answer(self):
        assert "Evidence was insufficient" in AGENT_SRC
        assert "provisional" in AGENT_SRC


class TestTheUiPath:
    def test_the_ui_regrades_after_expanding(self):
        """It showed a pre-expansion grade, which is a stale disclosure, not a disclosure."""
        _, _, tail = UI_SRC.partition("**Augmented to:**")
        assert tail, "the UI expansion block moved; re-read this"
        assert "agent.evaluate_evidence(evidence)" in tail.split("st.write")[1] + \
               tail[:400], "the UI must re-grade after augmenting the evidence"

    def test_the_ui_shows_the_post_expansion_value(self):
        assert "Evidence quality after expansion" in UI_SRC


class TestTheScannerBlindSpot:
    def test_the_ui_still_contains_the_comparison_that_fooled_the_scanner(self):
        """Not a defect — control flow. Pinned so the next reader knows why B missed cart.

        If this disappears the detector's stripping rule has nothing to strip here, and the
        note in docs/build/CLAIM_CLASSES.md should be revisited rather than silently kept.
        """
        assert 'if quality == "insufficient"' in UI_SRC
