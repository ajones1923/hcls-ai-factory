"""Cardiology intelligence: an answer must say how much evidence stood behind it.

Found by the claim-class sweep, 2026-11-04 (class B). Seven subjects shared the defect:
`evaluate_evidence` graded retrieval `sufficient` / `partial` / `insufficient`, the grade
decided whether to widen the search, and was then dropped -- so an answer built on almost
nothing was indistinguishable from one built on good evidence.

This engine produces cardiology decision support, including risk stratification.

Two fixes, the second not cosmetic: the verdict is carried to the caller, and the grade is
RECOMPUTED after the sub-question expansion. The original call graded the evidence that
triggered widening, which is not the evidence the answer ends up resting on.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.models import CardioResponse  # noqa: E402

SRC = (ROOT / "src" / "agent.py").read_text()


class TestTheFieldExists:
    def test_the_response_model_carries_the_verdict(self):
        assert "evidence_sufficiency" in CardioResponse.model_fields
        assert CardioResponse.model_fields["evidence_sufficiency"].default == "unknown"


class TestRunWiresItUp:
    def test_run_attaches_the_verdict(self):
        assert "response.evidence_sufficiency = quality" in SRC, (
            "run() computes the grade and must attach it; without this the grade is dropped")

    def test_the_grade_is_recomputed_after_widening(self):
        _, _, tail = SRC.partition("response.results.extend(sub_response)")
        assert tail, "the sub-question expansion block moved; re-read this"
        head = tail.split("evidence_sufficiency")[0]
        assert "self._evaluate_evidence_quality(response.results)" in head, (
            "the verdict must be recomputed after the expansion, not carried from before it")


class TestTheReaderIsTold:
    def test_an_insufficient_verdict_is_disclosed(self):
        assert "Evidence was insufficient" in SRC
        assert "provisional" in SRC
