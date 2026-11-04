"""Regression tests: no pediatric safety verdict without a cardiac liability model.

`pediatric_safety_assessment` issues a clinical-shaped verdict (`pediatric_safe`,
`recommendation: "Suitable for pediatric development"`) that `report_<run>.json` counts in
`summary.pediatric_safe_count`. Until 2026-11-04 its only `critical` rule -- the only thing that
could set `pediatric_safe=False` -- fired on `herg_ic50 < 10`, and `herg_ic50` came from an
invented expression in `stage_8_ranking`:

    if logp > 3.5 and mw > 350:
        estimated_herg = max(0.5, 30 - (logp * 4) - (mw - 350) * 0.02)

presented as a literature method: "Estimate hERG IC50 from physicochemical properties ... Ref:
Aronov (2005) J Med Chem 48:1341-1348". Aronov 2005 publishes no such equation. Worse, outside
that logP/MW gate no value was set at all, `molecule_data.get("herg_ic50", inf)` returned
infinity, no cardiac flag was raised, and the candidate was reported
`pediatric_safe: true, "Suitable for pediatric development"`. A candidate with dofetilide's
physicochemistry -- MW 441, logP 2.3, hERG IC50 about 10 nM, three orders of magnitude below the
threshold -- passed exactly that way. The function had no test of any kind.

The verdict is now three-valued: absence of a hERG prediction yields `None` (not assessed), never
`True`.
"""
import ast
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from src.pipeline import pediatric_safety_assessment  # noqa: E402

# The dofetilide-like candidate that used to pass: below the old logP/MW gate entirely.
DOFETILIDE_LIKE = {"molecular_weight": 441.0, "logP": 2.3, "tpsa": 75.0, "rotatable_bonds": 8}


class TestNoVerdictWithoutACardiacModel:
    def test_absent_herg_prediction_is_not_assessed_not_safe(self):
        a = pediatric_safety_assessment(DOFETILIDE_LIKE)
        assert a.pediatric_safe is None, "no hERG prediction means no verdict, not 'safe'"
        assert a.cardiac_liability_assessed is False
        assert a.unassessed_reasons
        assert "NOT ASSESSED" in a.recommendation

    def test_the_old_recommendation_string_is_never_issued_unassessed(self):
        a = pediatric_safety_assessment(DOFETILIDE_LIKE)
        assert "Suitable for pediatric development" not in a.recommendation

    def test_an_empty_molecule_is_also_unassessed(self):
        assert pediatric_safety_assessment({}).pediatric_safe is None


class TestVerdictsFromARealPrediction:
    def test_high_herg_probability_raises_the_critical_flag(self):
        a = pediatric_safety_assessment(dict(DOFETILIDE_LIKE,
                                             herg_inhibition_probability=0.93))
        assert a.pediatric_safe is False
        assert a.cardiac_liability_assessed is True
        assert any(f.flag == "HERG_CARDIAC_RISK" and f.severity == "critical"
                   for f in a.safety_flags)

    def test_low_herg_probability_permits_a_positive_verdict(self):
        a = pediatric_safety_assessment(dict(DOFETILIDE_LIKE,
                                             herg_inhibition_probability=0.04))
        assert a.pediatric_safe is True
        assert a.cardiac_liability_assessed is True

    def test_a_measured_ic50_is_still_accepted(self):
        assert pediatric_safety_assessment(
            dict(DOFETILIDE_LIKE, herg_ic50=0.01)).pediatric_safe is False
        assert pediatric_safety_assessment(
            dict(DOFETILIDE_LIKE, herg_ic50=80.0)).pediatric_safe is True

    def test_non_critical_flags_do_not_overturn_an_assessed_verdict(self):
        a = pediatric_safety_assessment({"molecular_weight": 620, "logP": 6.0,
                                         "herg_inhibition_probability": 0.1})
        assert a.pediatric_safe is True
        assert {f.flag for f in a.safety_flags} >= {"HIGH_MW_BBB_RISK", "HIGH_LIPOPHILICITY"}
        assert "additional pediatric safety evaluation" in a.recommendation


class TestTheInventedEquationIsGone:
    @staticmethod
    def _code():
        src = (ENGINE / "src" / "pipeline.py").read_text()
        lines = src.splitlines()
        drop = set()
        for node in ast.walk(ast.parse(src)):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
        return "\n".join(ln for i, ln in enumerate(lines, start=1)
                         if i not in drop and not ln.lstrip().startswith("#"))

    @pytest.mark.parametrize("fragment", [
        "estimated_herg",
        "30 - (logp * 4)",
        "Aronov",
        "48:1341-1348",
    ])
    def test_fragment_is_absent_from_executable_code(self, fragment):
        assert fragment not in self._code(), (
            f"{fragment!r} is back: a hand-written expression must not stand in for a hERG model"
        )

    def test_the_report_counts_the_unassessed_separately(self):
        assert "pediatric_not_assessed_count" in self._code()
