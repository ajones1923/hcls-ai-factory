"""Tests for drug discovery scoring and Lipinski filter logic."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import (
    DockingResult,
    GeneratedMolecule,
    MoleculeProperties,
    PipelineConfig,
    RankedCandidate,
)


class TestLipinskiViolations:
    """Test Lipinski Rule of Five violation counting.

    Rules: MW > 500, logP > 5, HBD > 5, HBA > 10 each add one violation.
    This mirrors the logic in pipeline.py stage_5_chemistry_qc().
    """

    @staticmethod
    def _count_violations(mw=400, logp=3.0, hbd=2, hba=5):
        """Count violations using the same logic as stage_5."""
        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        return violations

    def test_no_violations(self):
        assert self._count_violations(mw=400, logp=3, hbd=2, hba=5) == 0

    def test_mw_boundary_pass(self):
        assert self._count_violations(mw=500) == 0

    def test_mw_boundary_fail(self):
        assert self._count_violations(mw=501) == 1

    def test_logp_boundary_pass(self):
        assert self._count_violations(logp=5.0) == 0

    def test_logp_boundary_fail(self):
        assert self._count_violations(logp=5.1) == 1

    def test_hbd_boundary_pass(self):
        assert self._count_violations(hbd=5) == 0

    def test_hbd_boundary_fail(self):
        assert self._count_violations(hbd=6) == 1

    def test_hba_boundary_pass(self):
        assert self._count_violations(hba=10) == 0

    def test_hba_boundary_fail(self):
        assert self._count_violations(hba=11) == 1

    def test_all_violations(self):
        assert self._count_violations(mw=600, logp=6, hbd=7, hba=12) == 4


class TestDockNormalization:
    """Test docking score normalization formula.

    Formula: max(0.0, min(1.0, -dock_score / 12.0))
    Lower docking scores are better (more negative = stronger binding).
    Range: -12 kcal/mol (excellent) → 1.0, 0 kcal/mol (no binding) → 0.0
    When dock_result is None → 0.0 (no docking data = worst score).
    """

    @staticmethod
    def _normalize_dock(dock_score, has_result=True):
        """Mirror the normalization logic from pipeline.py stage_8_ranking."""
        if has_result:
            return max(0.0, min(1.0, -dock_score / 12.0))
        return 0.0

    def test_zero_dock_no_binding(self):
        # dock=0 with a result → -0/12 = 0.0 (no binding energy)
        assert self._normalize_dock(0) == 0.0

    def test_none_dock_no_result(self):
        # No docking result → 0.0 (worst score)
        assert self._normalize_dock(0, has_result=False) == 0.0

    def test_excellent_binding(self):
        # dock=-12 -> -(-12)/12 = 1.0
        assert self._normalize_dock(-12) == 1.0

    def test_strong_binding(self):
        # dock=-8 -> -(-8)/12 = 0.6667
        assert self._normalize_dock(-8) == pytest.approx(2 / 3)

    def test_moderate_binding(self):
        # dock=-6 -> -(-6)/12 = 0.5
        assert self._normalize_dock(-6) == pytest.approx(0.5)

    def test_weak_binding(self):
        # dock=-2 -> -(-2)/12 = 0.1667
        assert self._normalize_dock(-2) == pytest.approx(1 / 6)

    def test_positive_dock_clamps_to_zero(self):
        # dock=+2 -> -(+2)/12 = -0.167 -> clamped to 0
        assert self._normalize_dock(2) == 0.0

    def test_very_negative_dock_clamps_to_one(self):
        # dock=-15 -> -(-15)/12 = 1.25 -> clamped to 1.0
        assert self._normalize_dock(-15) == 1.0


class TestCompositeScoring:
    """Test the composite scoring formula from stage_8_ranking.

    composite = gen_w * gen_score + dock_w * dock_normalized + qed_w * qed_score
    dock_normalized = max(0.0, min(1.0, -dock_score / 12.0)) when has result, else 0.0
    """

    @staticmethod
    def _compute_composite(gen_score, dock_score, qed_score,
                            has_dock_result=True,
                            gen_w=0.3, dock_w=0.4, qed_w=0.3):
        if has_dock_result:
            dock_normalized = max(0.0, min(1.0, -dock_score / 12.0))
        else:
            dock_normalized = 0.0
        return gen_w * gen_score + dock_w * dock_normalized + qed_w * qed_score

    def test_default_weights_sum_to_one(self):
        config = PipelineConfig(target_gene="TEST")
        total = config.docking_weight + config.generation_weight + config.qed_weight
        assert abs(total - 1.0) < 0.001

    def test_all_perfect_scores(self):
        # gen=1.0, dock=-12 (normalized=1.0), qed=1.0
        composite = self._compute_composite(1.0, -12.0, 1.0)
        assert composite == pytest.approx(1.0)

    def test_all_zero_scores_no_dock(self):
        # gen=0, no dock result (normalized=0.0), qed=0
        composite = self._compute_composite(0.0, 0, 0.0, has_dock_result=False)
        assert composite == pytest.approx(0.0)

    def test_zero_dock_with_result(self):
        # gen=0, dock=0 with result (normalized=0.0), qed=0
        composite = self._compute_composite(0.0, 0, 0.0, has_dock_result=True)
        assert composite == pytest.approx(0.0)

    def test_balanced_mid_scores(self):
        composite = self._compute_composite(0.5, -6, 0.5)
        # dock: -(-6)/12 = 0.5
        expected = 0.3 * 0.5 + 0.4 * 0.5 + 0.3 * 0.5
        assert composite == pytest.approx(expected)

    def test_strong_docking_dominates(self):
        # Strong docking (-10) should produce high dock component
        composite = self._compute_composite(0.5, -10, 0.5)
        dock_norm = 10 / 12.0  # 0.8333
        expected = 0.3 * 0.5 + 0.4 * dock_norm + 0.3 * 0.5
        assert composite == pytest.approx(expected)


class TestPipelineConfigWeights:
    """Test PipelineConfig weight validation."""

    def test_default_weights(self):
        config = PipelineConfig(target_gene="VCP")
        assert config.docking_weight == 0.4
        assert config.generation_weight == 0.3
        assert config.qed_weight == 0.3

    def test_num_molecules_range(self):
        config = PipelineConfig(target_gene="VCP", num_molecules=1)
        assert config.num_molecules == 1

        config = PipelineConfig(target_gene="VCP", num_molecules=1000)
        assert config.num_molecules == 1000

    def test_num_molecules_too_low(self):
        with pytest.raises(ValueError):
            PipelineConfig(target_gene="VCP", num_molecules=0)

    def test_num_molecules_too_high(self):
        with pytest.raises(ValueError):
            PipelineConfig(target_gene="VCP", num_molecules=1001)

    def test_max_mw_default(self):
        config = PipelineConfig(target_gene="VCP")
        assert config.max_mw == 550

    def test_max_lipinski_violations_default(self):
        config = PipelineConfig(target_gene="VCP")
        assert config.max_lipinski_violations == 1
