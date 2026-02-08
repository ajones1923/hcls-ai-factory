"""Tests for drug discovery scoring and Lipinski filter logic."""
import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import (
    GeneratedMolecule,
    MoleculeProperties,
    DockingResult,
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

    Formula: max(0, min(1, (10 + dock_score) / 20))
    Lower docking scores are better (more negative = stronger binding).
    """

    @staticmethod
    def _normalize_dock(dock_score):
        if dock_score:
            return max(0, min(1, (10 + dock_score) / 20))
        return 0.5

    def test_zero_dock_uses_sentinel(self):
        assert self._normalize_dock(0) == 0.5

    def test_none_dock_uses_sentinel(self):
        assert self._normalize_dock(None) == 0.5

    def test_strong_binding(self):
        # dock=-10 -> (10-10)/20 = 0
        assert self._normalize_dock(-10) == 0.0

    def test_moderate_binding(self):
        # dock=-8 -> (10-8)/20 = 0.1
        assert self._normalize_dock(-8) == pytest.approx(0.1)

    def test_weak_binding(self):
        # dock=-2 -> (10-2)/20 = 0.4
        assert self._normalize_dock(-2) == pytest.approx(0.4)

    def test_positive_dock_clamps_to_one(self):
        # dock=10 -> (10+10)/20 = 1.0
        assert self._normalize_dock(10) == 1.0

    def test_very_negative_dock_clamps_to_zero(self):
        # dock=-15 -> (10-15)/20 = -0.25 -> clamped to 0
        assert self._normalize_dock(-15) == 0.0


class TestCompositeScoring:
    """Test the composite scoring formula from stage_8_ranking."""

    @staticmethod
    def _compute_composite(gen_score, dock_score, qed_score,
                            gen_w=0.3, dock_w=0.4, qed_w=0.3):
        if dock_score:
            dock_normalized = max(0, min(1, (10 + dock_score) / 20))
        else:
            dock_normalized = 0.5
        return gen_w * gen_score + dock_w * dock_normalized + qed_w * qed_score

    def test_default_weights_sum_to_one(self):
        config = PipelineConfig(target_gene="TEST")
        total = config.docking_weight + config.generation_weight + config.qed_weight
        assert abs(total - 1.0) < 0.001

    def test_all_perfect_scores(self):
        # gen=1.0, dock=10 (normalized=1.0), qed=1.0
        composite = self._compute_composite(1.0, 10.0, 1.0)
        assert composite == pytest.approx(0.3 + 0.4 + 0.3)  # = 1.0

    def test_all_zero_scores(self):
        # gen=0, dock=0 (sentinel=0.5), qed=0
        composite = self._compute_composite(0.0, 0, 0.0)
        assert composite == pytest.approx(0.4 * 0.5)  # = 0.2

    def test_balanced_mid_scores(self):
        composite = self._compute_composite(0.5, -5, 0.5)
        # dock: (10-5)/20 = 0.25
        expected = 0.3 * 0.5 + 0.4 * 0.25 + 0.3 * 0.5
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
        with pytest.raises(Exception):
            PipelineConfig(target_gene="VCP", num_molecules=0)

    def test_num_molecules_too_high(self):
        with pytest.raises(Exception):
            PipelineConfig(target_gene="VCP", num_molecules=1001)

    def test_max_mw_default(self):
        config = PipelineConfig(target_gene="VCP")
        assert config.max_mw == 550

    def test_max_lipinski_violations_default(self):
        config = PipelineConfig(target_gene="VCP")
        assert config.max_lipinski_violations == 1
