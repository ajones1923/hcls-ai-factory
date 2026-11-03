"""Regression tests: a DiffDock confidence must never be handled as a binding affinity.

DiffDock's confidence model is a binary classifier trained to predict whether a pose is
within 2 A RMSD of the true pose, so its output is a log-odds where HIGHER is better
(Corso et al., ICLR 2023). Until 2026-11-03 the cloud client wrote that value into
`docking_score` and the ranking stage normalized it with `-score / 12.0`, which both
inverted the sign and asserted a kcal/mol unit the number never had.
"""
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import DockingResult, RankedCandidate  # noqa: E402
from src.nim_clients import CloudDiffDockClient  # noqa: E402


def _parse(payload):
    client = CloudDiffDockClient.__new__(CloudDiffDockClient)  # no network/api key needed
    return CloudDiffDockClient._parse_cloud_response(client, payload, "CCO")


class TestCloudResponseUnits:
    """`_parse_cloud_response` must report a confidence as a confidence."""

    PAYLOAD = {
        "position_confidence": [-4.46, -3.77, -5.05],
        "ligand_positions": ["sdf-a", "sdf-b", "sdf-c"],
    }

    def test_confidence_is_not_written_into_docking_score(self):
        for pose in _parse(self.PAYLOAD):
            assert pose["docking_score"] is None, (
                "DiffDock predicts no binding energy; docking_score must stay None"
            )
            assert pose["score_type"] == "diffdock_confidence"

    def test_raw_confidence_is_preserved_and_calibrated(self):
        poses = _parse(self.PAYLOAD)
        raw = {p["confidence_score"] for p in poses}
        assert raw == {-4.46, -3.77, -5.05}
        for p in poses:
            # sigmoid(c) = P(pose within 2 A RMSD)
            assert abs(p["confidence"] - 1 / (1 + math.exp(-p["confidence_score"]))) < 1e-12
            assert 0.0 <= p["confidence"] <= 1.0

    def test_best_pose_first_means_highest_confidence(self):
        poses = _parse(self.PAYLOAD)
        # -3.77 is the most confident of the three; the old code sorted ascending by the
        # value stored in docking_score and so returned the WORST pose as poses[0].
        assert poses[0]["confidence_score"] == -3.77
        confidences = [p["confidence_score"] for p in poses]
        assert confidences == sorted(confidences, reverse=True)

    def test_rmsd_is_not_claimed_to_be_zero(self):
        # The API returns no RMSD; reporting 0.0 asserted a perfect match to a reference pose.
        assert all(p["rmsd"] is None for p in _parse(self.PAYLOAD))

    def test_missing_confidence_does_not_become_a_score(self):
        poses = _parse({"ligand_positions": ["sdf-a"]})
        assert poses[0]["confidence_score"] is None
        assert poses[0]["docking_score"] is None


class TestDockingResultModel:
    """The model must be able to say 'no affinity was predicted'."""

    def test_docking_score_may_be_absent(self):
        r = DockingResult(molecule_id="m1", structure_id="5FTK",
                          score_type="diffdock_confidence", confidence_score=-3.77,
                          confidence=0.023)
        assert r.docking_score is None

    def test_binding_energy_result_still_works(self):
        r = DockingResult(molecule_id="m1", structure_id="5FTK", docking_score=-8.5)
        assert r.score_type == "binding_energy"
        assert r.confidence_score is None


class TestRankingNormalization:
    """Stage 8's docking term must branch on the units it was given."""

    @staticmethod
    def _normalize(dock_result):
        """The normalization logic of pipeline.stage_8_ranking, isolated."""
        if dock_result is None:
            return 0.0
        if dock_result.score_type == "diffdock_confidence":
            return max(0.0, min(1.0, dock_result.confidence))
        if dock_result.docking_score is not None:
            return max(0.0, min(1.0, -dock_result.docking_score / 12.0))
        return 0.0

    def test_binding_energy_normalizes_on_kcal_per_mol(self):
        r = DockingResult(molecule_id="m", structure_id="s", docking_score=-6.0)
        assert self._normalize(r) == 0.5

    def test_confident_pose_is_not_scored_as_zero(self):
        # A DiffDock confidence of +2.0 is a GOOD pose. The old `-score / 12` clamped every
        # positive confidence to 0.0 -- the most confident poses ranked worst.
        conf = 1 / (1 + math.exp(-2.0))
        r = DockingResult(molecule_id="m", structure_id="s",
                          score_type="diffdock_confidence", confidence_score=2.0,
                          confidence=conf)
        assert self._normalize(r) > 0.8
        assert max(0.0, min(1.0, -2.0 / 12.0)) == 0.0  # what the old code produced

    def test_confidence_ordering_survives_normalization(self):
        def res(c):
            return DockingResult(molecule_id="m", structure_id="s",
                                 score_type="diffdock_confidence", confidence_score=c,
                                 confidence=1 / (1 + math.exp(-c)))
        assert self._normalize(res(-3.77)) > self._normalize(res(-5.05))

    def test_no_docking_result_scores_zero(self):
        assert self._normalize(None) == 0.0

    def test_ranked_candidate_carries_the_units(self):
        c = RankedCandidate(rank=1, molecule_id="m", smiles="CCO", target_gene="VCP",
                            score_type="diffdock_confidence", confidence_score=-3.77,
                            generation_score=0.8, composite_score=0.5)
        assert c.docking_score is None
        assert c.score_type == "diffdock_confidence"
