"""Regression tests: simulated docking must be impossible to serve as a measurement.

Audit compute#8 / C7. `NIM_ALLOW_MOCK_FALLBACK=true` -- a documented degradation mode -- swaps in
`MockDiffDockClient`, whose `dock()` returned `base_score = -8.0 + random.gauss(0, 2)` labelled
`score_type: "binding_energy"`, `confidence = 0.7 + random.gauss(0, 0.2)`,
`rmsd = random.uniform(0.5, 3.0)`, `hydrogen_bonds = random.randint(0, 4)` and contacts
`[f"ALA{randint(100,500)}", f"GLY{...}", f"ASP{...}"]` -- byte-for-byte the fabrication whose
removal from modules/drug_discovery.nf and run_pipeline.py is the documented fix, reintroduced
behind a client class. It stamped `_provenance="mock"`, but the ranking stage dropped the flag,
so the report showed a green "Docking: MockDiffDockClient" banner over `-9.17 kcal/mol` cells
with `dock_norm = 0.764` at weight 0.4. The word "mock" survived only inside a class name.

Audit compute#33 / C15 additionally: without RDKit, stage 5 did `mol.passed_qc = True` for every
molecule and stage 8 substituted `qed_score = 0.5` for all of them -- a drug-likeness filter that
examined nothing, reported as having run.

The honesty helpers written for this class -- `is_mock_result` / `contains_mock` / `assert_real`
-- had no caller anywhere in the repo. They have one now.
"""
import ast
import json
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from src.models import PipelineConfig  # noqa: E402
from src.nim_clients import MockDiffDockClient, contains_mock  # noqa: E402
from src.pipeline import DrugDiscoveryPipeline  # noqa: E402

RECEPTOR = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"


class TestTheMockClientDoesNotWearRealUnits:
    def test_no_pose_claims_a_binding_energy(self):
        poses = MockDiffDockClient().dock(protein_pdb=RECEPTOR, ligand_smiles="CCO", num_poses=4)
        assert poses
        for p in poses:
            assert p["score_type"] == "mock_simulated"
            assert "docking_score" not in p, "a simulated number must not wear kcal/mol"
            assert p["_provenance"] == "mock"

    def test_no_pose_invents_binding_residues_or_geometry(self):
        poses = MockDiffDockClient().dock(protein_pdb=RECEPTOR, ligand_smiles="CCO", num_poses=4)
        for p in poses:
            assert "contacts" not in p, "a client that performs no docking knows no contacts"
            assert "rmsd" not in p
            assert "hydrogen_bonds" not in p

    def test_contains_mock_sees_them(self):
        assert contains_mock(
            MockDiffDockClient().dock(protein_pdb=RECEPTOR, ligand_smiles="CCO", num_poses=2))


class TestTheReportBoundaryRefusesMock:
    @staticmethod
    def _mock_manager(monkeypatch):
        from unittest.mock import Mock
        manager = Mock()
        manager.molmim.generate.return_value = [{"smiles": "CCO", "score": 0.9}]
        manager.diffdock = MockDiffDockClient()
        manager.get_available_services.return_value = ["molmim", "diffdock"]
        return manager

    def test_a_mock_pose_makes_the_report_refuse_by_default(
        self, monkeypatch, sample_target, temp_output_dir, structure_cache
    ):
        monkeypatch.delenv("NIM_ALLOW_MOCK_FALLBACK", raising=False)
        config = PipelineConfig(target_gene="VCP", num_molecules=1,
                                output_dir=str(temp_output_dir),
                                structure_cache=str(structure_cache))
        pipeline = DrugDiscoveryPipeline(config, nim_manager=self._mock_manager(monkeypatch))
        with pytest.raises(RuntimeError, match="MOCK"):
            pipeline.run_pipeline(sample_target)
        assert not list(temp_output_dir.glob("report_*.json")), (
            "no report may be written from simulated docking without an explicit opt-in"
        )

    def test_with_the_opt_in_the_artifact_is_stamped_and_renamed(
        self, monkeypatch, sample_target, temp_output_dir, structure_cache
    ):
        monkeypatch.setenv("NIM_ALLOW_MOCK_FALLBACK", "true")
        config = PipelineConfig(target_gene="VCP", num_molecules=1,
                                output_dir=str(temp_output_dir),
                                structure_cache=str(structure_cache))
        pipeline = DrugDiscoveryPipeline(config, nim_manager=self._mock_manager(monkeypatch))
        pipeline.run_pipeline(sample_target)

        assert not list(temp_output_dir.glob("report_*.json"))
        mock_reports = list(temp_output_dir.glob("MOCK_report_*.json"))
        assert len(mock_reports) == 1
        report = json.loads(mock_reports[0].read_text())
        assert report["data_provenance"].startswith("MOCK")
        assert "SIMULATED" in report["data_provenance"]
        assert all(c["docking_provenance"] == "mock" for c in report["top_candidates"])

    def test_a_simulated_pose_earns_no_ranking_term(
        self, monkeypatch, sample_target, temp_output_dir, structure_cache
    ):
        monkeypatch.setenv("NIM_ALLOW_MOCK_FALLBACK", "true")
        config = PipelineConfig(target_gene="VCP", num_molecules=1,
                                output_dir=str(temp_output_dir),
                                structure_cache=str(structure_cache),
                                docking_weight=0.4, generation_weight=0.3, qed_weight=0.3)
        pipeline = DrugDiscoveryPipeline(config, nim_manager=self._mock_manager(monkeypatch))
        pipeline.run_pipeline(sample_target)
        c = pipeline.ranked_candidates[0]
        # 0.3*gen + 0.4*0 + 0.3*qed -- the mock contributes nothing.
        expected = 0.3 * c.generation_score + 0.3 * (c.qed_score or 0)
        assert c.composite_score == pytest.approx(round(expected, 4), abs=1e-4)


class TestQcCannotPassWhatItDidNotExamine:
    def test_stage_5_refuses_without_rdkit(
        self, monkeypatch, sample_config, mock_nim_manager, sample_target, temp_output_dir
    ):
        """`mol.passed_qc = True` for every molecule was a pass from a filter that never ran."""
        import builtins

        sample_config.output_dir = str(temp_output_dir)
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()
        pipeline.stage_3_structure_prep()
        pipeline.stage_4_molecule_generation()

        real_import = builtins.__import__

        def no_rdkit(name, *args, **kwargs):
            if name.startswith("rdkit"):
                raise ImportError("simulated: RDKit native extension unavailable")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", no_rdkit)
        with pytest.raises(RuntimeError, match="RDKit"):
            pipeline.stage_5_chemistry_qc()
        assert all(not m.passed_qc for m in pipeline.molecules)

    def test_a_missing_qed_is_not_replaced_by_a_constant(self):
        src = (ENGINE / "src" / "pipeline.py").read_text()
        lines = src.splitlines()
        drop = set()
        for node in ast.walk(ast.parse(src)):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
        code = "\n".join(ln for i, ln in enumerate(lines, start=1)
                         if i not in drop and not ln.lstrip().startswith("#"))
        assert "if mol.properties else 0.5" not in code
        assert "mol.passed_qc = True" not in code.replace(
            "mol.passed_qc = True\n                passed.append(mol)", "")
