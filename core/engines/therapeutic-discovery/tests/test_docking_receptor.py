"""Regression tests: stage 3 must not claim preparation, stage 7 must not dock an accession.

Audit compute#10 / C4: `stage_3_structure_prep` was a four-line comment ("In production, would:
1. Download PDB files 2. Remove water, ions 3. Add hydrogens 4. Identify binding site") followed
by `for struct in self.structures.structures: struct.prepared = True`. `report_<run>.json` then
carried `"prepared": true` for structures that had never been downloaded, let alone stripped of
water or protonated, and a consumer read that as docking-ready input. The identical bug was fixed
in hcls-orchestrator/run_pipeline.py, whose honesty test greps only that file -- so this copy
survived.

Audit compute#11 / C5: `stage_7_docking` passed `self.structures.primary_structure` -- a
four-character accession, "5FTK" -- as `protein_pdb`. The local client only substituted file
content when `Path(protein_pdb).exists()`, so the literal string was POSTed as the receptor; the
cloud client uploaded those four bytes to NVCF as `contentType="chemical/x-pdb"`. Either way
stage 8 weighted the result at 0.4 and stage 9 wrote it up as docking "against 5FTK".
"""
import ast
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from src.models import PipelineConfig  # noqa: E402
from src.nim_clients import resolve_protein_pdb  # noqa: E402
from src.pipeline import DrugDiscoveryPipeline  # noqa: E402


class TestStage3DoesNotClaimPreparation:
    def test_prepared_is_false_and_says_what_happened(
        self, sample_config, mock_nim_manager, sample_target, temp_output_dir
    ):
        sample_config.output_dir = str(temp_output_dir)
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()
        pipeline.stage_3_structure_prep()

        assert pipeline.structures.structures
        for struct in pipeline.structures.structures:
            assert struct.prepared is False
            assert "not performed" in (struct.preparation or "")

    def test_the_report_cannot_say_prepared_true(
        self, sample_config, mock_nim_manager, sample_target, temp_output_dir
    ):
        import json

        sample_config.output_dir = str(temp_output_dir)
        sample_config.num_molecules = 2
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        pipeline.run_pipeline(sample_target)
        report = json.loads(next(temp_output_dir.glob("report_*.json")).read_text())
        assert report["structures_used"]
        assert all(s["prepared"] is False for s in report["structures_used"])

    def test_the_hard_coded_assignment_is_gone_from_the_source(self):
        src = (ENGINE / "src" / "pipeline.py").read_text()
        lines = src.splitlines()
        drop = set()
        for node in ast.walk(ast.parse(src)):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
        code = "\n".join(ln for i, ln in enumerate(lines, start=1)
                         if i not in drop and not ln.lstrip().startswith("#"))
        assert "struct.prepared = True" not in code


class TestStage3LocatesCoordinates:
    def test_coordinates_found_in_the_cache_are_recorded(
        self, sample_config, mock_nim_manager, sample_target, temp_output_dir
    ):
        sample_config.output_dir = str(temp_output_dir)
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()
        pipeline.stage_3_structure_prep()
        assert any(s.structure_file for s in pipeline.structures.structures)

    def test_pdb_cache_dir_env_is_searched(
        self, sample_config, mock_nim_manager, temp_output_dir, monkeypatch, tmp_path
    ):
        """The engine's own StructureViewer caches under data/structures/pdb_cache, which this
        pipeline never consulted, and the project cache lives outside the repo entirely."""
        external = tmp_path / "external_cache"
        external.mkdir()
        (external / "1abc.pdb").write_text("ATOM      1  CA  ALA A   1   0.0 0.0 0.0\n")
        monkeypatch.setenv("PDB_CACHE_DIR", str(external))
        sample_config.output_dir = str(temp_output_dir)
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        assert pipeline._resolve_structure_file("1ABC") == external / "1abc.pdb"


class TestStage7RefusesAnAccession:
    def test_docking_without_coordinates_raises(
        self, mock_nim_manager, sample_target, temp_output_dir
    ):
        config = PipelineConfig(
            target_gene="VCP",
            num_molecules=1,
            output_dir=str(temp_output_dir),
            structure_cache=str(temp_output_dir / "empty"),   # nothing cached
        )
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim_manager)
        pipeline.stage_0_initialize(sample_target)
        pipeline.stage_1_normalize_target()
        pipeline.stage_2_structure_discovery()
        pipeline.stage_3_structure_prep()
        pipeline.stage_4_molecule_generation()
        pipeline.stage_5_chemistry_qc()

        with pytest.raises(RuntimeError) as e:
            pipeline.stage_7_docking()
        assert "accession" in str(e.value)
        assert pipeline.docking_results == [], "no pose may be recorded for a refused stage"

    def test_the_receptor_passed_to_the_client_is_coordinates_not_an_id(
        self, sample_config, mock_nim_manager, sample_target, temp_output_dir
    ):
        sample_config.output_dir = str(temp_output_dir)
        sample_config.num_molecules = 1
        pipeline = DrugDiscoveryPipeline(sample_config, nim_manager=mock_nim_manager)
        pipeline.run_pipeline(sample_target)

        kwargs = mock_nim_manager.diffdock.dock.call_args.kwargs
        receptor = kwargs["protein_pdb"]
        assert receptor not in {s.pdb_id for s in pipeline.structures.structures}
        assert "ATOM" in receptor


class TestTheClientsRefuseANonStructure:
    @pytest.mark.parametrize("bogus", ["5FTK", "PDB:5FTK", "", "   ", "not a structure"])
    def test_rejects_a_pdb_accession_as_a_receptor(self, bogus):
        with pytest.raises(ValueError):
            resolve_protein_pdb(bogus)

    def test_accepts_coordinate_text(self):
        text = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n"
        assert resolve_protein_pdb(text) == text

    def test_accepts_a_path_and_reads_it(self, tmp_path):
        f = tmp_path / "r.pdb"
        f.write_text("HETATM    1  ZN   ZN A   1       0.000   0.000   0.000\n")
        assert "HETATM" in resolve_protein_pdb(str(f))

    def test_both_clients_route_through_the_guard(self):
        src = (ENGINE / "src" / "nim_clients.py").read_text()
        # One call in DiffDockClient.dock, one in CloudDiffDockClient.dock.
        assert src.count("resolve_protein_pdb(protein_pdb)") >= 2
        assert "if Path(protein_pdb).exists():" not in src, (
            "an existence check alone let a 4-character accession through"
        )
