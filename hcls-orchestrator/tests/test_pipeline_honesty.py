"""Regression tests: the orchestrator must not fabricate stages it did not run.

Before 2026-11-03 `run_demo_pipeline()`:
  * fell back to molecules built by string-replacing the seed SMILES, with
    molecular_weight / logP / qed drawn from `random.uniform()`;
  * printed "Stage 7: Molecular Docking -- Running DiffDock..." and then filled
    docking_score from `random.uniform(-12, -6)` and confidence from
    `random.uniform(0.6, 0.95)`, with DiffDock never contacted;
  * hard-coded `prepared: True` on every structure although no preparation happens.

Run with:
    PYTHONPATH=<repo>/lib python -m pytest hcls-orchestrator/tests -q
"""
import importlib.util
from pathlib import Path

import pytest

ORCH = Path(__file__).resolve().parent.parent


def _executable_code(path: Path) -> str:
    """Source with comments AND docstrings stripped.

    The fabrications that were removed are quoted verbatim in the new comments and docstrings so
    a reader can see what went; those must not count as occurrences. String *values* must, since
    several of the fabrications were literals.
    """
    import ast

    src = path.read_text()
    lines = src.splitlines()
    drop = set()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                and isinstance(node.value.value, str):
            drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
    return "\n".join(ln for i, ln in enumerate(lines, start=1)
                     if i not in drop and not ln.lstrip().startswith("#"))


def _load():
    spec = importlib.util.spec_from_file_location("run_pipeline", ORCH / "run_pipeline.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rp():
    return _load()


class TestNoFabricatedDocking:
    def test_missing_structure_file_yields_no_results(self, rp):
        results, status = rp._run_docking(
            [{"id": "m1", "smiles": "CCO"}],
            {"gene": "VCP", "pdb_ids": ["ZZZZ"]},   # no such cached structure
        )
        assert results == [], "docking must return nothing, not placeholder scores"
        assert status["performed"] is False
        assert "no cached structure" in status["reason"]

    def test_status_explains_itself(self, rp):
        _, status = rp._run_docking([], {"gene": "VCP", "pdb_ids": []})
        assert status["reason"], "a skipped stage must record why"
        assert "SKIPPED" in status["message"]

    def test_no_random_fallback_remains_in_the_source(self):
        # Executable lines only: the removed calls are quoted in explanatory comments.
        code = "\n".join(ln for ln in (ORCH / "run_pipeline.py").read_text().splitlines()
                         if not ln.lstrip().startswith("#"))
        # The exact fabrications that were removed. Narrow, literal checks so this test
        # fails if the fallbacks are reintroduced in the same form.
        assert "random.uniform(-12, -6)" not in code
        assert "random.uniform(0.6, 0.95)" not in code
        assert "random.uniform(0.35, 0.55)" not in code
        assert "random.uniform(420, 520)" not in code
        # No `random` call survives anywhere in the orchestrator's executable code.
        assert "random." not in code

    def test_docking_is_not_announced_before_it_is_attempted(self, rp):
        # The stage header is now built from the real status, so "Running DiffDock..." can
        # no longer be printed for a run in which nothing was docked.
        src = (ORCH / "run_pipeline.py").read_text()
        assert 'print_stage(7, "Molecular Docking", "Running DiffDock...")' not in src


class TestStructurePreparationIsNotClaimed:
    def test_find_structure_file_returns_none_for_unknown_id(self, rp):
        assert rp._find_structure_file("ZZZZ") is None

    def test_prepared_flag_is_not_hard_coded_true(self):
        src = (ORCH / "run_pipeline.py").read_text()
        assert "'prepared': True" not in src


class TestGenerationFailsLoudly:
    def test_report_records_the_docking_status(self):
        src = (ORCH / "run_pipeline.py").read_text()
        # report.json must carry the docking provenance so a consumer of the file can tell
        # whether a docking term contributed to the ranking.
        assert "'docking': docking_status" in src

    def test_html_report_no_longer_claims_bionemo(self):
        src = (ORCH / "run_pipeline.py").read_text()
        assert "Powered by NVIDIA BioNeMo" not in src


class TestPortalDoesNotSimulateSuccess:
    """The portal's Launch Pipeline button was a sleep loop ending in balloons."""

    @staticmethod
    def _code():
        return "\n".join(
            ln for ln in (ORCH / "portal" / "app.py").read_text().splitlines()
            if not ln.lstrip().startswith("#")
        )

    def test_no_simulated_completion(self):
        code = self._code()
        assert "Pipeline completed successfully!" not in code
        assert "st.balloons()" not in code

    def test_no_sleep_standing_in_for_work(self):
        code = self._code()
        assert "time.sleep" not in code

    def test_launch_actually_invokes_the_orchestrator(self):
        code = self._code()
        assert "run_pipeline.py" in code
        assert "subprocess.run" in code
        # A non-zero exit must be reported as a failure, not swallowed.
        assert "Pipeline FAILED" in code

    def test_nim_status_is_probed_not_hard_coded(self):
        code = self._code()
        assert '("MolMIM", "localhost:8001", False)' not in code
        assert "Mock Mode" not in code


class TestProvenanceSurvivesToTheReport:
    """The generator labels its output; the report used to drop every label.

    `molecule_generator` sets `properties["provenance"] = "hand_written_analogue_library"`,
    `properties["score_is_placeholder"] = True` and `generation_method="analogue_library"`.
    `_mol_to_dict` kept only id/smiles/target_gene/score/properties and the candidate rows kept
    no provenance field at all, so report.json's `top_candidates` and report.html's "Top Drug
    Candidates" presented a hand-written CB-5083 analogue library ranked by a fixed
    `0.9 - 0.05*i` placeholder, with nothing marking either fact.
    """

    class _Mol:
        name = "VCP-Inh-A1"
        smiles = "CCO"
        target_gene = "VCP"
        score = 0.85
        generation_method = "analogue_library"
        properties = {"molecular_weight": 46.07,
                      "provenance": "hand_written_analogue_library",
                      "score_is_placeholder": True}

    def test_mol_to_dict_carries_the_labels(self, rp):
        d = rp._mol_to_dict(self._Mol())
        assert d["provenance"] == "hand_written_analogue_library"
        assert d["score_is_placeholder"] is True
        assert d["generation_method"] == "analogue_library"

    def test_the_html_report_renders_provenance_and_the_placeholder_warning(self, rp, tmp_path):
        report = {
            "target": {"gene": "VCP"},
            "docking": {"performed": False, "reason": "no structure"},
            "summary": {"molecules_generated": 1, "top_candidates": 1},
            "top_candidates": [{
                "rank": 1, "molecule_id": "VCP-Inh-A1", "qed_score": 0.372,
                "composite_score": 0.611, "provenance": "hand_written_analogue_library",
                "score_is_placeholder": True,
            }],
        }
        out = tmp_path / "report.html"
        rp.generate_html_report(report, out)
        html = out.read_text()
        assert "hand_written_analogue_library" in html
        assert "placeholder rank" in html


class TestQedIsComputedNotDefaulted:
    """`properties.get('qed', 0.4)` printed "0.400" for every candidate, to three decimals.

    `_calculate_properties` never computes `qed`, so the ranking always took the default. With
    no docking the weights are 0.5/0.5, which made half of every composite score that constant.
    """

    def test_qed_is_computed_from_the_smiles(self, rp):
        qed = rp._qed({"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "properties": {}})
        assert qed is not None
        assert qed != 0.4, "0.4 was the hard-coded default, not a computed QED"

    def test_distinct_molecules_get_distinct_qed(self, rp):
        a = rp._qed({"smiles": "CCO", "properties": {}})
        b = rp._qed({"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "properties": {}})
        assert a != b

    def test_an_unparseable_smiles_yields_none_not_a_default(self, rp):
        assert rp._qed({"smiles": "not-a-molecule", "properties": {}}) is None
        assert rp._qed({"smiles": "", "properties": {}}) is None

    def test_the_default_is_gone_from_the_source(self):
        assert "get('qed', 0.4)" not in _executable_code(ORCH / "run_pipeline.py")

    def test_a_missing_qed_cell_says_so(self, rp, tmp_path):
        report = {
            "target": {"gene": "VCP"},
            "docking": {"performed": False, "reason": "no structure"},
            "summary": {},
            "top_candidates": [{"rank": 1, "molecule_id": "m1", "qed_score": None,
                                "composite_score": 0.5}],
        }
        out = tmp_path / "report.html"
        rp.generate_html_report(report, out)
        assert "not computed" in out.read_text()


class TestMockDockingIsNotReportedAsPerformed:
    """NIM_ALLOW_MOCK_FALLBACK=true made report.json say docking was performed by
    "MockDiffDockClient" and report.html show the green banner over random scores."""

    def test_mock_poses_mark_the_stage_as_skipped(self, rp, tmp_path, monkeypatch):
        structures = tmp_path / "structures"
        structures.mkdir()
        (structures / "5ftk.pdb").write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n")
        monkeypatch.setattr(rp, "DATA_DIR", tmp_path)

        engine_src = ORCH.parent / "core/engines/therapeutic-discovery/src"
        import sys
        sys.path.insert(0, str(ORCH.parent / "core/engines/therapeutic-discovery"))
        assert engine_src.exists()
        monkeypatch.setenv("NIM_ALLOW_MOCK_FALLBACK", "true")

        from src.nim_clients import MockDiffDockClient, MockMolMIMClient, NIMServiceManager

        manager = NIMServiceManager()
        manager._molmim = MockMolMIMClient()
        manager._diffdock = MockDiffDockClient()
        monkeypatch.setattr("src.nim_clients.create_nim_clients", lambda *a, **k: manager)

        results, status = rp._run_docking(
            [{"id": "m1", "smiles": "CCO"}], {"gene": "VCP", "pdb_ids": ["5FTK"]})
        assert results == [], "simulated poses must not enter the ranking"
        assert status["performed"] is False
        assert "mock" in status["reason"].lower()

    def test_the_guard_is_wired_in_the_source(self):
        # The honesty helpers had no caller anywhere in the repo.
        assert "contains_mock" in _executable_code(ORCH / "run_pipeline.py")


class TestQcDoesNotPassWhatItDidNotExamine:
    """`properties.get('molecular_weight', 500) < 550` let a molecule with no computed MW take
    the literal default 500 and pass an MW filter that never examined it."""

    def test_the_default_is_gone(self):
        assert "get('molecular_weight', 500)" not in _executable_code(
            ORCH / "run_pipeline.py")

    def test_unexamined_molecules_are_counted_separately(self):
        code = (ORCH / "run_pipeline.py").read_text()
        assert "were NOT " in code and "examined" in code


class TestNextflowModulesDoNotFabricate:
    """The MolMIM/DiffDock processes produced random numbers under those tool names."""

    @staticmethod
    def _dd():
        """Non-comment lines only — the removed fabrications are quoted in the new comments.
        Both `#` (inside the Python process bodies) and ` *` (the Groovy block header)."""
        out = []
        for ln in (ORCH / "modules" / "drug_discovery.nf").read_text().splitlines():
            s = ln.lstrip()
            if s.startswith("#") or s.startswith("*") or s.startswith("//") or s.startswith("/*"):
                continue
            out.append(ln)
        return "\n".join(out)

    def test_no_random_molecule_or_docking_generation(self):
        dd = self._dd()
        assert "random.uniform" not in dd
        assert "random.randint" not in dd

    def test_versions_yml_no_longer_claims_the_tools(self):
        dd = self._dd()
        assert "molmim: 1.0" not in dd
        assert "diffdock: 1.0" not in dd

    def test_unimplemented_stages_fail_rather_than_being_ignored(self):
        dd = self._dd()
        assert "sys.exit(1)" in dd
        # `errorStrategy ... 'ignore'` let a failed stage emit nothing and the workflow
        # continue as though it had succeeded.
        assert "'ignore'" not in dd
