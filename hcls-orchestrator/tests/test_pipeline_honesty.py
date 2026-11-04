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
