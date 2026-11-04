"""Regression tests: the flagship PDF must report a run, never invent one.

Until 2026-11-04 `generate_vcp_report_enhanced.VCPReportGeneratorEnhanced` emitted a hard-coded
result set under the names of real tools, and both UI export buttons reached it with no data
argument at all:

  * "Top Ranked Drug Candidates" was four literal dicts -- VCP-AI-001..004, docking -8.62 to
    -10.95 kcal/mol, QED/MW/logP and composite scores -- under the prose "NVIDIA BioNeMo MolMIM
    generated novel molecules ... docked against VCP using DiffDock";
  * the variant section attributed a variant (chr9:35,065,263 G>A, rs188935092) and a
    pathogenicity panel (AlphaMissense 0.89, CADD 28.5, ClinVar PATHOGENIC) to "NVIDIA
    Parabricks 4.6 on DGX Spark ... identifying" it;
  * the structure section asserted binding-site evidence: ALA464/GLY479/ASP320/GLY215, pocket
    volume ~450 A^3, druggability score 0.92;
  * the executive summary quoted those docking scores back and printed a "Pipeline Performance"
    table of stage timings nobody measured;
  * the target section credited a search over "3.5 million genomic evidence embeddings".

The file imports no NIM client, no docking code, no RDKit descriptor and no ADMET predictor:
none of it was ever computed. Select target EGFR in the :8505 UI and the exported
EGFR_Drug_Candidate_Report.pdf was the VCP report, with a Parabricks/MolMIM/DiffDock provenance
line and none of the session's own molecules.
"""
import ast
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from generate_vcp_report_enhanced import (  # noqa: E402
    ReportDataUnavailable,
    VCPReportGeneratorEnhanced,
)

# One ranked candidate, the shape run_cloud_nim_report.py writes.
RUN = {
    "pipeline": {
        "target_gene": "EGFR",
        "generator": "CloudMolMIMClient",
        "docking_method": None,
        "num_generated": 1,
        "num_passed_qc": 1,
        "num_docked": 0,
        "seed_compound": "Erlotinib",
        "seed_smiles": "C1=CC2=C(C=C1OCCOC)C(=NC=N2)NC3=CC=CC(=C3)C#C",
    },
    "top_candidates": [{
        "rank": 1,
        "id": "EGFR-TEST-001",
        "smiles": "CCO",
        "composite_score": 0.5,
        "score_type": None,
        "docking_score": None,
        "properties": {"qed": 0.41, "mw": 46.07, "logp": -0.14, "hbd": 1, "hba": 1,
                       "provenance": "hand_written_analogue_library"},
    }],
}


@pytest.fixture
def cwd(tmp_path, monkeypatch):
    """The generator writes an image cache relative to the working directory."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


class TestItRefusesWithoutARun:
    def test_generate_raises_when_no_results_were_supplied(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"))
        with pytest.raises(ReportDataUnavailable):
            gen.generate()
        assert not (cwd / "r.pdf").exists(), "no PDF may be written without a run to report"

    def test_molecules_section_raises_and_says_where_a_real_run_comes_from(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"))
        with pytest.raises(ReportDataUnavailable) as e:
            gen._create_molecules_section()
        assert "run_cloud_nim_report" in str(e.value)

    def test_summary_section_raises_too(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"))
        with pytest.raises(ReportDataUnavailable):
            gen._create_summary_section()

    def test_variant_section_states_its_absence_instead_of_inventing_one(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"))
        text = " ".join(getattr(el, "text", "") for el in gen._create_variant_section())
        assert "NOT INCLUDED" in text
        assert "rs188935092" not in text
        assert "0.89" not in text


class TestItRendersTheRunItIsGiven:
    def test_a_supplied_run_produces_a_pdf(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"), results=RUN)
        out = gen.generate()
        assert Path(out).exists() and Path(out).stat().st_size > 0

    def test_candidate_rows_come_from_the_run_and_carry_provenance(self, cwd):
        gen = VCPReportGeneratorEnhanced(output_path=str(cwd / "r.pdf"), results=RUN)
        rows = [c for el in gen._create_molecules_section()
                for c in getattr(el, "_cellvalues", []) or []]
        flat = [str(cell) for row in rows for cell in row]
        assert "EGFR-TEST-001" in flat
        assert "hand_written_analogue_library" in flat
        # None of the withdrawn candidates may appear for a run that did not contain them.
        assert not any("VCP-AI-00" in c for c in flat)

    def test_an_undocked_candidate_is_never_given_a_kcal_per_mol_number(self):
        assert VCPReportGeneratorEnhanced._dock_cell({}) == "not docked"
        assert VCPReportGeneratorEnhanced._dock_cell(
            {"score_type": "diffdock_confidence", "pose_confidence": 0.42}
        ) == "pose conf. 0.420"
        assert "kcal/mol" in VCPReportGeneratorEnhanced._dock_cell({"docking_score": -9.1})

    def test_a_missing_number_is_reported_as_not_computed_not_defaulted(self):
        assert VCPReportGeneratorEnhanced._num(None, 3) == "not computed"
        assert VCPReportGeneratorEnhanced._num(0.4, 3) == "0.400"


class TestTheFabricatedLiteralsAreGone:
    @staticmethod
    def _code(name):
        """Executable code only: comments and docstrings stripped.

        The removed literals are deliberately quoted in the new comments and docstrings so a
        reader can see what was taken out, so those must not count as occurrences — but string
        *values* must, because the fabrications were table cells.
        """
        src = (ENGINE / name).read_text()
        lines = src.splitlines()
        tree = ast.parse(src)
        drop = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
        return "\n".join(
            ln for i, ln in enumerate(lines, start=1)
            if i not in drop and not ln.lstrip().startswith("#")
        )

    @pytest.mark.parametrize("name", ["generate_vcp_report_enhanced.py", "generate_vcp_report.py"])
    @pytest.mark.parametrize("literal", [
        "VCP-AI-001", "VCP-AI-004",          # the invented candidate IDs
        "-8.62", "-10.95", "0.444",          # the invented docking values and composite score
        "ALA464", "GLY479",                  # invented binding-site residues
        "450 Å³",                  # invented pocket volume
        "Druggability Score",                # invented druggability metric
        "3.5 million",                       # unmeasured embedding count
        "120-240 min", "30-110 minutes",     # unmeasured stage timings
        "35,065,263", "rs188935092",         # the hard-coded demo variant
        "28.5",                              # the invented CADD score
    ])
    def test_literal_is_absent(self, name, literal):
        assert literal not in self._code(name), (
            f"{literal!r} is still emitted by {name}; it was never computed"
        )

    @pytest.mark.parametrize("claim", [
        "MolMIM generated novel molecules",
        "All candidates satisfy Lipinski",
        "favorable ADMET predictions",
        "Parabricks 4.6 on DGX Spark processed",
    ])
    def test_tool_claim_is_absent(self, claim):
        assert claim not in self._code("generate_vcp_report_enhanced.py")

    def test_the_sibling_generator_is_withdrawn_not_merely_unused(self):
        from generate_vcp_report import VCPReportGenerator
        with pytest.raises(RuntimeError) as e:
            VCPReportGenerator()
        assert "withdrawn" in str(e.value)


class TestTheUiPassesItsOwnMolecules:
    def test_export_button_no_longer_constructs_the_generator_without_data(self):
        src = (ENGINE / "app" / "discovery_ui.py").read_text()
        code = "\n".join(ln for ln in src.splitlines() if not ln.lstrip().startswith("#"))
        assert "VCPReportGeneratorEnhanced(output_path=output_path)" not in code
        assert "_session_results(target, molecules)" in code

    def test_the_session_payload_declares_that_nothing_was_docked(self):
        code = (ENGINE / "app" / "discovery_ui.py").read_text()
        assert "'num_docked': 0" in code
        assert "'docking_method': None" in code
