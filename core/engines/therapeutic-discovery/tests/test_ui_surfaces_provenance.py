"""Regression tests: the :8505 surface must not present a hand-written library as model output.

Audit crosscut I1. `therapeutic-discovery-engine` is registered `live` on :8505, the page
subtitle read "Target Hypothesis -> Structural Evidence -> BioNeMo Molecule Generation" and the
step heading "Step 3: BioNeMo Molecule Generation / AI-generated drug candidates" -- but
`discovery_ui.py` constructed `MoleculeGenerator()` (use_bionemo=True, require_real=False),
`_check_bionemo_available()` returns False on this box, and the call fell through to
`_get_vcp_inhibitor_analogues()`: a hard-coded nine-entry SMILES list with `score = 0.85 - i*0.05`,
a fabricated descending rank.

Two separate defects, both fixed here:

  1. The library was appended for ANY `target_gene`, each molecule stamped with the gene the user
     asked for. Request EGFR and you received VCP-Inh-A1..H1 labelled target_gene="EGFR".
  2. The library *is* honest -- `properties["provenance"] = "hand_written_analogue_library"` and
     `properties["score_is_placeholder"] = True` -- and the UI read neither field
     (`grep -c "provenance\\|score_is_placeholder" discovery_ui.py` was 0), rendering a bare
     "(Score: 0.85)" under a BioNeMo heading.
"""
import ast
from pathlib import Path
import sys

import pytest

ENGINE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ENGINE))

from src.molecule_generator import MoleculeGenerator  # noqa: E402

SEED = "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5"
UI = ENGINE / "app" / "discovery_ui.py"

try:
    from rdkit import Chem  # noqa: F401
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@pytest.fixture
def gen(tmp_path):
    return MoleculeGenerator(use_bionemo=False, output_dir=tmp_path)


class TestTheVcpLibraryIsNotStampedWithAnotherGene:
    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    def test_a_non_vcp_request_returns_only_the_seed(self, gen):
        mols = gen._generate_rdkit_analogues(SEED, "EGFR", 10)
        assert len(mols) == 1, "no analogue library exists for EGFR; none may be invented"
        assert mols[0].generation_method == "seed"
        assert mols[0].properties["library_available"] is False

    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    def test_no_vcp_inhibitor_is_ever_labelled_with_another_gene(self, gen):
        for gene in ("EGFR", "KRAS", "BRAF"):
            for mol in gen.generate_from_seed(SEED, gene, num_molecules=10):
                assert not (mol.name or "").startswith("VCP-Inh-"), (
                    f"{mol.name} was returned for {gene}"
                )

    @pytest.mark.skipif(not HAS_RDKIT, reason="requires RDKit")
    def test_the_vcp_path_still_works_and_stays_labelled(self, gen):
        mols = gen._generate_rdkit_analogues(SEED, "VCP", 5)
        assert len(mols) > 1
        for mol in mols[1:]:
            assert mol.properties["provenance"] == "hand_written_analogue_library"
            assert mol.properties["score_is_placeholder"] is True


class TestScoreLabelling:
    def test_a_placeholder_score_says_so(self):
        import importlib.util

        # discovery_ui imports streamlit at module scope; test the pure helper via its source.
        src = UI.read_text()
        tree = ast.parse(src)
        fn = next(n for n in tree.body
                  if isinstance(n, ast.FunctionDef) and n.name == "_score_label")
        ns = {}
        exec(compile(ast.Module(body=[fn], type_ignores=[]), str(UI), "exec"), ns)
        label = ns["_score_label"]

        class M:
            score = 0.85
            properties = {"score_is_placeholder": True,
                          "provenance": "hand_written_analogue_library"}

        assert "placeholder rank" in label(M())
        assert "hand-written" in label(M())

        class R:
            score = 0.85
            properties = {"provenance": "cloud_molmim"}

        assert label(R()) == "(Score: 0.85)"
        assert importlib.util.find_spec is not None  # keep the import meaningful


class TestTheUiNoLongerClaimsBioNeMo:
    @staticmethod
    def _code():
        """Executable code only; the removed wording is quoted in the new comments."""
        src = UI.read_text()
        lines = src.splitlines()
        drop = set()
        for node in ast.walk(ast.parse(src)):
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) \
                    and isinstance(node.value.value, str):
                drop.update(range(node.lineno, (node.end_lineno or node.lineno) + 1))
        return "\n".join(ln for i, ln in enumerate(lines, start=1)
                         if i not in drop and not ln.lstrip().startswith("#"))

    def test_the_page_does_not_advertise_bionemo_generation(self):
        code = self._code()
        assert "BioNeMo Molecule Generation" not in code
        assert "AI-generated drug candidates" not in code

    def test_the_generator_is_constructed_on_the_labeled_path(self):
        code = self._code()
        assert "MoleculeGenerator()" not in code, (
            "the default use_bionemo=True only warns and falls through to the library"
        )
        assert "MoleculeGenerator(use_bionemo=False)" in code

    def test_the_ui_reads_the_honesty_fields_the_library_sets(self):
        code = self._code()
        assert "score_is_placeholder" in code
        assert "provenance" in code
