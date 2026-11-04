"""The repo scripts that iterate subjects must know about every subject.

Two of them take the roster as input and index a dict by what they find in it, so a new kind of
subject -- or a new agent -- is a KeyError or a silent omission rather than a report.

That is not hypothetical. `scripts/run_all_tests.py` grew an `orchestrator` row; `BASE` in
`check_test_depth.py` had three kinds; and the test-depth ratchet crashed with
`KeyError: 'orchestrator'` on every CI run. It looked like a working gate for as long as nothing
triggered CI, because `ci.yml` fires only on push-to-main and PR-to-main -- 35 pushes to a branch
ran it zero times. Opening the PR is what surfaced it.

These are import-and-inspect checks: no Milvus, no agent dependencies, no subprocess.
"""
import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "scripts"


def load(name: str):
    path = SCRIPTS / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def roster():
    return load("run_all_tests.py")


class TestTheTestDepthRatchetKnowsEveryKind:
    def test_every_roster_kind_has_a_source_directory(self, roster):
        depth = load("check_test_depth.py")
        kinds = {kind for kind, _ in roster.ROSTER}
        # run_all_tests appends the orchestrator outside ROSTER, so take it from the harness too.
        missing = {k for k in kinds | {"orchestrator"} if k not in depth.BASE}
        assert not missing, (
            f"check_test_depth.BASE has no entry for {sorted(missing)}; the ratchet will raise "
            f"KeyError on that row and the whole CI step fails before measuring anything")

    def test_every_source_directory_resolves(self):
        depth = load("check_test_depth.py")
        for kind, rel in depth.BASE.items():
            assert (ROOT / rel).is_dir(), f"BASE[{kind!r}] = {rel!r} does not exist"

    def test_the_orchestrator_row_measures_a_nonzero_line_count(self):
        """A path that resolves but matches no file gives per100 = 0 and silently exempts it."""
        depth = load("check_test_depth.py")
        assert depth.source_loc("orchestrator", "hcls-orchestrator") > 0


class TestTheCorpusWeightCheckerKnowsEverySubject:
    """Every subject either declares a workflow weight table or is recorded as not having one."""

    @staticmethod
    def _subjects_with_a_weight_table():
        import re
        pat = re.compile(r"^(?:WORKFLOW_COLLECTION_(?:WEIGHTS|BOOST)|STAGE_COLLECTION_BOOST)\s*[:=]",
                         re.M)
        found = set()
        for base in ("core/agents", "core/engines", "core/disease-programs"):
            for d in sorted((ROOT / base).glob("*")):
                if not d.is_dir():
                    continue
                for py in list(d.glob("src/*.py")) + list(d.glob("agent/src/*.py")):
                    if pat.search(py.read_text(errors="ignore")):
                        found.add(d.name)
                        break
        return found

    def test_no_subject_with_a_weight_table_is_missing_from_the_map(self):
        checker = load("check_corpus_weights.py")
        declared = set(checker.WEIGHT_TABLES) | set(checker.NO_WEIGHT_TABLE)
        missing = self._subjects_with_a_weight_table() - declared
        assert not missing, (
            f"{sorted(missing)} declare a collection weight table that "
            f"check_corpus_weights.py does not read, so their coverage goes unreported")

    def test_the_subjects_recorded_as_having_no_table_really_have_none(self):
        checker = load("check_corpus_weights.py")
        wrong = set(checker.NO_WEIGHT_TABLE) & self._subjects_with_a_weight_table()
        assert not wrong, f"{sorted(wrong)} do declare a weight table after all"

    def test_a_subject_is_not_in_both_lists(self):
        checker = load("check_corpus_weights.py")
        assert not (set(checker.WEIGHT_TABLES) & set(checker.NO_WEIGHT_TABLE))

    @pytest.mark.parametrize("subject", [
        "clinical-trial", "cart", "pharmacogenomics", "neurology",
        "rare-disease-diagnostic", "single-cell", "cardiology",
    ])
    def test_each_mapped_table_exists_where_the_map_says(self, subject):
        """Checked by reading the file, not importing it -- these modules need their own venv."""
        checker = load("check_corpus_weights.py")
        rel, mod, attr = checker.WEIGHT_TABLES[subject]
        path = ROOT / "core" / rel / (mod.replace(".", "/") + ".py")
        assert path.exists(), f"{path} does not exist"
        assert f"\n{attr}" in path.read_text(errors="ignore"), (
            f"{attr} is not defined in {path}")
