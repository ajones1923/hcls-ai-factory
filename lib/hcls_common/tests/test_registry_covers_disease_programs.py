"""Every component directory must map to a registered capability — programs included.

Found by the claim-integrity review of `tuberous-sclerosis`, 2026-11-04, while check 1
reported "no `live` capability declares an endpoint" for a program whose capability was at
that moment answering on :8560.

`validate_registry.py` loops over `dirs("engines") + dirs("agents")`. `core/disease-programs`
was not in that list, so the flagship program was never checked, had no `COVERAGE` entry, and
was invisible to everything built on the map: `subject_facts.py` matched it zero capabilities,
and the 85 generated documentation assets stated that the disease program registers nothing.

The gate could not see the directory it was meant to be gating — the same shape as
`validate_registry` checking supervised-port-to-registry but not the reverse (nine `live`
capabilities serving nothing), and as the cardiology seeder's success-on-failure returns.
Pinning the directory list is what stops the second disease program repeating it.
"""
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[3]
VALIDATOR = REPO / "scripts" / "validate_registry.py"
REGISTRY = REPO / "lib" / "hcls_common" / "capabilities.json"

pytestmark = pytest.mark.skipif(
    not (VALIDATOR.exists() and REGISTRY.exists()),
    reason="run from a source checkout; scripts/ or capabilities.json not found")


def coverage() -> dict[str, list[str]]:
    """Parse COVERAGE out of the validator, the same way subject_facts.py does."""
    m = re.search(r"COVERAGE\s*[:=][^{]*\{(.*?)\n\}", VALIDATOR.read_text(), re.S)
    assert m, "COVERAGE table not found in validate_registry.py"
    out = {}
    for line in m.group(1).splitlines():
        km = re.match(r'\s*"([^"]+)":\s*\[(.*?)\]', line)
        if km:
            out[km.group(1)] = re.findall(r'"([^"]+)"', km.group(2))
    return out


def component_dirs(kind: str) -> list[str]:
    base = REPO / "core" / kind
    if not base.is_dir():
        return []
    return sorted(f"core/{kind}/{p.name}" for p in base.iterdir()
                  if p.is_dir() and not p.name.startswith("."))


class TestEveryComponentIsCovered:
    @pytest.mark.parametrize("kind", ["engines", "agents", "disease-programs"])
    def test_every_directory_has_a_coverage_entry(self, kind):
        missing = [c for c in component_dirs(kind) if c not in coverage()]
        assert not missing, (
            f"{missing} have no COVERAGE entry, so they are invisible to subject_facts.py "
            "and to every document generated from it")

    def test_the_validator_actually_checks_disease_programs(self):
        """The omission was in the loop, not in the table. Pin the loop."""
        src = VALIDATOR.read_text()
        loop = re.search(r"for comp in (dirs\(.*?\)(?:\s*\+\s*dirs\(.*?\))*)\s*:", src)
        assert loop, "could not find the coverage loop in validate_registry.py"
        assert "disease-programs" in loop.group(1), (
            "validate_registry.py iterates engines and agents only; a disease program with "
            "no COVERAGE entry would pass unnoticed, which is how the flagship went "
            "unmapped. Loop found: " + loop.group(1))


class TestTheFlagshipSpecifically:
    def test_tuberous_sclerosis_maps_to_its_registered_capability(self):
        assert coverage().get("core/disease-programs/tuberous-sclerosis") == [
            "tuberous-sclerosis-engine"]

    def test_subject_facts_resolves_that_capability(self):
        """The end the defect was visible from: zero capabilities on the flagship."""
        import json
        reg = {c["id"] for c in json.loads(REGISTRY.read_text())["capabilities"]}
        mapped = coverage().get("core/disease-programs/tuberous-sclerosis", [])
        assert mapped and set(mapped) <= reg, (
            f"{mapped} is not registered in capabilities.json")
