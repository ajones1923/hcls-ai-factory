"""HPO ontology validator (src/utils/hpo) and its wiring into the Phenome Mapper."""
import pytest

from src.utils.hpo import get_ontology

onto = get_ontology()
needs_onto = pytest.mark.skipif(not onto.available, reason="HPO release not present (run scripts/fetch_ontology.sh)")


@needs_onto
def test_valid_id_is_authoritative_and_relabels():
    # a valid id with a wrong label -> id trusted, label corrected to the official term
    r = onto.resolve("HP:0009718", "SEGA brain tumor")
    assert r["hpo_id"] == "HP:0009718"
    assert r["label"] == "Subependymal giant-cell astrocytoma"
    assert r["status"] == "relabeled"


@needs_onto
def test_bad_id_recovered_from_label():
    r = onto.resolve("HP:9999999", "Cortical tubers")
    assert r["hpo_id"] == "HP:0009717"
    assert r["status"] == "recovered"


@needs_onto
def test_hallucinated_term_is_unknown():
    r = onto.resolve("HP:9999998", "a phenotype that does not exist anywhere")
    assert r["status"] == "unknown"


@needs_onto
def test_cohort_phenome_is_fully_grounded():
    from app._engine import get_engine, featured
    orch, _ = get_engine()
    for pid in featured().values():
        prof = orch.store.projection(pid).get("hpo_profile") or {}
        assert prof.get("ontology_release", "unavailable") != "unavailable"
        terms = prof["hpo_terms"]
        assert terms, f"{pid} has no HPO terms"
        # every structured term resolves cleanly (no 'unknown'/'unverified') and ids are canonical
        for t in terms:
            assert t["validation"] in {"verified", "relabeled", "remapped", "recovered"}
            assert onto.label_for(t["hpo_id"]) == t["label"]


def test_passthrough_when_unavailable(monkeypatch):
    # the engine must still run if the release is absent
    if onto.available:
        from src.utils import hpo
        stub = hpo.HPOOntology.__new__(hpo.HPOOntology)
        stub.available = False
        stub._id_to_label, stub._name_to_id, stub._alt = {}, {}, {}
        r = stub.resolve("HP:0009718", "anything")
        assert r["status"] == "unverified" and r["hpo_id"] == "HP:0009718"
