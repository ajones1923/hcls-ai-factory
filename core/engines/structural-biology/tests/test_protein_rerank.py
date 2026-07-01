"""B3: exact Smith-Waterman re-rank (Biopython local alignment)."""
from protein_rerank import sw_align, rerank
from protein_search import ProteinSearchIndex, InMemoryBackend
from esm2_embeddings import ESM2Embedder

AAS = "ACDEFGHIKLMNPQRSTVWY"
def fake_embed(seq):
    v = [seq.count(c) for c in AAS]
    n = (sum(x * x for x in v)) ** 0.5 or 1.0
    return [x / n for x in v]


class TestSW:
    def test_identical_full_identity(self):
        r = sw_align("MKTAYIAKQR", "MKTAYIAKQR")
        assert r["pct_identity"] == 1.0 and r["matches"] == 10

    def test_divergent_lower_score(self):
        # local SW aligns only the conserved core ("MKT") so %identity can be 1.0;
        # the SW *score* is the real discriminator (fewer aligned residues -> lower score)
        full = sw_align("MKTAYIAKQR", "MKTAYIAKQR")["sw_score"]
        partial = sw_align("MKTAYIAKQR", "MKTGGGGGGG")["sw_score"]
        assert partial < full

    def test_dissimilar_zero_identity(self):
        r = sw_align("KKKKKKKKKK", "EEEEEEEEEE")   # K and E are never identical
        assert r["pct_identity"] == 0.0 and r["matches"] == 0

    def test_rerank_sorts_by_sw(self):
        hits = [{"name": "a", "cosine": 0.9}, {"name": "b", "cosine": 0.95}]
        lut = {"a": "MKTAYIAKQR", "b": "GGGGGGGGGG"}
        out = rerank("MKTAYIAKQR", hits, lut)
        assert out[0]["name"] == "a"            # exact match wins on SW despite lower cosine


class TestRerankSearch:
    def test_search_rerank_adds_fields_and_orders(self):
        emb = ESM2Embedder(_embed_fn=fake_embed, dim=20)
        idx = ProteinSearchIndex(embedder=emb, backend=InMemoryBackend(dim=20))
        idx.add([{"id": 1, "name": "exact", "sequence": "MKTAYIAKQRQI"},
                 {"id": 2, "name": "near", "sequence": "MKTAYIAKQRQL"},
                 {"id": 3, "name": "far", "sequence": "GGGGWWWWPPPP"}])
        res = idx.search("MKTAYIAKQRQI", top_k=3, rerank=True)
        assert res[0]["name"] == "exact" and res[0]["pct_identity"] == 1.0
        assert all("sw_score" in r for r in res)
