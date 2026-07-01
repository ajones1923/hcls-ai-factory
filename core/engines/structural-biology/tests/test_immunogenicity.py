"""B10: MHCflurry immunogenicity scan (windowing/burden logic; predictor monkeypatched)."""
import immunogenicity


def test_scan_burden(monkeypatch):
    def fake_predict(peptides, alleles):
        return [{"peptide": p, "best_allele": alleles[0],
                 "presentation_score": 0.9 if p == "NLVPMVATV" else 0.1} for p in peptides]
    monkeypatch.setattr(immunogenicity, "predict", fake_predict)
    scan = immunogenicity.scan_sequence("NLVPMVATVAAAAAAAA", ["HLA-A*02:01"], length=9)
    assert scan["n_strong_binders"] == 1
    assert 0.0 < scan["immunogenic_burden"] < 1.0
    assert scan["top"][0]["peptide"] == "NLVPMVATV"


def test_empty_for_short_sequence(monkeypatch):
    monkeypatch.setattr(immunogenicity, "predict", lambda p, a: [])
    scan = immunogenicity.scan_sequence("MKT", ["HLA-A*02:01"], length=9)
    assert scan["n_peptides"] == 0 and scan["immunogenic_burden"] == 0.0
