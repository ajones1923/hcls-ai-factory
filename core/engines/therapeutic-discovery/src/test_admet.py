"""Tests for real ADMET prediction (C2). Model injected; RDKit validation is real."""
import pytest

from admet_predictor import (
    ADMETPredictor, SmilesError, summarize, validate_smiles,
)

# canned outputs shaped like the pretrained suite's keys
CLEAN = {"hERG": 0.10, "AMES": 0.20, "DILI": 0.15,
         "BBB_Martins": 0.30, "Solubility_AqSolDB": -2.5, "Clearance_Hepatocyte_AZ": 20.0}
RISKY = {"hERG": 0.85, "AMES": 0.70, "DILI": 0.60, "CYP3A4_Veith": 0.90, "BBB_Martins": 0.80}


class TestSmiles:
    def test_canonicalizes_valid(self):
        assert validate_smiles("OCC") == "CCO"          # ethanol, canonicalized

    def test_rejects_invalid(self):
        with pytest.raises(SmilesError):
            validate_smiles("C(C(C")                     # unbalanced — RDKit returns None


class TestSummarize:
    def test_clean_has_no_flags(self):
        r = summarize("CCO", CLEAN)
        assert r.n_flags == 0 and r.verdict == "clean"
        # reportable endpoints surfaced
        assert any("BBB" in k for k in r.highlights) and any("Solubility" in k for k in r.highlights)

    def test_risky_flags_and_verdict(self):
        r = summarize("c1ccccc1", RISKY)
        labels = " ".join(r.flags).lower()
        assert "herg" in labels and "ames" in labels and "liver" in labels and "cyp3a4" in labels
        assert r.n_flags == 4 and r.verdict == "high-risk"

    def test_caution_for_one_or_two_flags(self):
        r = summarize("CCO", {"hERG": 0.9, "AMES": 0.1, "DILI": 0.1})
        assert r.n_flags == 1 and r.verdict == "caution"


class TestPredictor:
    def test_predict_uses_injected_model_and_validates(self):
        p = ADMETPredictor(_model=lambda smi: CLEAN)
        res = p.predict("OCC")
        assert res.served_by == "test" and res.smiles == "CCO" and res.verdict == "clean"

    def test_predict_rejects_bad_smiles_before_model(self):
        called = []
        p = ADMETPredictor(_model=lambda smi: called.append(smi) or CLEAN)
        with pytest.raises(SmilesError):
            p.predict("not_a_molecule)))")
        assert not called          # validation happened before the model

    def test_to_dict_full_includes_predictions(self):
        p = ADMETPredictor(_model=lambda smi: RISKY)
        d = p.predict("c1ccccc1").to_dict(full=True)
        assert "predictions" in d and d["n_flags"] == 4
