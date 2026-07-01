"""
Real ADMET / toxicity prediction (C2) for the Therapeutic Discovery Engine.

Replaces the Lipinski/QED *drug-likeness heuristics* with trained models: a pretrained
ADMET suite (ChemProp D-MPNN models on the standard ADMET benchmarks) predicting ~45
absorption / distribution / metabolism / excretion / toxicity endpoints from a SMILES.

The pure logic — SMILES validation (RDKit), flag rules, and the developability verdict —
is separated from the heavy model so it is unit-testable with an injected predictor.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable


class SmilesError(ValueError):
    """Raised for an invalid SMILES string."""


def validate_smiles(smiles: str) -> str:
    """Validate + canonicalize a SMILES with RDKit."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles((smiles or "").strip())
    if mol is None:
        raise SmilesError(f"invalid SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol)


# Flag rules: concept -> (substring to find in the model's output keys, threshold, human label).
# Substring matching keeps this robust to the suite's exact endpoint key names.
FLAG_RULES: list[tuple[str, float, str]] = [
    ("hERG", 0.5, "cardiotoxicity risk (hERG)"),
    ("AMES", 0.5, "mutagenicity risk (Ames)"),
    ("DILI", 0.5, "drug-induced liver injury risk"),
    ("Carcinogen", 0.5, "carcinogenicity risk"),
    ("CYP2D6_Veith", 0.5, "CYP2D6 inhibition (DDI risk)"),
    ("CYP3A4_Veith", 0.5, "CYP3A4 inhibition (DDI risk)"),
]
# Endpoints worth reporting even when not a hard flag (context-dependent).
REPORT_KEYS = ["BBB", "Solubility", "Caco2", "Clearance", "Half_Life", "Bioavailability"]


def _find(preds: dict[str, Any], substr: str) -> tuple[str, float] | None:
    for k, v in preds.items():
        if substr.lower() in k.lower():
            try:
                return k, float(v)
            except (TypeError, ValueError):
                return None
    return None


@dataclass
class ADMETResult:
    smiles: str
    predictions: dict[str, Any]      # full model output
    flags: list[str]                 # triggered safety/DDI flags
    n_flags: int
    highlights: dict[str, float]     # curated reportable endpoints
    verdict: str                     # clean / caution / high-risk
    served_by: str                   # "admet-ai" (real) | "test"

    def to_dict(self, full: bool = False) -> dict[str, Any]:
        d = {"smiles": self.smiles, "flags": self.flags, "n_flags": self.n_flags,
             "highlights": self.highlights, "verdict": self.verdict, "served_by": self.served_by}
        if full:
            d["predictions"] = self.predictions
        return d


def summarize(smiles: str, preds: dict[str, Any], served_by: str = "admet-ai") -> ADMETResult:
    """Turn raw ADMET predictions into flags + a developability verdict (deterministic)."""
    flags = []
    for substr, thresh, label in FLAG_RULES:
        hit = _find(preds, substr)
        if hit and hit[1] >= thresh:
            flags.append(label)
    highlights = {}
    for substr in REPORT_KEYS:
        hit = _find(preds, substr)
        if hit:
            highlights[hit[0]] = round(hit[1], 4)
    verdict = "clean" if not flags else ("caution" if len(flags) <= 2 else "high-risk")
    return ADMETResult(smiles=smiles, predictions=preds, flags=flags, n_flags=len(flags),
                       highlights=highlights, verdict=verdict, served_by=served_by)


class ADMETPredictor:
    """Lazy wrapper over the pretrained ADMET suite. `_model` injectable for tests."""

    def __init__(self, _model: Callable[[str], dict] | None = None) -> None:
        self._model = _model        # callable(smiles)->preds dict for tests; None => real suite
        self._loaded = None

    def _load(self):
        if self._loaded is None:
            from admet_ai import ADMETModel
            self._loaded = ADMETModel()
        return self._loaded

    def predict(self, smiles: str) -> ADMETResult:
        canon = validate_smiles(smiles)
        if self._model is not None:                       # test path
            preds, served = self._model(canon), "test"
        else:                                             # real path
            preds, served = self._load().predict(smiles=canon), "admet-ai"
        return summarize(canon, preds, served_by=served)
