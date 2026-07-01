"""
Clinical layer (PRD §3 FR-CG-2/4/5; master paper §15).

Deterministic structured phenotypes (HPO), longitudinal SEGA series, TAND signals,
surveillance gaps, and synthetic notes per patient — assembled into the `cohort` dict
the agents consume. Offline it uses templates; with TSC_ANTHROPIC_API_KEY set, the
note/imaging text is the hook for frontier-model generation (src/cohort/notes_gen,
imaging_reports). The structured signals (which the current agents read) are real here.
"""
from __future__ import annotations

import random

from src.cohort.notes_gen import build_notes
from src.cohort.spec import PatientSpec

# HPO terms by organ system — codes verified against the HPO release (src/utils/hpo;
# the Phenome Mapper re-validates every term at runtime). Canonical TSC major-feature codes.
HPO = {
    "cortical_tuber": {"hpo_id": "HP:0009717", "label": "Cortical tubers"},
    "sega": {"hpo_id": "HP:0009718", "label": "Subependymal giant-cell astrocytoma"},
    "seizure": {"hpo_id": "HP:0001250", "label": "Seizure"},
    "infantile_spasms": {"hpo_id": "HP:0012469", "label": "Infantile spasms"},
    "aml": {"hpo_id": "HP:0006772", "label": "Renal angiomyolipoma"},
    "hypomelanotic": {"hpo_id": "HP:0009719", "label": "Hypomelanotic macule"},
    "angiofibroma": {"hpo_id": "HP:0009720", "label": "Adenoma sebaceum"},
    "id": {"hpo_id": "HP:0001249", "label": "Intellectual disability"},
    "asd": {"hpo_id": "HP:0000717", "label": "Autism"},
}


def _phenotypes(spec: PatientSpec, rng: random.Random) -> list[dict]:
    terms = [HPO["cortical_tuber"], HPO["seizure"], HPO["hypomelanotic"]]
    if spec.age < 3 or rng.random() < 0.3:
        terms.append(HPO["infantile_spasms"])
    if spec.age >= 8 or rng.random() < 0.5:
        terms.append(HPO["aml"])
    if "sega_series_cm" in spec.extras:
        terms.append(HPO["sega"])
    if spec.severity in ("moderate-severe", "severe"):
        terms.append(rng.choice([HPO["id"], HPO["asd"]]))
    out = []
    for t in terms:
        out.append({**t, "onset": "infancy" if spec.age < 5 else "childhood",
                    "evidence_span": "synthetic clinical note"})
    return out


def _tand_clusters(spec: PatientSpec) -> dict:
    """The patient's TAND cluster load — explicit for featured patients (presence of the
    key is honored even when empty, so a quiet patient stays quiet), severity-derived
    otherwise. Drives both the note-embedded discourse markers and the structured signals."""
    if "tand_clusters" in spec.extras:
        return dict(spec.extras["tand_clusters"])
    if spec.severity == "severe":
        return {"behavioral": 2, "academic": 1}
    return {}


# NB: TAND signal now lives in the note spans (the real, offset-grounded source). The
# structured `tand_signals` channel is reserved for genuine structured-record flags / the
# LLM discourse pass, and is left empty here so note spans are not double-counted.


def generate_clinical(spec: PatientSpec, vcf_path: str, seed: int) -> dict:
    rng = random.Random(f"{seed}:{spec.patient_id}")
    phenos = _phenotypes(spec, rng)
    genomics = {
        "gene": spec.gene,
        "variant": (f"{spec.variant['cdna']} {spec.variant['protein']}" if spec.variant else None),
        "vaf": spec.vaf if spec.gene else None,
        "expected_acmg": spec.expected_acmg,
        "vcf_path": vcf_path,
    }
    if spec.zygosity == "mosaic":
        genomics["tissue_vaf"] = spec.vaf            # affected-tissue VAF (mosaic recovery)
        genomics["tissue"] = "resected tuber (epilepsy surgery)"
    cohort: dict = {
        "patient_id": spec.patient_id,
        "featured": spec.featured,
        "demographics": {"age": spec.age, "sex": spec.sex, "severity": spec.severity},
        "genomics": genomics,
        "phenotypes": phenos,
        "tand_signals": [],   # see note above — signal is span-grounded in `notes`
        "surveillance_gaps": (["abdominal MRI overdue"] if rng.random() < 0.25 else []),
        "notes": build_notes(spec, phenos, _tand_clusters(spec)),
        "rag_sources": [],
    }
    for k in ("sega_series_cm", "sega_months", "sega_location", "aml_max_cm", "aml_series_cm",
              "aml_months", "therapy", "refractory_epilepsy", "seizures",
              "egfr_series", "egfr_months", "seizure_freq_series", "seizure_months"):
        if k in spec.extras:
            cohort[k] = spec.extras[k]
    if spec.extras.get("trial_eligible"):
        cohort["trial_matches"] = [{
            "nct": "NCT-synthetic", "match": "requires_clarification",
            "criterion": "stable regimen ≥30 days — verify",
        }]
    return cohort
