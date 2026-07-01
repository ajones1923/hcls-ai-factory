"""
Cohort builder (PRD §3 FR-CG-1/6/7; scripts/regen_cohort.py).

Assembles the deterministic, version-controlled 50-patient cohort: per-patient VCF +
patient.json (the agent-facing record), plus a manifest with content hashes. Same seed
-> byte-identical cohort -> identical cohort_hash (NFR-R-3). Everything is watermarked
synthetic. The faithful BAMSurgeon/Parabricks substrate is the W2 upgrade (RunPod).
"""
from __future__ import annotations

import json
from pathlib import Path

from config.settings import settings
from src.cohort.clinical import generate_clinical
from src.cohort.genomic import write_vcf
from src.cohort.spec import build_roster
from src.utils.provenance import stable_hash


def build_cohort(out_dir: str | Path | None = None, seed: int | None = None) -> dict:
    out_dir = Path(out_dir) if out_dir else settings.COHORT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    roster = build_roster(seed)
    used_seed = seed if seed is not None else _default_seed()

    patients = []
    composition: dict[str, int] = {}
    featured: dict[str, str] = {}
    for spec in roster:
        pdir = out_dir / spec.patient_id
        pdir.mkdir(parents=True, exist_ok=True)
        write_vcf(spec, pdir / "variants.vcf", seed=used_seed)
        cohort = generate_clinical(spec, f"{spec.patient_id}/variants.vcf", used_seed)
        record = json.dumps(cohort, sort_keys=True, ensure_ascii=False)
        (pdir / "patient.json").write_text(record)
        phash = stable_hash(cohort)

        key = f"{spec.gene or 'NMI'}_{spec.zygosity}"
        composition[key] = composition.get(key, 0) + 1
        if spec.featured:
            featured[spec.featured] = spec.patient_id
        patients.append({
            "patient_id": spec.patient_id, "gene": spec.gene, "zygosity": spec.zygosity,
            "vaf": spec.vaf, "featured": spec.featured, "patient_hash": phash,
        })

    manifest = {
        "watermark": "SYNTHETIC — TSC Intelligence Engine demonstration cohort",
        "seed": used_seed,
        "n_patients": len(patients),
        "composition": composition,
        "featured": featured,
        "patients": patients,
        "cohort_hash": stable_hash([p["patient_hash"] for p in patients]),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def _default_seed() -> int:
    import yaml

    with open(settings.DEMO_CONFIG_PATH) as f:
        return yaml.safe_load(f)["cohort"]["seed"]
