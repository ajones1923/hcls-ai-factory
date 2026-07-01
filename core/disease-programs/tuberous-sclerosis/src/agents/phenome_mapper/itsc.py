"""
ITSC surveillance-gap analyzer (PRD §3 FR-PM-7; master paper §9).

A deterministic rule engine encoding a simplified 2021 International TSC Consensus
surveillance schedule. Given a patient's phenotype labels and age, it lists the
recommended surveillance items and flags those the record marks overdue. The full
schedule (age-stratified intervals, organ-specific nuance) is a W3-W4 refinement;
the contract and the gap-detection logic are real here.
"""
from __future__ import annotations

# type, recommended interval (months), and the phenotype-label substrings that trigger it.
# An empty trigger list means "recommended for all TSC patients".
SCHEDULE = [
    {"type": "brain MRI", "interval_months": 12, "triggers": ["tuber", "astrocytoma", "nodule"]},
    {"type": "renal MRI", "interval_months": 24, "triggers": ["angiomyolipoma", "renal", "cyst"]},
    {"type": "EEG", "interval_months": 12, "triggers": ["seizure", "spasm", "epilep"]},
    {"type": "neuropsychiatric screen (TAND-L)", "interval_months": 12, "triggers": []},
    {"type": "ophthalmology", "interval_months": 12, "triggers": []},
    {"type": "dermatology", "interval_months": 12, "triggers": ["macule", "angiofibroma", "adenoma sebaceum", "shagreen"]},
    {"type": "echocardiogram", "interval_months": 36, "triggers": ["rhabdomyoma"]},
]


def analyze(hpo_labels: list[str], overdue_hints: list[str] | None = None) -> list[dict]:
    labels = " ".join(hpo_labels).lower()
    overdue_text = " ".join(overdue_hints or []).lower()
    out = []
    for item in SCHEDULE:
        triggered = (not item["triggers"]) or any(t in labels for t in item["triggers"])
        if not triggered:
            continue
        # overdue if the record hints this surveillance type is overdue
        key = item["type"].split()[0].lower()
        status = "overdue" if key in overdue_text or item["type"].lower() in overdue_text else "due_per_schedule"
        out.append({
            "type": item["type"], "interval_months": item["interval_months"], "status": status,
            "rationale": ("recommended for all TSC patients" if not item["triggers"]
                          else f"indicated by phenotype ({', '.join(item['triggers'])})"),
        })
    return out
