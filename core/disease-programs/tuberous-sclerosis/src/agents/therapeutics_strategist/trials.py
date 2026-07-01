"""
Synthetic ClinicalTrials.gov snapshot + eligibility-criterion matcher (PRD §3 FR-TR-2/3).

Real per-criterion evaluation against the patient profile: each criterion resolves to
eligible / ineligible / requires_clarification, defaulting to requires_clarification when
the record is insufficient (honest). The live ClinicalTrials.gov parser is the W5-W6
upgrade; the matching contract is real here.
"""
from __future__ import annotations

SNAPSHOT = [
    {
        "nct": "NCT-SYN-001",
        "title": "Next-generation selective mTORC1 inhibitor for treatment-resistant TSC epilepsy",
        "inclusion": [
            {"text": "Confirmed TSC1/TSC2 molecular diagnosis", "key": "molecular_dx"},
            {"text": "Age 6-30 years", "key": "age", "min": 6, "max": 30},
            {"text": "Treatment-resistant focal seizures", "key": "refractory_epilepsy"},
            {"text": "Stable anti-seizure regimen for >=30 days", "key": "stable_regimen"},
        ],
        "exclusion": [{"text": "Prior next-generation mTORC1 inhibitor exposure", "key": "prior_mtorc1"}],
    },
    {
        "nct": "NCT-SYN-002",
        "title": "mTOR inhibitor for TSC-associated renal angiomyolipoma >=3 cm",
        "inclusion": [
            {"text": "Confirmed TSC diagnosis", "key": "molecular_dx"},
            {"text": "Renal angiomyolipoma >= 3 cm", "key": "aml_ge_3cm"},
        ],
        "exclusion": [],
    },
]


def _eval(crit: dict, profile: dict) -> str:
    k = crit["key"]
    if k == "molecular_dx":
        return "eligible" if profile.get("molecular_dx_confirmed") else "requires_clarification"
    if k == "age":
        a = profile.get("age")
        if a is None:
            return "requires_clarification"
        return "eligible" if crit["min"] <= a <= crit["max"] else "ineligible"
    if k == "refractory_epilepsy":
        return {True: "eligible", False: "ineligible", None: "requires_clarification"}[
            profile.get("refractory_epilepsy")
        ]
    if k == "aml_ge_3cm":
        aml = profile.get("aml_max_cm")
        if aml is None:
            return "requires_clarification"
        return "eligible" if aml >= 3.0 else "ineligible"
    if k in ("stable_regimen", "prior_mtorc1"):
        return "requires_clarification"   # not in the record -> clinician verifies
    return "requires_clarification"


def match_trials(profile: dict) -> list[dict]:
    out = []
    for trial in SNAPSHOT:
        crits = [{"text": c["text"], "status": _eval(c, profile)} for c in trial["inclusion"]]
        excl = [{"text": c["text"], "status": _eval(c, profile)} for c in trial["exclusion"]]
        statuses = [c["status"] for c in crits]
        if "ineligible" in statuses:
            verdict = "ineligible"
        elif "requires_clarification" in statuses:
            verdict = "requires_clarification"
        else:
            verdict = "eligible"
        out.append({"nct": trial["nct"], "title": trial["title"], "match": verdict,
                    "inclusion": crits, "exclusion": excl})
    return out
