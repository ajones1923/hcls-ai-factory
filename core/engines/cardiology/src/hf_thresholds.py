"""Single source of truth for the heart-failure boundaries and recommendation
classes this engine applies.

Four surfaces of this engine used to carry their own copies of these values and
disagreed about three of them, so the same patient got different advice from
`/gdmt/optimize`, `/workflow/heart-failure` and the route's built-in fallback:

* the HFrEF LVEF boundary (`< 40` in `api/routes/cardio_clinical.py` against
  `<= 40` in `GDMTOptimizer`, `HeartFailureWorkflow` and `risk_calculators`), so
  an LVEF of exactly 40.0% dropped from four Class I "INITIATE" recommendations
  to "Consider";
* the QRS duration at which CRT becomes Class 2a (120 ms in `GDMTOptimizer`,
  130 ms in `HeartFailureWorkflow`), so an LVEF 32% / QRS 125 ms / LBBB patient
  got a Class 2a CRT recommendation from one endpoint and nothing from the other;
* the recommendation class for SGLT2i in HFmrEF and HFpEF (Class I in
  `GDMTOptimizer`, Class 2a in `HeartFailureWorkflow` and the route) — the
  difference between "should" and "is reasonable" at the point of care.

Every module now imports these names instead of restating the numbers.

Sources
-------
LVEF boundaries and the SGLT2i class: 2022 AHA/ACC/HFSA Guideline for the
Management of Heart Failure (Heidenreich PA, Bozkurt B, Aguilar D, et al.,
*Circulation* 2022;145:e895-e1032; *J Am Coll Cardiol* 2022;79(17):e263-e421,
doi:10.1016/j.jacc.2021.12.012).  Table 4 classifies HFrEF as LVEF <=40%,
HFmrEF as LVEF 41-49% and HFpEF as LVEF >=50%; the ACC's summary of the
guideline quotes both verbatim: "HFrEF (HF with reduced ejection fraction [EF]):
LVEF <=40%" and "SGLT2i can also be beneficial in patients with HFmrEF and HFpEF
(Class of Recommendation 2a)"
(https://www.acc.org/latest-in-cardiology/ten-points-to-remember/2022/03/29/19/53/2022-aha-acc-hfsa-heart-failure-guideline-gl-hf,
retrieved 2026-11-03).

CRT QRS thresholds: the same guideline, device-therapy recommendations
(section 7.4.1).  The publisher's full text is not retrievable from here --
ahajournals.org and jacc.org both return HTTP 403 -- so the 120 ms lower edge of
the Class 2a LBBB band rests on reproductions of the guideline's own
recommendation text rather than on the article itself, and remains an open
question in CLINICAL_REVIEW.md section 12.  The 130 ms it replaces is the 2021
ESC figure (McDonagh TA et al., *Eur Heart J* 2021;42:3599-3726), not a value of
the guideline the code cites.
"""

from __future__ import annotations

from src.models import GuidelineClass

# ── LVEF boundaries (2022 AHA/ACC/HFSA Table 4) ─────────────────────────────
# HFrEF: LVEF <= 40.  HFmrEF: 41-49, i.e. > 40 and < 50.  HFpEF: >= 50.
HFREF_LVEF_CUTOFF: float = 40.0
HFPEF_LVEF_CUTOFF: float = 50.0
# The highest integer EF inside the HFmrEF band.  Kept for reporting only --
# classification compares against HFPEF_LVEF_CUTOFF so that a non-integer EF
# such as 49.5% stays in HFmrEF instead of falling into a band whose banner
# reads "LVEF >=50%".
HFMREF_LVEF_UPPER_INTEGER: float = 49.0

# ── CRT QRS duration thresholds, ms (2022 AHA/ACC/HFSA section 7.4.1) ───────
CRT_QRS_CLASS_I_MS: int = 150
CRT_QRS_CLASS_IIA_MS: int = 120

# ── SGLT2i recommendation class outside HFrEF ──────────────────────────────
# Class 2a in both HFmrEF and HFpEF.  In HFrEF it is Class 1.
SGLT2I_CLASS_HFREF: GuidelineClass = GuidelineClass.CLASS_I
SGLT2I_CLASS_HFMREF: GuidelineClass = GuidelineClass.CLASS_IIA
SGLT2I_CLASS_HFPEF: GuidelineClass = GuidelineClass.CLASS_IIA
# Human-readable forms for recommendation prose.
SGLT2I_CLASS_HFMREF_TEXT: str = "Class 2a"
SGLT2I_CLASS_HFPEF_TEXT: str = "Class 2a"

__all__ = [
    "HFREF_LVEF_CUTOFF",
    "HFPEF_LVEF_CUTOFF",
    "HFMREF_LVEF_UPPER_INTEGER",
    "CRT_QRS_CLASS_I_MS",
    "CRT_QRS_CLASS_IIA_MS",
    "SGLT2I_CLASS_HFREF",
    "SGLT2I_CLASS_HFMREF",
    "SGLT2I_CLASS_HFPEF",
    "SGLT2I_CLASS_HFMREF_TEXT",
    "SGLT2I_CLASS_HFPEF_TEXT",
]
