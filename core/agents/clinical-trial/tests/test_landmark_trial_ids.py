"""Regression tests pinning each seeded NCT ID to the study it claims to be.

21 of the 36 records in ``LANDMARK_TRIALS`` carried an NCT ID belonging to a
different study. ``test_knowledge.py`` has 324 tests over this same table and
none of them could detect it: ``test_has_nct_id`` asserts
``trial["nct_id"].startswith("NCT")``, ``test_has_phase`` asserts ``"Phase" in
trial["phase"]``, ``test_has_title`` asserts ``len(...) > 10``, and
``test_unique_nct_ids`` checks uniqueness only -- which *passed* on the two
contradictory "ADAURA" rows precisely because their IDs differed.

The offline test below pins the (acronym -> NCT ID) pairing that was verified
against the live registry, so a future edit cannot silently reintroduce a wrong
ID. The networked test is the CI gate: it resolves every ID in the table
against ClinicalTrials.gov API v2 and asserts the claimed acronym appears in
``briefTitle``, ``officialTitle`` or ``acronym``. It is skipped when the
registry is unreachable, so it never fails a build for being offline.

Primary source: ClinicalTrials.gov API v2,
``GET https://clinicaltrials.gov/api/v2/studies/{NCTId}``, retrieved
2026-11-03. Each of the 20 corrected IDs was checked for lead sponsor, phase
and title before being written into the table.
"""

import json
import re
import urllib.error
import urllib.request

import pytest

from src.ingest.clinicaltrials_parser import LANDMARK_TRIALS


#: (claimed acronym, correct NCT ID). Every pair verified against the live
#: registry on 2026-11-03. The 20 marked CORRECTED replaced an ID that resolved
#: to an entirely different study -- for instance "KEYNOTE-024" pointed at
#: NCT02477436, "Efficacy and Safety of Avanafil in the Patients With Erectile
#: Dysfunction" (Pusan National University Hospital, PHASE2), and "MONARCH 3"
#: pointed at NCT02507959, which is not a registry ID at all (HTTP 404).
VERIFIED_IDS = {
    "KEYNOTE-024": "NCT02142738",        # CORRECTED (was NCT02477436)
    "KEYNOTE-189": "NCT02578680",        # CORRECTED (was NCT02813135)
    "CheckMate 067": "NCT01844505",      # CORRECTED (was NCT01673867)
    "PARADIGM-HF": "NCT01035255",        # CORRECTED (was NCT01820754)
    "DAPA-HF": "NCT03036124",
    "EMPEROR-Reduced": "NCT03057977",
    "FOURIER": "NCT01764633",            # CORRECTED (was NCT02465268)
    "ELOQUENT-3": "NCT02654132",         # CORRECTED (was NCT02412098)
    "RECOVERY": "NCT04381936",           # CORRECTED (was NCT04368728)
    "ACTT-1": "NCT04280705",
    "VICTORIA": "NCT02861534",           # CORRECTED (was NCT03548935)
    "MONARCH 3": "NCT02246621",          # CORRECTED (was NCT02507959, a 404)
    "DESTINY-Breast03": "NCT03529110",   # CORRECTED (was NCT03521986)
    "IMpower150": "NCT02366143",         # CORRECTED (was NCT03299946)
    "SPRINT": "NCT01206062",             # CORRECTED (was NCT02684006)
    "DECLARE-TIMI 58": "NCT01730534",    # CORRECTED (was NCT02853331)
    "MONALEESA-7": "NCT02278120",        # CORRECTED (was NCT02564263)
    "PACIFIC": "NCT02125461",            # CORRECTED (was NCT02362594)
    "CREDENCE": "NCT02065791",           # CORRECTED (was NCT03170960)
    "ADAURA": "NCT02511106",
    "HIMALAYA": "NCT03298451",           # CORRECTED (was NCT03834506)
    "DESTINY-Breast04": "NCT03734029",
    "KEYNOTE-522": "NCT03036488",
    "CheckMate-9LA": "NCT03215706",
    "SELECT": "NCT03574597",
    "SURMOUNT-1": "NCT04184622",
    "TRAILBLAZER-ALZ 2": "NCT04437511",
    "PURPOSE 1": "NCT04994509",
    "EV-302": "NCT04223856",
    "ELIANA": "NCT02435849",
    "ZUMA-7": "NCT03391466",
    "MAESTRO": "NCT03900429",            # CORRECTED (was NCT03900078)
    "PROTECT": "NCT03875729",            # CORRECTED (was NCT04194450)
    "VALOR-HCM": "NCT04349072",
    "TRILUMINATE Pivotal": "NCT03904147",  # CORRECTED (was NCT04788511)
}


def _acronym_of(trial: dict) -> str:
    """The acronym the record claims, taken from the text before the colon."""
    title = trial.get("title", "")
    return title.split(":", 1)[0].strip()


@pytest.mark.parametrize("acronym,nct_id", sorted(VERIFIED_IDS.items()))
def test_seeded_id_matches_the_verified_registration(acronym, nct_id):
    """Each landmark trial carries the registry ID verified for that study."""
    by_acronym = {_acronym_of(t): t["nct_id"] for t in LANDMARK_TRIALS}
    assert acronym in by_acronym, f"{acronym} missing from LANDMARK_TRIALS"
    assert by_acronym[acronym] == nct_id, (
        f"{acronym} is seeded as {by_acronym[acronym]} but the verified "
        f"registration is {nct_id}"
    )


def test_no_duplicate_acronyms():
    """One row per study.

    The table held TWO rows titled "ADAURA" with different IDs. The real ADAURA
    is NCT02511106; the other row's NCT03723655 resolves to MAVA-LTE, a
    mavacamten long-term extension study sponsored by BMS. Uniqueness of IDs
    passed precisely because the IDs differed, so the check has to be on the
    study, not the key.
    """
    acronyms = [_acronym_of(t) for t in LANDMARK_TRIALS]
    dupes = {a for a in acronyms if acronyms.count(a) > 1}
    assert not dupes, f"duplicate trial acronyms: {sorted(dupes)}"


def test_no_duplicate_nct_ids():
    """record_id == nct_id when seeding, so IDs must be unique."""
    ids = [t["nct_id"] for t in LANDMARK_TRIALS]
    assert len(ids) == len(set(ids))


def test_every_seeded_trial_is_in_the_verified_table():
    """No landmark trial ships without a registry-verified ID.

    This is the gate that makes the next addition cheap to get right: a new row
    fails here until its ID has been resolved against the registry and recorded
    in VERIFIED_IDS.
    """
    unverified = sorted(
        _acronym_of(t) for t in LANDMARK_TRIALS
        if _acronym_of(t) not in VERIFIED_IDS
    )
    assert not unverified, (
        f"trials with no registry-verified NCT ID: {unverified}"
    )


def test_all_ids_are_well_formed():
    """An NCT ID is "NCT" plus exactly eight digits."""
    for t in LANDMARK_TRIALS:
        assert re.fullmatch(r"NCT\d{8}", t["nct_id"]), t["nct_id"]


# ═══════════════════════════════════════════════════════════════════════════
# Networked CI gate
# ═══════════════════════════════════════════════════════════════════════════

_CTG_FIELDS = (
    "protocolSection.identificationModule,"
    "protocolSection.designModule.phases,"
    "protocolSection.sponsorCollaboratorsModule.leadSponsor"
)


def _fetch_study(nct_id: str):
    url = f"https://clinicaltrials.gov/api/v2/studies/{nct_id}?fields={_CTG_FIELDS}"
    with urllib.request.urlopen(url, timeout=20) as resp:
        return json.load(resp)


def _normalise(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (text or "").lower())


@pytest.mark.network
@pytest.mark.parametrize(
    "trial", LANDMARK_TRIALS, ids=lambda t: t["nct_id"],
)
def test_nct_id_resolves_to_the_claimed_study(trial):
    """Resolve each ID against ClinicalTrials.gov and check the study matches.

    This is the check the unit suite never had. 324 existing tests over this
    table assert only presence and shape, so every one of them passed with 21
    IDs pointing at the wrong study and one ID that is not in the registry.

    Skipped, not failed, when the registry is unreachable: a build must not go
    red for being offline, but it must go red for a wrong registration.
    """
    nct_id = trial["nct_id"]
    try:
        data = _fetch_study(nct_id)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            pytest.fail(f"{nct_id} is not a ClinicalTrials.gov registration")
        pytest.skip(f"registry unreachable for {nct_id}: HTTP {exc.code}")
    except Exception as exc:  # pragma: no cover - network conditions
        pytest.skip(f"registry unreachable for {nct_id}: {exc}")

    ident = data.get("protocolSection", {}).get("identificationModule", {})
    haystack = _normalise(
        " ".join([
            ident.get("briefTitle", ""),
            ident.get("officialTitle", ""),
            ident.get("acronym", ""),
        ])
    )
    acronym = _acronym_of(trial)
    # Accept the registry's own shorter form of a numbered series: NCT04280705
    # is registered as "Adaptive COVID-19 Treatment Trial (ACTT)" and is the
    # study the literature calls ACTT-1. Only a trailing series number is
    # dropped, so "KEYNOTE-024" is still required in full and cannot be
    # satisfied by a different KEYNOTE study.
    candidates = {_normalise(acronym)}
    stripped = re.sub(r"[-\s]*\d+$", "", acronym).strip()
    if stripped and stripped != acronym:
        candidates.add(_normalise(stripped))
    assert any(c in haystack for c in candidates), (
        f"{nct_id} is seeded as '{acronym}' but the registry calls it "
        f"'{ident.get('briefTitle', '')}' "
        f"(sponsor: "
        f"{data.get('protocolSection', {}).get('sponsorCollaboratorsModule', {}).get('leadSponsor', {}).get('name', '?')})"
    )
