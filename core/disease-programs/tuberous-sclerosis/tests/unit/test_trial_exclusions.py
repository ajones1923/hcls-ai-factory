"""An exclusion criterion must be able to make a patient ineligible.

`match_trials` evaluated each trial's exclusion criteria and returned them in the response,
but only the inclusion statuses fed the verdict, so no exclusion could ever change the match.
A patient with a documented bar read "eligible".

The eligibility words also invert for an exclusion: `_eval` answers "is this criterion
satisfied", and a satisfied exclusion means the condition is PRESENT and the patient is
excluded. Exclusion statuses are therefore reported as present / absent / unverified.
"""
from src.agents.therapeutics_strategist.trials import SNAPSHOT, match_trials


def _by_nct(matches):
    return {m["nct"]: m for m in matches}


def test_unverified_exclusion_blocks_a_clean_eligible():
    """You cannot clear someone against a bar nobody checked."""
    trial = next(t for t in SNAPSHOT if t["exclusion"])
    profile = {
        "molecular_dx_confirmed": True,
        "age": 18,
        "refractory_epilepsy": True,
        "aml_max_cm": 4.0,
    }
    m = _by_nct(match_trials(profile))[trial["nct"]]
    assert m["match"] != "eligible"
    assert any(c["status"] == "unverified" for c in m["exclusion"])


def test_exclusion_statuses_do_not_use_eligibility_words():
    """"eligible" on an exclusion row is unreadable; it must say present/absent/unverified."""
    profile = {"molecular_dx_confirmed": True, "age": 18,
               "refractory_epilepsy": True, "aml_max_cm": 4.0}
    for m in match_trials(profile):
        for crit in m["exclusion"]:
            assert crit["status"] in ("present", "absent", "unverified"), crit


def test_a_present_exclusion_makes_the_patient_ineligible(monkeypatch):
    """With the bar documented as present, the verdict must be ineligible."""
    import src.agents.therapeutics_strategist.trials as trials_mod

    real_eval = trials_mod._eval

    def fake_eval(crit, profile):
        # Report the one exclusion key as satisfied, i.e. the condition is present.
        if crit["key"] == "prior_mtorc1":
            return "eligible"
        return real_eval(crit, profile)

    monkeypatch.setattr(trials_mod, "_eval", fake_eval)
    trial = next(t for t in SNAPSHOT if t["exclusion"])
    profile = {"molecular_dx_confirmed": True, "age": 18,
               "refractory_epilepsy": True, "aml_max_cm": 4.0}
    m = _by_nct(trials_mod.match_trials(profile))[trial["nct"]]
    assert m["match"] == "ineligible"


def test_a_trial_without_exclusions_is_unaffected():
    profile = {"molecular_dx_confirmed": True, "age": 18,
               "refractory_epilepsy": True, "aml_max_cm": 4.0}
    m = _by_nct(match_trials(profile))["NCT-SYN-002"]
    assert m["match"] == "eligible"
    assert m["exclusion"] == []
