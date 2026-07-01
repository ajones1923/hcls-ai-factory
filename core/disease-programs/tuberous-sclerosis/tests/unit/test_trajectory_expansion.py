"""Multi-quantity trajectory: eGFR + seizures, graded crossing, survival probability,
Bayesian slope, surveillance cadence (FR-TM-1/2/4/5)."""
from app._engine import featured, get_engine


def _quants(tag):
    orch, _ = get_engine()
    tr = orch.store.projection(featured()[tag]).get("trajectory") or {}
    return {q["quantity"]: q for q in tr.get("quantities", [])}


def test_patient_c_models_egfr_and_seizures():
    q = _quants("C")
    assert {"AML", "eGFR", "seizure_frequency"} <= set(q)
    egfr = q["eGFR"]
    assert egfr["direction"] == "decreasing" and egfr["unit"] == "mL/min/1.73m2"
    # declining renal function heads toward the 60 floor within ~18 months
    assert egfr["crossing_grade"] in {"likely", "possible"}
    assert egfr["months_to_threshold"] is not None


def test_graded_crossing_and_survival_probability():
    q = _quants("C")
    seiz = q["seizure_frequency"]
    # rising refractory seizures -> high probability of crossing escalation threshold
    assert seiz["crossing_grade"] == "likely"
    assert seiz["crossing_probability"]["m18"] >= 0.5
    # a stable quantity is graded unlikely with low probability
    qb = _quants("B")
    assert qb["seizure_frequency"]["crossing_grade"] == "unlikely"
    assert qb["seizure_frequency"]["crossing_probability"]["m18"] < 0.15


def test_surveillance_cadence_tightens_for_likely_crossing():
    egfr = _quants("C")["eGFR"]
    sr = egfr["surveillance_recommendation"]
    assert sr["recommended_interval_months"] < sr["itsc_floor_months"]   # tightened
    # a stable quantity holds the floor
    sb = _quants("B")["seizure_frequency"]["surveillance_recommendation"]
    assert sb["recommended_interval_months"] == sb["itsc_floor_months"]


def test_bayesian_slope_population_prior_for_sega_weak_for_others():
    qb = _quants("B")
    assert qb["SEGA"]["bayesian"]["population_informed"] is True
    assert _quants("C")["eGFR"]["bayesian"]["population_informed"] is False


def test_sega_lesion_view_preserved_for_backcompat():
    orch, _ = get_engine()
    tr = orch.store.projection(featured()["B"]).get("trajectory") or {}
    lesions = tr["lesions"]
    assert lesions and lesions[0]["lesion"] == "SEGA"   # legacy consumers still see SEGA/AML
    assert lesions[0]["crosses_in_12_18mo_window"] is True
