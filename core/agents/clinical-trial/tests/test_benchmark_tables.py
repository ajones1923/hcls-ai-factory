"""Tests pinning the benchmark tables to the primary sources they cite.

Two hand-transcribed reference tables were materially wrong against the
documents named in their own comments:

* ``ProtocolDesignWorkflow._HISTORICAL_SUCCESS`` (F13) -- every value differed
  from BIO Figure 2, cardiology Phase II by 81% relative, and the table carried
  a ``phase_iv`` row that has no counterpart in the figure at all.
* ``DiversityAssessmentWorkflow._US_CENSUS`` (F14) -- labelled "US Census 2020
  approximate proportions" and not the 2020 Census; multiracial understated
  four-fold.

Primary sources, both retrieved 2026-11-03 and read directly:

* BIO / Informa Pharma Intelligence / QLS Advisors, "Clinical Development
  Success Rates and Contributing Factors 2011-2020", Figure 2 ("Phase
  transition success rates by disease area"), p.8. PDF from
  https://www.cphi-online.com/46/resourcefile/11/36/30/2021%20Clinical%20Development%20Success%20Rates%202011-2020%20v17.pdf
* U.S. Census Bureau, "Improved Race and Ethnicity Measures Reveal U.S.
  Population Is Much More Multiracial", 12 August 2021,
  https://www.census.gov/library/stories/2021/08/improved-race-ethnicity-measures-reveal-united-states-population-much-more-multiracial.html
"""

import pytest

from src.clinical_workflows import (
    DiversityAssessmentWorkflow,
    ProtocolDesignWorkflow,
)


# ═══════════════════════════════════════════════════════════════════════════
# F13 — BIO Figure 2 phase transition success rates
# ═══════════════════════════════════════════════════════════════════════════

#: BIO Figure 2, 'POS' column, verbatim. (category, transition, rate).
#: The pre-fix value is given where it differed.
BIO_FIGURE_2 = [
    # category          transition          BIO     (was)
    ("oncology",        "phase_i",          0.488),  # 0.63
    ("oncology",        "phase_ii",         0.246),  # 0.29
    ("oncology",        "phase_iii",        0.477),  # 0.36
    ("oncology",        "nda_to_approval",  0.920),
    ("cardiology",      "phase_i",          0.500),  # 0.72
    ("cardiology",      "phase_ii",         0.210),  # 0.38 -- +81% relative
    ("cardiology",      "phase_iii",        0.552),  # 0.52
    ("cardiology",      "nda_to_approval",  0.825),
    ("neurology",       "phase_i",          0.477),  # 0.68
    ("neurology",       "phase_ii",         0.268),  # 0.25
    ("neurology",       "phase_iii",        0.531),  # 0.30
    ("neurology",       "nda_to_approval",  0.867),
    ("immunology",      "phase_i",          0.552),  # 0.70 (BIO: Autoimmune)
    ("immunology",      "phase_ii",         0.314),  # 0.35
    ("immunology",      "phase_iii",        0.653),  # 0.48
    ("immunology",      "nda_to_approval",  0.941),
    ("infectious_disease", "phase_i",       0.578),  # 0.71
    ("infectious_disease", "phase_ii",      0.384),  # 0.40
    ("infectious_disease", "phase_iii",     0.640),  # 0.58
    ("infectious_disease", "nda_to_approval", 0.929),
    ("metabolic",       "phase_ii",         0.450),
    ("hematology",      "phase_iii",        0.768),
    ("respiratory",     "phase_ii",         0.219),
    ("default",         "phase_i",          0.520),  # 0.67 (BIO: All indications)
    ("default",         "phase_ii",         0.289),  # 0.33 -- the one value
                                                     #   that was already right
    ("default",         "phase_iii",        0.578),  # 0.45
    ("default",         "nda_to_approval",  0.906),
]


@pytest.mark.parametrize("category,transition,expected", BIO_FIGURE_2)
def test_success_rate_matches_bio_figure_2(category, transition, expected):
    """Each rate equals the published BIO Figure 2 value exactly."""
    table = ProtocolDesignWorkflow._HISTORICAL_SUCCESS
    assert category in table, category
    assert table[category][transition] == pytest.approx(expected, abs=0.0005), (
        f"{category}/{transition}: {table[category][transition]} != {expected}"
    )


def test_cardiology_phase_ii_is_not_overstated():
    """The worst single value: 0.38 against a measured 0.210.

    Clinical/programme consequence: a cardiology team runs
    ProtocolDesignWorkflow on a Phase 2 HFpEF programme. The agent reports a
    38% probability of advancing -- citing industry benchmarks -- when the
    measured decade-long rate is 21%. The programme is resourced on roughly
    double the true likelihood: a Phase 3 site network is pre-committed, the
    Phase 2 is under-powered to the point where "advance" is judged on a
    favourable trend, and the confirmatory trial enrols thousands of patients
    into a study with a far lower prior of benefit than the sponsor was told.
    """
    rate = ProtocolDesignWorkflow._HISTORICAL_SUCCESS["cardiology"]["phase_ii"]
    assert rate == pytest.approx(0.210, abs=0.0005)


def test_neurology_phase_iii_is_not_understated():
    """0.30 against a measured 0.531 -- the error running the other way.

    A neurology programme that would have advanced is killed on a reported
    30% probability of success, and its patients lose access.
    """
    rate = ProtocolDesignWorkflow._HISTORICAL_SUCCESS["neurology"]["phase_iii"]
    assert rate == pytest.approx(0.531, abs=0.0005)


def test_no_phase_iv_success_rate_is_invented():
    """BIO's fourth column is NDA/BLA -> Approval, not Phase IV.

    Both tables carried a ``phase_iv`` entry of 0.82-0.90 with no counterpart
    in Figure 2: the code had taken a submission-to-approval rate (90.6% for
    all indications) and relabelled it a Phase IV probability of success. That
    is a wrong-field read, not a rounding error.
    """
    for category, rates in ProtocolDesignWorkflow._HISTORICAL_SUCCESS.items():
        assert "phase_iv" not in rates, category
        assert "nda_to_approval" in rates, category


def test_rare_disease_has_no_invented_rate():
    """BIO Figure 2 does not break out rare disease, so no rate is cited.

    The old table had rare_disease at 0.74/0.42/0.55/0.90 -- the most
    optimistic row in the whole table, and sourced to nothing. BIO's
    rare-disease analysis is a separate section that was not retrieved.
    """
    table = ProtocolDesignWorkflow._HISTORICAL_SUCCESS
    assert "rare_disease" not in table
    assert "rare_disease" in ProtocolDesignWorkflow._HISTORICAL_SUCCESS_UNSOURCED


def test_the_table_names_its_source():
    """The citation names the report, edition, figure and page, so the next
    audit can diff it."""
    cite = ProtocolDesignWorkflow._HISTORICAL_SUCCESS_CITATION
    assert "2011-2020" in cite
    assert "Figure 2" in cite
    assert "BIO" in cite


def test_unknown_phase_yields_no_probability_of_success():
    """An unparsed phase must not land on a real rate.

    ``.get(phase, 0.40)`` routed "Phase 1", "PHASE1", "Early Phase 1" and
    ``["PHASE1","PHASE2"]`` to 0.40 -- near the old Phase II value -- so a
    first-in-human study was scored with a rate that was not its own.
    """
    wf = ProtocolDesignWorkflow()
    # execute() is called directly: ProtocolDesignWorkflow.preprocess has its
    # own defect (F12) -- its valid_phases set lacks early_phase_1 and
    # not_applicable, so it silently WARNS a CTG "Phase 1" into phase_ii
    # before execute ever sees it. That is out of scope for this branch; this
    # test pins the scoring path's own behaviour on an unmapped phase.
    result = wf.execute({
        "indication": "oncology", "phase": "early_phase_1",
        "comparator": "placebo", "mechanism_of_action": "novel agent",
        "design_type": "double_blind", "num_arms": 2,
    })
    joined = " ".join(result.findings)
    assert "No phase transition success rate is available" in joined
    assert result.confidence is None


def test_known_phase_still_reports_its_rate():
    """The normal path still works and cites the source in the findings."""
    wf = ProtocolDesignWorkflow()
    inputs = wf.preprocess({
        "indication": "oncology", "phase": "phase_ii",
        "comparator": "placebo", "mechanism_of_action": "novel agent",
    })
    result = wf.execute(inputs)
    joined = " ".join(result.findings)
    assert "24.6%" in joined
    assert "Figure 2" in joined


def test_design_adjustment_is_labelled_a_heuristic():
    """The adjusted number is not published and must not read as if it were."""
    wf = ProtocolDesignWorkflow()
    inputs = wf.preprocess({
        "indication": "oncology", "phase": "phase_ii",
        "comparator": "placebo", "mechanism_of_action": "novel agent",
    })
    result = wf.execute(inputs)
    joined = " ".join(result.findings)
    assert "not published rates" in joined


def test_no_ungrounded_ich_or_fda_citations():
    """The reference list cited ICH E6/E8/E9 and an FDA guidance that does not
    exist, on every result, deriving nothing from any of them.

    The old list ended with an f-string, ``f"FDA Guidance: {category.title()}
    Trial Design (current)"``, which fabricates a document title per
    therapeutic area.
    """
    wf = ProtocolDesignWorkflow()
    inputs = wf.preprocess({
        "indication": "oncology", "phase": "phase_ii",
        "comparator": "placebo", "mechanism_of_action": "novel agent",
    })
    result = wf.execute(inputs)
    joined = " ".join(result.guideline_references)
    assert "FDA Guidance: Oncology Trial Design" not in joined
    assert "ICH E9" not in joined
    assert "BIO" in joined


# ═══════════════════════════════════════════════════════════════════════════
# F14 — 2020 Census proportions
# ═══════════════════════════════════════════════════════════════════════════

#: 2020 Census, race on the "alone" basis plus Hispanic or Latino as a separate
#: ethnicity axis. Pre-fix value in the comment where it differed.
CENSUS_2020 = [
    ("white",             0.616),  # 0.576 (a non-Hispanic-alone figure)
    ("black",             0.124),  # 0.134
    ("hispanic",          0.187),  # 0.189
    ("asian",             0.060),  # 0.060 -- the one value already right
    ("native_american",   0.011),  # 0.013
    ("pacific_islander",  0.002),  # 0.003
    ("multiracial",       0.102),  # 0.025 -- understated FOUR-FOLD
]


@pytest.mark.parametrize("group,expected", CENSUS_2020)
def test_census_proportion_matches_the_2020_census(group, expected):
    """Each proportion equals the published 2020 Census share."""
    table = DiversityAssessmentWorkflow._US_CENSUS
    assert table[group] == pytest.approx(expected, abs=0.0005), (
        f"{group}: {table[group]} != {expected}"
    )


def test_multiracial_is_not_understated_fourfold():
    """0.025 against an actual 0.102.

    Two or More Races went from 9 million (2.9%) in 2010 to 33.8 million
    (10.2%) in 2020 -- the fastest-growing group in the US population.

    Consequence of the old value: a trial that enrolled 2.5% multiracial
    participants was scored as PERFECTLY representative and the workflow
    reported no gap, so a Diversity Action Plan is built on a target that
    under-represents that group four-fold -- while the agent's output is cited
    as conforming to FDA diversity guidance. The mirror case penalised a trial
    that enrolled 10% (essentially proportional) by scoring the 7.7-point
    "excess" as under-representation of everyone else.
    """
    assert DiversityAssessmentWorkflow._US_CENSUS["multiracial"] == pytest.approx(
        0.102, abs=0.0005
    )


def test_census_table_is_not_forced_to_sum_to_one():
    """Hispanic or Latino is an ethnicity that cuts across every race category.

    The old table was forced to sum to exactly 1.000, which no consistent
    Census basis does, and mixed a non-Hispanic-White-alone figure against
    alone-or-in-combination figures elsewhere.
    """
    total = sum(DiversityAssessmentWorkflow._US_CENSUS.values())
    assert total > 1.0, total


def test_census_table_names_its_source_and_basis():
    cite = DiversityAssessmentWorkflow._US_CENSUS_CITATION
    assert "2020 Census" in cite
    assert "alone" in cite
    assert "do not sum to 1" in cite


def test_pri_is_not_computed_against_the_general_population_for_an_indication():
    """A named indication must not be scored against the US population.

    ``reference_pct = epi.get("race_prevalence", self._US_CENSUS)`` silently
    substituted general-population proportions whenever indication
    epidemiology was absent. For sickle cell disease, Tay-Sachs or lupus the
    affected population looks nothing like the US population, so the PRI
    computed against it is meaningless -- and it was reported as a number
    anyway.
    """
    wf = DiversityAssessmentWorkflow()
    inputs = wf.preprocess({
        "indication": "sickle cell disease",
        "enrollment_data": {
            "total_enrolled": 100,
            "demographics": {"black": 90, "white": 10},
        },
    })
    result = wf.execute(inputs)
    joined = " ".join(result.findings)
    assert "not a valid reference for a named indication" in joined
    assert "Population Representativeness Index (PRI):" not in joined


def test_pri_is_computed_when_indication_epidemiology_is_supplied():
    """With a real reference the PRI is computed and the reference is named."""
    wf = DiversityAssessmentWorkflow()
    inputs = wf.preprocess({
        "indication": "sickle cell disease",
        "enrollment_data": {
            "total_enrolled": 100,
            "demographics": {"black": 90, "white": 10},
        },
        "indication_epidemiology": {
            "race_prevalence": {"black": 0.9, "white": 0.1},
        },
    })
    result = wf.execute(inputs)
    joined = " ".join(result.findings)
    assert "Population Representativeness Index (PRI):" in joined
    assert "indication-specific epidemiology" in joined
