"""The ACMG classification must be read from the interpretation, never defaulted into.

`/variants/interpret` used to scan the reply for the labels in the order
("Pathogenic", "Likely Pathogenic", "Likely Benign", "Benign") with a case-folded substring
test and take the first hit. Two consequences:

  * "This variant is not pathogenic" was classified **Pathogenic**.
  * "Likely Pathogenic" was classified **Pathogenic** — a different clinical call.

And because the prompt instructs the model to "Provide final classification: Pathogenic,
Likely Pathogenic, VUS, ...", almost every reply contained the word, so the endpoint was close
to hardcoded. The criteria block then asserted PVS1, PS1 and PM2 as met from that label alone.
"""
from api.routes.diagnostic_clinical import (
    _extract_cited_criteria,
    _parse_acmg_classification,
)


class TestNegationIsNotPathogenic:
    def test_plain_negation(self):
        assert _parse_acmg_classification("This variant is not pathogenic.") != "Pathogenic"

    def test_negation_returns_none_rather_than_guessing(self):
        assert _parse_acmg_classification("This variant is not pathogenic.") is None

    def test_non_hyphen_form(self):
        assert _parse_acmg_classification("Evidence suggests a non-pathogenic role.") is None


class TestLabelPrecision:
    def test_likely_pathogenic_is_not_collapsed(self):
        assert _parse_acmg_classification("Final classification: Likely Pathogenic") == "Likely Pathogenic"

    def test_likely_benign_is_not_collapsed(self):
        assert _parse_acmg_classification("Final classification: Likely Benign") == "Likely Benign"

    def test_explicit_pathogenic_still_reads(self):
        assert _parse_acmg_classification("Final classification: Pathogenic") == "Pathogenic"

    def test_vus(self):
        assert _parse_acmg_classification("Classification: VUS") == "VUS"

    def test_verdict_line_beats_an_earlier_mention(self):
        answer = (
            "Criteria for Pathogenic were considered but PM2 is weak.\n"
            "Final classification: VUS"
        )
        assert _parse_acmg_classification(answer) == "VUS"


class TestAmbiguityIsNotAVerdict:
    def test_prompt_echo_does_not_classify(self):
        """The prompt's own instruction, echoed back, is not a classification."""
        answer = (
            "I will classify per ACMG/AMP 2015. Options are Pathogenic, Likely Pathogenic, "
            "VUS, Likely Benign, or Benign. Further functional data is required."
        )
        assert _parse_acmg_classification(answer) is None

    def test_empty(self):
        assert _parse_acmg_classification("") is None


class TestCriteriaComeFromTheText:
    def test_cited_code_is_returned(self):
        crits = _extract_cited_criteria("PVS1 applies: this is a canonical splice variant.")
        assert [c.code for c in crits] == ["PVS1"]
        assert crits[0].met is True

    def test_not_met_is_respected(self):
        crits = _extract_cited_criteria("PM2 is not met; the variant is present in gnomAD.")
        assert crits[0].code == "PM2"
        assert crits[0].met is False

    def test_uncited_criteria_are_not_invented(self):
        """A Pathogenic verdict alone must not produce met criteria."""
        assert _extract_cited_criteria("Final classification: Pathogenic") == []

    def test_evidence_quotes_the_source_sentence(self):
        crits = _extract_cited_criteria("PS1 applies: same amino acid change as a known pathogenic variant.")
        assert "same amino acid change" in crits[0].evidence
