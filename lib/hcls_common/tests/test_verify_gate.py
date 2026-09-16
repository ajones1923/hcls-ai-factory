"""Tests for the HCLS verify-gate / honesty register (ORCH-11)."""
from hcls_common.verify_gate import (
    honesty_check, verify_text, is_publishable,
)


# ── deterministic honesty register (no LLM) ─────────────────────────────────
class TestHonestyRegister:
    def test_blocks_fda_clearance_claim(self):
        v = honesty_check("This assay is FDA-approved for diagnosing the variant.")
        assert any(x["severity"] == "block" for x in v)

    def test_blocks_cure_and_absolute_claims(self):
        assert any(x["severity"] == "block"
                   for x in honesty_check("This therapy cures the disease in 100% of patients."))

    def test_blocks_diagnostic_certainty(self):
        assert any(x["severity"] == "block"
                   for x in honesty_check("The result confirms the diagnosis of the disorder."))

    def test_blocks_efficacy_overclaim(self):
        # "clinically proven to treat/prevent" is an efficacy overclaim that must not slip through
        # (the bare verb "treat" evades the treat\w+ clinical trigger).
        assert any(x["severity"] == "block"
                   for x in honesty_check("Clinically proven to treat tuberous sclerosis."))
        assert any(x["severity"] == "block"
                   for x in honesty_check("This drug is proven to prevent seizures."))

    def test_warns_on_missing_disclaimer(self):
        v = honesty_check("The pathogenic variant suggests a treatment option for the patient.")
        assert any(x["severity"] == "warn" for x in v)
        assert not any(x["severity"] == "block" for x in v)

    def test_clean_text_with_disclaimer_passes(self):
        txt = ("For research use only and not a clinical diagnosis: the variant is "
               "classified likely pathogenic; clinician review required.")
        assert honesty_check(txt) == []


# ── combined gate verdicts ──────────────────────────────────────────────────
class TestVerifyText:
    def test_overclaim_blocks(self):
        v = verify_text("Our pipeline is FDA-approved and cures the disease.")
        assert v["status"] == "blocked"
        assert not is_publishable(v)

    def test_clean_verified(self):
        v = verify_text("For research use only: the candidate ranks highest by docking score.")
        assert v["status"] == "verified"
        assert is_publishable(v)

    def test_missing_disclaimer_is_draft(self):
        v = verify_text("The patient's variant is pathogenic and indicates a therapy.")
        assert v["status"] == "draft"
        assert not is_publishable(v)


# ── adversarial LLM verification (mocked client — no spend) ──────────────────
class _FakeLLM:
    """Returns canned JSON: 1 supported + 1 refuted (hallucinated) claim."""
    def generate_json(self, prompt, system_prompt="", max_tokens=1024):
        if "extract" in system_prompt.lower():
            return {"claims": ["CB-5083 inhibits VCP", "Docking score was -42.0 kcal/mol"]}
        return {"verdicts": [
            {"claim_index": 1, "verdict": "supported", "why": "evidence states CB-5083 is a VCP inhibitor"},
            {"claim_index": 2, "verdict": "refuted", "why": "evidence shows -10.5, not -42.0"},
        ]}


class TestAdversarialVerify:
    def test_refuted_claim_blocks(self):
        evidence = [{"text": "CB-5083 is an ATP-competitive VCP/p97 inhibitor; best docking -10.5 kcal/mol."}]
        v = verify_text(
            "CB-5083 inhibits VCP with a docking score of -42.0 kcal/mol. For research use only.",
            evidence=evidence, llm=_FakeLLM())
        assert v["summary"]["llm_checked"] is True
        assert v["summary"]["supported"] == 1
        assert v["summary"]["refuted"] == 1
        assert v["status"] == "blocked"          # the hallucinated stat is caught
        assert any(f["verdict"] == "refuted" for f in v["flagged"])
        assert not is_publishable(v)


# ── whose claim is it? (subject scoping, 2026-09-15) ─────────────────────────
class TestSubjectScoping:
    """A clearance claim about THIS platform must block; the same words about a named
    third-party therapy are ordinary clinical fact and must not be withheld.

    This distinction is what made it safe to enforce the gate by default. Without it,
    "Tisagenlecleucel is FDA-approved for paediatric r/r B-ALL" -- true, and exactly what a
    CAR-T agent exists to tell you -- would have been withheld as an overclaim.
    """

    def test_self_referential_clearance_blocks(self):
        for txt in ("This platform is FDA-approved for patients.",
                    "This assay is FDA-cleared for diagnosing the variant.",
                    "Our test is CE-marked for patient use."):
            assert any(x["severity"] == "block" for x in honesty_check(txt)), txt

    def test_third_party_clearance_is_a_warning_not_a_block(self):
        v = honesty_check("Tisagenlecleucel is FDA-approved for paediatric r/r B-ALL patients.")
        assert v, "should still be reported"
        assert not any(x["severity"] == "block" for x in v)
        assert any("third party" in x["message"] for x in v)

    def test_unsafe_regardless_of_subject_still_blocks(self):
        """A cure claim, absolute certainty and zero-risk are unsafe whoever the subject is."""
        for txt in ("Pembrolizumab cures the disease in patients.",
                    "This drug is certain to work for every patient.",
                    "The therapy carries zero risk for the patient."):
            assert any(x["severity"] == "block" for x in honesty_check(txt)), txt


# ── diagnostic certainty: about a patient, or about medicine? (2026-09-16) ────
class TestDiagnosticCertaintyScoping:
    """The flagship TSC question was WITHHELD on ~half of identical runs.

    "Which genes cause tuberous sclerosis and which pathway is dysregulated?" is a textbook
    genetics question, and a correct answer naturally says how the diagnosis is established.
    The rule fired on that and blocked the answer. A gate that withholds the right answer half
    the time teaches people to route around it.

    The distinction is the subject: a claim of certainty about A PATIENT (or about this
    platform's own output) is the thing the rule exists to stop. A statement about how diagnosis
    works in general is ordinary clinical education.
    """

    def test_patient_specific_certainty_still_blocks(self):
        for txt in ("This result confirms the diagnosis of tuberous sclerosis.",
                    "The patient's diagnosis is confirmed by the variant.",
                    "This establishes the diagnosis for this patient.",
                    "In this case the findings confirm the diagnosis.",
                    "Our analysis provides a definitive diagnosis."):
            v = honesty_check(txt)
            assert any(x["severity"] == "block" for x in v), txt

    def test_general_clinical_statements_degrade_to_warn(self):
        for txt in ("Genetic testing confirms the diagnosis in approximately 85% of cases; "
                    "for research use only, clinician review required.",
                    "A definitive diagnosis requires molecular confirmation of TSC1 or TSC2; "
                    "decision support only."):
            v = honesty_check(txt)
            assert not any(x["severity"] == "block" for x in v), txt
            assert any("general statement" in x["message"] for x in v), txt

    def test_the_flagship_answer_is_not_withheld(self):
        """The exact shape of answer that was being withheld.

        Asserts the contract the API gate actually applies — it withholds on a `block` finding or
        a refuted claim, and publishes a `warn` with the disclaimer attached. (`is_publishable`
        is a stricter, separate notion: it requires a clean `verified`.)
        """
        txt = ("Tuberous sclerosis complex is caused by loss-of-function variants in TSC1 or "
               "TSC2, which dysregulate the mTORC1 pathway. Genetic testing confirms the "
               "diagnosis in most patients. For research use only — decision support, not a "
               "diagnosis; clinician review required.")
        assert not any(x["severity"] == "block" for x in honesty_check(txt))
        assert verify_text(txt)["status"] != "blocked"

    def test_cure_and_absolutes_are_unaffected_by_scoping(self):
        """Scoping must not leak into the rules that are unsafe whatever the subject."""
        for txt in ("Everolimus cures tuberous sclerosis.",
                    "This therapy works in 100% of cases.",
                    "The treatment carries zero risk."):
            assert any(x["severity"] == "block" for x in honesty_check(txt)), txt
