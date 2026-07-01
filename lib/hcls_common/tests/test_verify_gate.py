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
