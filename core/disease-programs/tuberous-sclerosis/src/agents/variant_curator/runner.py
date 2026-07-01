"""
Agent 1 — TSC-Variant Curator (PRD §3 FR-VC-*; master paper §8). Act One centerpiece.

Four-step pipeline (the canonical tiering example):
  1. variant calling      — deterministic (Parabricks/GATK, mosaic-aware)   [W2 upgrade]
  2. annotation           — deterministic (snpEff/VEP + ClinVar/gnomAD/LOVD) [curated stub here]
  3. evidence aggregation — Sonnet (narrative + provenance)
  4. ACMG-AMP synthesis   — Opus proposes per-criterion reasoning; the deterministic
                            combinatorial validator (acmg.classify) is AUTHORITATIVE.

Recovers low-VAF somatic mosaic variants standard blood testing misses, flags mosaic
from VAF, recommends orthogonal validation (ddPCR), and produces a DRAFT for
board-certified molecular-geneticist review — never autonomous.
"""
from __future__ import annotations

from pathlib import Path

from config.settings import settings
from src.agents.base import Agent, AgentContext, AgentOutput
from src.agents.variant_curator.acmg import classify
from src.agents.variant_curator.annotation import annotate, assign_criteria
from src.agents.variant_curator.schema import (
    CriterionOut, VariantCuratorOutput, VariantInterpretation,
)
from src.cohort.genomic import count_filtered_artifacts, parse_variants
from src.utils.model_router import get_router
from src.utils.prompts import load_prompt

MOSAIC_VAF_CEILING = 0.30
_CONFIDENCE = {
    "Pathogenic": "high", "Likely Pathogenic": "moderate-high",
    "Variant of Uncertain Significance": "low",
    "Likely Benign": "moderate", "Benign": "high",
}


def _infer_consequence(variant_text: str) -> str:
    t = (variant_text or "").lower()
    if "frameshift" in t or "del" in t or "dup" in t or "ins" in t:
        return "frameshift"
    if "ter" in t or "*" in t or "nonsense" in t:
        return "nonsense"
    if "splice" in t:
        return "splice"
    return "missense"


def _variants_for(ctx: AgentContext) -> tuple[list[dict], str, int]:
    """Return (PASS variant records, source, count of strand-bias artifacts rejected)."""
    g = ctx.cohort.get("genomics", {})
    vcf_rel = g.get("vcf_path")
    if vcf_rel:
        path = Path(settings.COHORT_DIR) / vcf_rel
        if path.exists():
            return parse_variants(path), "vcf", count_filtered_artifacts(path)
    if g.get("gene"):  # inline fallback (e.g., dry_run_demo without a built cohort)
        parts = (g.get("variant") or "").split()
        hgvsc = parts[0] if parts else None
        return [{
            "gene": g["gene"], "hgvsc": hgvsc, "hgvsp": None,
            "consequence": _infer_consequence(g.get("variant", "")),
            "af": g.get("tissue_vaf") or g.get("vaf") or 0.5,
        }], "inline", 0
    return [], "none", 0


def _artifact_assessment(strand_balance: float | None) -> str:
    """Read-level QC: a balanced alt-strand distribution passes; a strongly skewed one is
    flagged. PASS records reaching the curator are balanced by construction (FR-VC-4)."""
    if strand_balance is None:
        return "not assessed"
    return "pass" if 0.2 <= strand_balance <= 0.8 else "strand-bias-suspected"


class VariantCurator(Agent):
    name = "variant_curator"

    def run(self, ctx: AgentContext) -> AgentOutput:
        router = get_router()
        variants, source, artifacts_filtered = _variants_for(ctx)
        phenotype_consistent = bool(ctx.cohort.get("phenotypes"))

        interps: list[VariantInterpretation] = []
        for v in variants:
            af = v.get("af")
            mosaic = af is not None and af < MOSAIC_VAF_CEILING and ctx.cohort.get("genomics", {}).get("tissue_vaf") is not None
            ann = annotate(v)
            crits = assign_criteria(ann, mosaic=mosaic, phenotype_consistent=phenotype_consistent)
            classification, rule = classify(crits)
            sb = v.get("strand_balance")
            interps.append(VariantInterpretation(
                gene=v.get("gene"), hgvsc=v.get("hgvsc"), hgvsp=v.get("hgvsp"),
                consequence=ann["consequence"] or None, vaf=af, mosaic=mosaic,
                recovered=mosaic,  # low-VAF tissue variant standard blood testing/pipelines miss
                acmg_classification=classification, acmg_rule=rule,
                acmg_criteria=[CriterionOut(code=c.code, bucket=c.bucket, rationale=c.rationale) for c in crits],
                confidence=_CONFIDENCE.get(classification, "low"),
                ddpcr_recommended=mosaic,
                depth=v.get("dp"), alt_reads=v.get("alt_reads"), strand_balance=sb,
                artifact_assessment=_artifact_assessment(sb),
            ))

        # rank: Pathogenic > Likely Pathogenic > VUS > ... ; mosaic recoveries first
        order = {"Pathogenic": 0, "Likely Pathogenic": 1, "Variant of Uncertain Significance": 2,
                 "Likely Benign": 3, "Benign": 4}
        interps.sort(key=lambda i: (order.get(i.acmg_classification, 9), not i.recovered))
        primary = interps[0] if interps else None

        # LLM steps: evidence aggregation (Sonnet) + per-criterion synthesis narrative (Opus).
        # Classification is the deterministic validator's; these add the narrative + provenance.
        _, prov_agg = router.call(
            self.name, "evidence_aggregation",
            system="Aggregate ClinVar/gnomAD/LOVD-TSC + literature evidence for the candidate variant.",
            prompt=f"{ctx.patient_id}: {primary.model_dump() if primary else 'no variant'}",
            prompt_template_version="evidence_aggregation/v0",
        )
        narrative, prov_acmg = router.call_structured(
            self.name, "acmg_synthesis",
            system=load_prompt("variant_curator", "acmg_synthesis"),
            prompt=f"criteria={[c.model_dump() for c in (primary.acmg_criteria if primary else [])]}, "
                   f"rule={primary.acmg_rule if primary else None}, "
                   f"classification={primary.acmg_classification if primary else None}",
            schema_hint={"per_criterion": [{"code": "PVS1", "reasoning": "..."}], "summary": "..."},
            prompt_template_version="acmg_synthesis/v0",
        )

        note = None
        if primary is None:
            note = ("No reportable variant on the available sample (e.g., negative blood test). "
                    "If clinical TSC is suspected, recommend disease-affected tissue sequencing "
                    "for somatic mosaic recovery (the no-mutation-identified gap).")

        out = VariantCuratorOutput(
            primary=primary, all_variants=interps, source=source, note=note,
            synthesis_narrative=narrative, artifacts_filtered=artifacts_filtered,
            draft_report_ref=f"tsc-reports/{ctx.patient_id}/variant_draft.md",
        )
        return self.ok(self.name, out.model_dump(), [prov_agg, prov_acmg])
