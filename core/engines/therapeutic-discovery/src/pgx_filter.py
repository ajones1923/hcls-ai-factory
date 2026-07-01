"""
Pharmacogenomic Drug Candidate Filter.

Receives PGx results from the Precision Biomarker Agent and filters/re-ranks
drug candidates based on the patient's metabolizer phenotype. This closes
the loop between biomarker intelligence and drug discovery.

The key insight: it's not enough to find molecules that bind — we need
molecules that THIS SPECIFIC PATIENT can metabolize safely.

Part of the HCLS AI Factory: Patient DNA → Drug Candidates pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, unique
from typing import Any

from loguru import logger


@unique
class MetabolismRisk(str, Enum):
    """Risk level for a drug candidate based on PGx profile."""
    SAFE = "safe"
    CAUTION = "caution"
    CONTRAINDICATED = "contraindicated"
    UNKNOWN = "unknown"


# Mapping of CYP enzymes to their known substrates and risk implications
# Based on CPIC guidelines and FDA drug labels
CYP_SUBSTRATE_MAP: dict[str, dict[str, list[str]]] = {
    "CYP2D6": {
        "substrates": ["codeine", "tramadol", "tamoxifen", "fluoxetine", "paroxetine",
                       "venlafaxine", "atomoxetine", "dextromethorphan", "metoprolol",
                       "propafenone", "flecainide", "ondansetron"],
        "prodrugs": ["codeine", "tramadol", "tamoxifen"],  # Require activation
    },
    "CYP2C19": {
        "substrates": ["clopidogrel", "omeprazole", "escitalopram", "voriconazole",
                       "diazepam", "phenytoin", "proguanil"],
        "prodrugs": ["clopidogrel", "proguanil"],
    },
    "CYP2C9": {
        "substrates": ["warfarin", "phenytoin", "celecoxib", "flurbiprofen",
                       "losartan", "glipizide", "tolbutamide"],
        "prodrugs": ["losartan"],
    },
    "CYP3A5": {
        "substrates": ["tacrolimus", "cyclosporine", "sirolimus"],
        "prodrugs": [],
    },
    "DPYD": {
        "substrates": ["fluorouracil", "capecitabine", "tegafur"],
        "prodrugs": ["capecitabine"],
    },
    "TPMT": {
        "substrates": ["azathioprine", "mercaptopurine", "thioguanine"],
        "prodrugs": [],
    },
    "SLCO1B1": {
        "substrates": ["simvastatin", "atorvastatin", "rosuvastatin", "pravastatin",
                       "pitavastatin", "methotrexate"],
        "prodrugs": [],
    },
}

# Phenotype risk mapping
PHENOTYPE_RISK: dict[str, dict[str, MetabolismRisk]] = {
    "poor_metabolizer": {
        "substrate": MetabolismRisk.CONTRAINDICATED,  # Drug accumulation
        "prodrug": MetabolismRisk.CONTRAINDICATED,     # No activation
    },
    "intermediate_metabolizer": {
        "substrate": MetabolismRisk.CAUTION,
        "prodrug": MetabolismRisk.CAUTION,
    },
    "normal_metabolizer": {
        "substrate": MetabolismRisk.SAFE,
        "prodrug": MetabolismRisk.SAFE,
    },
    "rapid_metabolizer": {
        "substrate": MetabolismRisk.SAFE,
        "prodrug": MetabolismRisk.SAFE,
    },
    "ultra_rapid_metabolizer": {
        "substrate": MetabolismRisk.CAUTION,    # Too fast clearance
        "prodrug": MetabolismRisk.CONTRAINDICATED,  # Toxic activation rate
    },
}


@dataclass
class PGxFilterResult:
    """Result of PGx filtering for a single drug candidate."""
    candidate_id: str
    original_rank: int
    adjusted_rank: int
    metabolism_risk: MetabolismRisk
    affected_enzymes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    recommendation: str = ""
    rank_adjustment: int = 0  # Positive = demoted, negative = promoted


@dataclass
class PGxFilterSummary:
    """Summary of PGx filtering across all candidates."""
    total_candidates: int = 0
    safe_count: int = 0
    caution_count: int = 0
    contraindicated_count: int = 0
    unknown_count: int = 0
    results: list[PGxFilterResult] = field(default_factory=list)
    patient_pgx_profile: dict[str, str] = field(default_factory=dict)


class PGxDrugFilter:
    """Filter and re-rank drug candidates based on pharmacogenomic profile.

    This is the bridge between the Precision Biomarker Agent's PGx analysis
    and the Drug Discovery Pipeline's candidate ranking. It implements the
    "closed loop" architecture where drug candidates are personalized to
    the patient's genome.

    Usage:
        filter = PGxDrugFilter(pgx_profile)
        summary = filter.filter_candidates(ranked_candidates)
        # summary.results contains re-ranked candidates with PGx annotations
    """

    def __init__(self, pgx_profile: dict[str, Any] | None = None):
        """
        Parameters
        ----------
        pgx_profile : dict
            PGx results from Biomarker Agent. Expected format:
            {
                "gene_results": [
                    {"gene": "CYP2D6", "phenotype": "Poor Metabolizer", ...},
                    {"gene": "CYP2C19", "phenotype": "Normal Metabolizer", ...},
                ]
            }
        """
        self.pgx_profile = pgx_profile or {}
        self._phenotype_map = self._build_phenotype_map()

    def _build_phenotype_map(self) -> dict[str, str]:
        """Build gene -> normalized phenotype mapping from profile."""
        result = {}
        for gene_result in self.pgx_profile.get("gene_results", []):
            gene = gene_result.get("gene", "")
            phenotype = gene_result.get("phenotype", "")
            if gene and phenotype:
                normalized = phenotype.lower().replace(" ", "_").replace("-", "_")
                result[gene] = normalized
        return result

    def assess_candidate(
        self,
        candidate: dict[str, Any],
        rank: int,
    ) -> PGxFilterResult:
        """Assess a single drug candidate against the PGx profile.

        Parameters
        ----------
        candidate : dict
            Drug candidate with at minimum a 'smiles' or 'name' field.
            May also include 'target_pathway', 'metabolism_enzymes'.
        rank : int
            Current ranking position (1 = best).
        """
        candidate_id = candidate.get("name", candidate.get("smiles", f"candidate_{rank}"))
        result = PGxFilterResult(
            candidate_id=str(candidate_id),
            original_rank=rank,
            adjusted_rank=rank,
            metabolism_risk=MetabolismRisk.UNKNOWN,
        )

        # Check known metabolism enzymes if specified
        metabolism_enzymes = candidate.get("metabolism_enzymes", [])

        worst_risk = MetabolismRisk.SAFE

        for gene, phenotype in self._phenotype_map.items():
            if gene not in CYP_SUBSTRATE_MAP:
                continue

            gene_info = CYP_SUBSTRATE_MAP[gene]
            risk_map = PHENOTYPE_RISK.get(phenotype, {})

            # Check if candidate is likely metabolized by this enzyme
            # For novel candidates, we flag based on the enzyme's known substrate class
            is_substrate = gene in metabolism_enzymes
            is_prodrug = candidate.get("is_prodrug", False)

            if is_substrate or not metabolism_enzymes:
                # If no specific enzyme info, we conservatively flag all
                # poor/ultra-rapid metabolizers
                if phenotype in ("poor_metabolizer", "ultra_rapid_metabolizer"):
                    drug_type = "prodrug" if is_prodrug else "substrate"
                    risk = risk_map.get(drug_type, MetabolismRisk.UNKNOWN)

                    if risk.value > worst_risk.value if hasattr(risk, 'value') else False:
                        worst_risk = risk

                    result.affected_enzymes.append(gene)
                    result.warnings.append(
                        f"{gene} {phenotype.replace('_', ' ')}: "
                        f"{'avoid' if risk == MetabolismRisk.CONTRAINDICATED else 'use caution'}"
                    )

        # Determine overall risk
        if result.affected_enzymes:
            # Use the worst risk found
            risks = []
            for gene in result.affected_enzymes:
                phenotype = self._phenotype_map.get(gene, "")
                risk_map = PHENOTYPE_RISK.get(phenotype, {})
                drug_type = "prodrug" if candidate.get("is_prodrug", False) else "substrate"
                risk = risk_map.get(drug_type, MetabolismRisk.UNKNOWN)
                risks.append(risk)

            if MetabolismRisk.CONTRAINDICATED in risks:
                result.metabolism_risk = MetabolismRisk.CONTRAINDICATED
                result.rank_adjustment = 100  # Heavily demote
                result.recommendation = "AVOID: Patient metabolizer phenotype contraindicates this compound"
            elif MetabolismRisk.CAUTION in risks:
                result.metabolism_risk = MetabolismRisk.CAUTION
                result.rank_adjustment = 10  # Moderate demotion
                result.recommendation = "CAUTION: Dose adjustment may be needed based on PGx profile"
            else:
                result.metabolism_risk = MetabolismRisk.SAFE
                result.recommendation = "PGx-compatible: no known metabolic concerns"
        else:
            if self._phenotype_map:
                result.metabolism_risk = MetabolismRisk.SAFE
                result.recommendation = "PGx-compatible: no known metabolic pathway conflicts"
            else:
                result.metabolism_risk = MetabolismRisk.UNKNOWN
                result.recommendation = "PGx data not available: standard dosing assumed"

        result.adjusted_rank = rank + result.rank_adjustment
        return result

    def filter_candidates(
        self,
        candidates: list[dict[str, Any]],
    ) -> PGxFilterSummary:
        """Filter and re-rank all drug candidates.

        Returns candidates re-ranked with PGx-safe compounds promoted
        and contraindicated compounds demoted.
        """
        summary = PGxFilterSummary(
            total_candidates=len(candidates),
            patient_pgx_profile=dict(self._phenotype_map),
        )

        for i, candidate in enumerate(candidates):
            result = self.assess_candidate(candidate, rank=i + 1)
            summary.results.append(result)

        # Re-rank by adjusted rank
        summary.results.sort(key=lambda r: r.adjusted_rank)
        for i, result in enumerate(summary.results):
            result.adjusted_rank = i + 1

        # Count by risk level
        for result in summary.results:
            if result.metabolism_risk == MetabolismRisk.SAFE:
                summary.safe_count += 1
            elif result.metabolism_risk == MetabolismRisk.CAUTION:
                summary.caution_count += 1
            elif result.metabolism_risk == MetabolismRisk.CONTRAINDICATED:
                summary.contraindicated_count += 1
            else:
                summary.unknown_count += 1

        logger.info(
            f"PGx filter: {summary.safe_count} safe, {summary.caution_count} caution, "
            f"{summary.contraindicated_count} contraindicated out of {summary.total_candidates}"
        )

        return summary
