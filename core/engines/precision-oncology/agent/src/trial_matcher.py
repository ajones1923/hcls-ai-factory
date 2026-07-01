"""
Oncology Intelligence Agent - Trial Matcher
Deterministic + semantic clinical trial matching for oncology patients.

Combines hard-filter deterministic matching (cancer type, recruiting status)
with semantic vector search (biomarker profiles, stage context) to produce
ranked, explained trial recommendations for Molecular Tumor Board review.

Author: Adam Jones
Date: February 2026
License: Apache 2.0
"""

import logging
import re
from typing import Any, Dict, List, Optional

from src.models import CaseSnapshot

logger = logging.getLogger(__name__)

# Regex for validating filter values (prevents Milvus injection)
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _.\-/,]+$")


class TrialMatcher:
    """Matches oncology patients to clinical trials using a hybrid approach.

    Strategy:
      1. Deterministic filter on cancer_type (fuzzy) and open statuses
      2. Semantic search using embedded query of cancer_type + biomarkers + stage
      3. Composite scoring: biomarker match + phase weight + status weight
      4. Age-based eligibility filtering
      5. Structured explanation for each match
    """

    COLLECTION_NAME = "onco_trials"

    # Phase weighting: higher phase = more weight for ranking
    PHASE_WEIGHTS = {
        "Phase 3": 1.0,
        "Phase 2/3": 0.9,
        "Phase 2": 0.8,
        "Phase 1/2": 0.7,
        "Phase 1": 0.6,
        "Phase 4": 0.5,  # post-marketing, less relevant for novel treatment
    }

    # Recruiting status weighting
    STATUS_WEIGHTS = {
        "Recruiting": 1.0,
        "Active, not recruiting": 0.6,
        "Enrolling by invitation": 0.8,
        "Not yet recruiting": 0.4,
    }

    # Open statuses considered for deterministic search
    OPEN_STATUSES = {
        "Recruiting", "Active, not recruiting",
        "Enrolling by invitation", "Not yet recruiting",
    }

    # Cancer type aliases for fuzzy matching
    _CANCER_ALIASES = {
        "nsclc": ["non-small cell lung cancer", "nsclc", "non small cell lung", "lung adenocarcinoma", "lung squamous"],
        "breast": ["breast cancer", "breast", "triple negative breast", "tnbc", "her2+ breast"],
        "colorectal": ["colorectal cancer", "colorectal", "colon cancer", "colon", "crc", "rectal cancer"],
        "melanoma": ["melanoma", "cutaneous melanoma", "skin cancer"],
        "pancreatic": ["pancreatic cancer", "pancreatic", "pancreatic adenocarcinoma", "pdac"],
        "ovarian": ["ovarian cancer", "ovarian", "epithelial ovarian"],
        "prostate": ["prostate cancer", "prostate", "castration-resistant prostate", "crpc"],
        "glioblastoma": ["glioblastoma", "gbm", "glioblastoma multiforme", "glioma"],
        "aml": ["acute myeloid leukemia", "aml"],
        "cml": ["chronic myeloid leukemia", "cml"],
        "bladder": ["bladder cancer", "bladder", "urothelial carcinoma", "urothelial"],
        "renal": ["renal cell carcinoma", "renal", "kidney cancer", "rcc"],
        "gastric": ["gastric cancer", "gastric", "stomach cancer"],
        "hepatocellular": ["hepatocellular carcinoma", "hcc", "liver cancer"],
        "head_and_neck": ["head and neck", "hnscc", "head and neck squamous cell carcinoma"],
        "thyroid": ["thyroid cancer", "thyroid"],
        "esophageal": ["esophageal cancer", "esophageal"],
        "cholangiocarcinoma": ["cholangiocarcinoma", "bile duct cancer"],
        "sclc": ["small cell lung cancer", "sclc"],
    }

    def __init__(self, collection_manager, embedder):
        """Initialize the trial matcher.

        Args:
            collection_manager: Milvus collection manager for vector storage.
            embedder: Embedding model wrapper (e.g., BGE-small-en-v1.5).
        """
        self.collection_manager = collection_manager
        self.embedder = embedder

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def match_trials(
        self,
        cancer_type: str,
        biomarkers: Dict[str, Any],
        stage: str,
        age: Optional[int] = None,
        top_k: int = 10,
    ) -> List[Dict]:
        """Find and rank clinical trials matching the patient profile.

        Args:
            cancer_type: Cancer type (e.g. "NSCLC", "CRC", "Breast").
            biomarkers: Dict of biomarker results, e.g.
                {"EGFR": "L858R", "PD-L1_TPS": 80, "TMB": 14.2}.
            stage: Clinical stage (e.g. "IIIB", "IV").
            age: Optional patient age for eligibility filtering.
            top_k: Maximum number of trials to return.

        Returns:
            Ranked list of trial match dicts, each containing:
              trial_id, title, phase, match_score, matched_criteria,
              unmatched_criteria, explanation.
        """
        # Step 1: Deterministic filter -- cancer type + recruiting status
        deterministic_hits = self._deterministic_search(cancer_type, top_k=top_k * 3)
        logger.info(
            "Deterministic search for '%s': %d hits", cancer_type, len(deterministic_hits)
        )

        # Step 2: Semantic search -- embed patient profile query
        query_text = self._build_eligibility_query(cancer_type, biomarkers, stage)
        semantic_hits = self._semantic_search(query_text, top_k=top_k * 3)
        logger.info("Semantic search: %d hits", len(semantic_hits))

        # Merge results (union by trial_id, keep best score)
        merged = self._merge_results(deterministic_hits, semantic_hits)

        # Step 3: Score each trial
        patient_data = {
            "cancer_type": cancer_type,
            "biomarkers": biomarkers,
            "stage": stage,
            "age": age,
        }
        scored_trials = []
        for trial in merged:
            composite_score = self._compute_composite_score(trial, patient_data)
            trial["match_score"] = round(composite_score, 4)
            scored_trials.append(trial)

        # Sort by composite score descending
        scored_trials.sort(key=lambda t: t["match_score"], reverse=True)
        scored_trials = scored_trials[:top_k]

        # Step 4: Generate explanations
        results = []
        for trial in scored_trials:
            explained = self._explain_match(trial, patient_data, trial["match_score"])
            results.append(explained)

        logger.info(
            "Returning %d matched trials for %s stage %s", len(results), cancer_type, stage
        )
        return results

    def match_for_case(self, case: CaseSnapshot, top_k: int = 10) -> List[Dict]:
        """Convenience method to match trials for an existing CaseSnapshot.

        Extracts biomarkers from the case, including gene-level markers derived
        from the variant list.

        Args:
            case: A CaseSnapshot object.
            top_k: Maximum number of trials to return.

        Returns:
            Ranked list of trial match dicts.
        """
        # Combine explicit biomarkers with variant-derived markers
        biomarkers = dict(case.biomarkers) if case.biomarkers else {}
        for v in case.variants:
            gene = v.get("gene", "")
            if gene and v.get("actionability", "VUS") != "VUS":
                variant_str = v.get("variant", v.get("hgvs", ""))
                biomarkers[gene] = variant_str

        return self.match_trials(
            cancer_type=case.cancer_type,
            biomarkers=biomarkers,
            stage=case.stage,
            top_k=top_k,
        )

    # ------------------------------------------------------------------
    # Query building
    # ------------------------------------------------------------------

    def _build_eligibility_query(
        self, cancer_type: str, biomarkers: Dict[str, Any], stage: str
    ) -> str:
        """Create a natural-language search query from patient data.

        Args:
            cancer_type: Cancer type string.
            biomarkers: Biomarker dict.
            stage: Clinical stage.

        Returns:
            Query string suitable for embedding.
        """
        parts = [f"{cancer_type} clinical trial"]

        if stage:
            parts.append(f"stage {stage}")

        # Add biomarker context
        for marker, value in biomarkers.items():
            if isinstance(value, (int, float)):
                parts.append(f"{marker} {value}")
            else:
                parts.append(f"{marker} {value}")

        return " ".join(parts)

    # ------------------------------------------------------------------
    # Search methods
    # ------------------------------------------------------------------

    def _resolve_cancer_aliases(self, cancer_type: str) -> List[str]:
        """Resolve a cancer type input to a list of known aliases for fuzzy matching."""
        ct_lower = cancer_type.strip().lower()
        # Check if input matches any alias group
        for _canonical, aliases in self._CANCER_ALIASES.items():
            if ct_lower in aliases or ct_lower == _canonical:
                return aliases
        # Fall back to the input itself
        return [ct_lower]

    def _deterministic_search(self, cancer_type: str, top_k: int = 30) -> List[Dict]:
        """Search onco_trials with deterministic filters.

        Uses fuzzy cancer type matching via aliases and searches across
        all open trial statuses (Recruiting, Active, Enrolling by invitation,
        Not yet recruiting).
        """
        all_results = []
        cancer_aliases = self._resolve_cancer_aliases(cancer_type)
        status_list = list(self.OPEN_STATUSES)

        for alias in cancer_aliases:
            for status in status_list:
                try:
                    if not _SAFE_FILTER_RE.match(alias):
                        logger.warning("Rejected unsafe cancer alias: %r", alias)
                        continue
                    if not _SAFE_FILTER_RE.match(status):
                        logger.warning("Rejected unsafe status filter: %r", status)
                        continue
                    filter_expr = (
                        f'cancer_type == "{alias}" and status == "{status}"'
                    )
                    results = self.collection_manager.query(
                        collection_name=self.COLLECTION_NAME,
                        filter_expr=filter_expr,
                        output_fields=[
                            "trial_id", "title", "phase", "status", "cancer_type",
                            "criteria", "text", "biomarker_criteria", "sponsor",
                        ],
                        limit=top_k,
                    )
                    all_results.extend([dict(r) for r in results])
                except Exception as exc:
                    logger.debug("Deterministic search for '%s'/'%s' failed: %s", alias, status, exc)

        # Deduplicate by trial_id
        seen = {}
        for r in all_results:
            tid = r.get("trial_id", "")
            if tid and tid not in seen:
                seen[tid] = r
        return list(seen.values())[:top_k]

    def _semantic_search(self, query_text: str, top_k: int = 30) -> List[Dict]:
        """Search onco_trials using semantic embedding similarity."""
        try:
            embedding = self.embedder.embed(query_text)
            results = self.collection_manager.search(
                collection_name=self.COLLECTION_NAME,
                query_vector=embedding,
                top_k=top_k,
                output_fields=[
                    "trial_id", "title", "phase", "status", "cancer_type",
                    "criteria", "text", "biomarker_criteria", "sponsor",
                ],
            )
            return [dict(r) for r in results]
        except Exception as exc:
            logger.warning("Semantic trial search failed: %s", exc)
            return []

    def _merge_results(
        self, deterministic: List[Dict], semantic: List[Dict]
    ) -> List[Dict]:
        """Merge deterministic and semantic results, deduplicating by trial_id."""
        seen = {}
        for trial in deterministic:
            tid = trial.get("trial_id", "")
            if tid:
                trial["_source"] = "deterministic"
                seen[tid] = trial

        for trial in semantic:
            tid = trial.get("trial_id", "")
            if tid and tid not in seen:
                trial["_source"] = "semantic"
                seen[tid] = trial
            elif tid and tid in seen:
                # Keep the semantic score if present
                if "score" in trial:
                    seen[tid]["_semantic_score"] = trial["score"]

        return list(seen.values())

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def _compute_composite_score(self, trial: Dict, patient_data: Dict) -> float:
        """Compute a composite match score for a trial.

        Components:
          - biomarker_match_score (0-1): fraction of patient biomarkers
            mentioned in trial criteria
          - phase_weight (0-1): Phase 3 > Phase 2 > Phase 1
          - status_weight (0-1): Recruiting > Active > Not yet recruiting
          - semantic_score (0-1): if available from vector search
        """
        criteria_text = trial.get("criteria", "") + " " + trial.get("biomarker_criteria", "")

        biomarker_score = self._score_biomarker_match(
            criteria_text, patient_data.get("biomarkers", {})
        )

        phase = trial.get("phase", "")
        phase_weight = self.PHASE_WEIGHTS.get(phase, 0.5)

        status = trial.get("status", "")
        status_weight = self.STATUS_WEIGHTS.get(status, 0.3)

        semantic_score = trial.get("_semantic_score", trial.get("score", 0.5))

        # Age penalty: reduce score if patient age is outside trial criteria
        age_penalty = self._compute_age_penalty(
            criteria_text, patient_data.get("age")
        )

        # Weighted combination
        composite = (
            0.40 * biomarker_score
            + 0.25 * semantic_score
            + 0.20 * phase_weight
            + 0.15 * status_weight
        )
        composite *= age_penalty
        return composite

    def _compute_age_penalty(
        self, criteria_text: str, age: Optional[int]
    ) -> float:
        """Apply a scoring penalty if patient age falls outside trial criteria.

        Looks for patterns like "Age >= 18" or "Age <= 75" or "18-75 years"
        in the criteria text.

        Returns:
            1.0 (no penalty) if age is within range or not specified,
            0.5 if age is out of range.
        """
        if age is None or not criteria_text:
            return 1.0

        criteria_lower = criteria_text.lower()

        # Pattern: "age >= 18", "age ≥ 18", "minimum age: 18"
        min_age_match = re.search(r"(?:age\s*(?:>=|≥|>)\s*(\d+))|(?:minimum\s*age[:\s]*(\d+))", criteria_lower)
        max_age_match = re.search(r"(?:age\s*(?:<=|≤|<)\s*(\d+))|(?:maximum\s*age[:\s]*(\d+))", criteria_lower)

        # Pattern: "18-75 years" or "18 to 75 years"
        range_match = re.search(r"(\d+)\s*(?:-|to)\s*(\d+)\s*years", criteria_lower)

        min_age = None
        max_age = None

        if min_age_match:
            min_age = int(min_age_match.group(1) or min_age_match.group(2))
        if max_age_match:
            max_age = int(max_age_match.group(1) or max_age_match.group(2))
        if range_match:
            min_age = min_age or int(range_match.group(1))
            max_age = max_age or int(range_match.group(2))

        if min_age is not None and age < min_age:
            return 0.5
        if max_age is not None and age > max_age:
            return 0.5

        return 1.0

    def _score_biomarker_match(
        self, trial_biomarker_criteria: str, patient_biomarkers: Dict[str, Any]
    ) -> float:
        """Fuzzy matching of patient biomarkers against trial criteria text.

        For each patient biomarker key/value pair, checks if it appears in
        the trial criteria text (case-insensitive). Returns fraction matched.

        Args:
            trial_biomarker_criteria: Combined criteria text from the trial.
            patient_biomarkers: Patient biomarker dict.

        Returns:
            Float between 0.0 and 1.0 representing match fraction.
        """
        if not patient_biomarkers or not trial_biomarker_criteria:
            return 0.0

        criteria_lower = trial_biomarker_criteria.lower()
        matches = 0
        total = len(patient_biomarkers)

        for marker, value in patient_biomarkers.items():
            marker_lower = marker.lower().replace("_", " ").replace("-", " ")
            value_str = str(value).lower()

            # Check if marker name appears in criteria
            if marker_lower in criteria_lower or value_str in criteria_lower:
                matches += 1
            # Also check common aliases
            elif marker_lower.replace(" ", "-") in criteria_lower:
                matches += 1

        return matches / total if total > 0 else 0.0

    # ------------------------------------------------------------------
    # Explanation generation
    # ------------------------------------------------------------------

    def _explain_match(self, trial: Dict, patient_data: Dict, score: float) -> Dict:
        """Generate a structured explanation for why a trial matches.

        Args:
            trial: Trial dict from search results.
            patient_data: Patient data dict with cancer_type, biomarkers, stage.
            score: Composite match score.

        Returns:
            Dict with trial details and structured explanation.
        """
        criteria_text = trial.get("criteria", "") + " " + trial.get("biomarker_criteria", "")
        criteria_lower = criteria_text.lower()

        matched_criteria = []
        unmatched_criteria = []

        biomarkers = patient_data.get("biomarkers", {})
        for marker, value in biomarkers.items():
            marker_str = f"{marker}={value}"
            if marker.lower() in criteria_lower or str(value).lower() in criteria_lower:
                matched_criteria.append(marker_str)
            else:
                unmatched_criteria.append(marker_str)

        # Check cancer type match
        cancer_type = patient_data.get("cancer_type", "")
        if cancer_type.lower() in criteria_lower or cancer_type.lower() in trial.get("cancer_type", "").lower():
            matched_criteria.insert(0, f"Cancer type: {cancer_type}")
        else:
            unmatched_criteria.insert(0, f"Cancer type: {cancer_type} (not explicitly listed)")

        # Check stage
        stage = patient_data.get("stage", "")
        if stage and stage.lower() in criteria_lower:
            matched_criteria.append(f"Stage {stage}")

        # Check age eligibility
        age = patient_data.get("age")
        if age is not None:
            age_penalty = self._compute_age_penalty(criteria_text, age)
            if age_penalty < 1.0:
                unmatched_criteria.append(f"Age {age} (may be outside eligibility range)")
            else:
                matched_criteria.append(f"Age {age}")

        # Build explanation text
        explanation_parts = []
        if matched_criteria:
            explanation_parts.append(f"Matched: {', '.join(matched_criteria)}.")
        if unmatched_criteria:
            explanation_parts.append(f"Not confirmed: {', '.join(unmatched_criteria)}.")

        phase = trial.get("phase", "Unknown")
        explanation_parts.append(f"Trial is {phase}, status: {trial.get('status', 'Unknown')}.")

        return {
            "trial_id": trial.get("trial_id", ""),
            "title": trial.get("title", ""),
            "phase": phase,
            "status": trial.get("status", ""),
            "sponsor": trial.get("sponsor", ""),
            "match_score": score,
            "matched_criteria": matched_criteria,
            "unmatched_criteria": unmatched_criteria,
            "explanation": " ".join(explanation_parts),
        }
