"""Radiomics feature ingest pipeline for Imaging Intelligence Agent.

Takes RadiomicsEngine output, generates text summaries, embeds via
BGE-small-en-v1.5, and stores in the imaging_radiomics Milvus collection.

Follows the same BaseIngestPipeline pattern as:
  - src/ingest/literature_parser.py (PubMedImagingIngestPipeline)
  - src/ingest/finding_parser.py (FindingIngestPipeline)

Author: Adam Jones
Date: April 2026
"""

import hashlib
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field

from src.radiomics_engine import RadiomicsEngine

from .base import BaseIngestPipeline


# ═══════════════════════════════════════════════════════════════════════
# Pydantic Model for Radiomics Records
# ═══════════════════════════════════════════════════════════════════════


class RadiomicsRecord(BaseModel):
    """A single radiomics feature record for Milvus storage.

    Maps to the imaging_radiomics collection schema defined in
    src/collections.py (RADIOMICS_FIELDS).
    """
    id: str = Field(..., max_length=200, description="Unique record identifier")
    patient_id: str = Field("", max_length=50)
    study_date: str = Field("", max_length=20, description="YYYY-MM-DD")
    region_label: str = Field("", max_length=100, description="Anatomical region name")
    modality: str = Field("", max_length=10, description="ct, mri, etc.")
    feature_class: str = Field("", max_length=50, description="firstorder, shape, glcm, etc.")
    feature_count: int = Field(0, ge=0)
    mean_intensity: float = Field(0.0)
    volume_mm3: float = Field(0.0)
    sphericity: float = Field(0.0)
    entropy: float = Field(0.0)
    feature_summary: str = Field("", max_length=5000, description="Text summary for RAG search")

    def to_embedding_text(self) -> str:
        """Generate text for BGE-small-en-v1.5 embedding.

        Combines the feature summary with structured metadata to produce
        a rich text representation suitable for semantic search.
        """
        parts = []
        if self.feature_summary:
            parts.append(self.feature_summary)
        if self.region_label:
            parts.append(f"Region: {self.region_label}")
        if self.modality:
            parts.append(f"Modality: {self.modality}")
        if self.feature_class:
            parts.append(f"Feature class: {self.feature_class}")
        if self.volume_mm3 > 0:
            parts.append(f"Volume: {self.volume_mm3:.1f} mm3")
        if self.mean_intensity != 0.0:
            parts.append(f"Mean intensity: {self.mean_intensity:.1f}")
        return " ".join(parts) if parts else "Radiomics feature extraction result"


# ═══════════════════════════════════════════════════════════════════════
# Radiomics Ingest Pipeline
# ═══════════════════════════════════════════════════════════════════════


class RadiomicsIngestPipeline(BaseIngestPipeline):
    """Ingest pipeline for radiomics feature data.

    Takes RadiomicsEngine extraction results, converts them into
    RadiomicsRecord models, and stores embeddings in the
    imaging_radiomics Milvus collection.

    Unlike other ingest pipelines that fetch from external APIs, this
    pipeline receives pre-computed feature data from the RadiomicsEngine.
    The fetch() method accepts feature dicts directly.

    Usage:
        engine = RadiomicsEngine(feature_classes=["firstorder", "shape", "glcm"])
        features = engine.extract_all_regions("image.nii.gz", "mask.nii.gz")

        pipeline = RadiomicsIngestPipeline(
            collection_manager=manager,
            embedder=embedder,
            radiomics_engine=engine,
        )
        count = pipeline.run(
            image_path="image.nii.gz",
            mask_path="mask.nii.gz",
            patient_id="PATIENT-001",
            modality="ct",
        )
    """

    COLLECTION_NAME = "imaging_radiomics"

    def __init__(
        self,
        collection_manager,
        embedder,
        radiomics_engine: Optional[RadiomicsEngine] = None,
    ):
        """Initialize the radiomics ingest pipeline.

        Args:
            collection_manager: ImagingCollectionManager for Milvus operations.
            embedder: Embedding model with encode() method (BGE-small-en-v1.5).
            radiomics_engine: RadiomicsEngine instance.  If None, a default
                mock-mode engine is created.
        """
        super().__init__(collection_manager, embedder)
        self.radiomics_engine = radiomics_engine or RadiomicsEngine(mock=True)

    def fetch(
        self,
        image_path: str = "",
        mask_path: str = "",
        patient_id: str = "",
        modality: str = "ct",
        study_date: str = "",
        region_names: Optional[Dict[int, str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        """Extract radiomics features from image and mask.

        Args:
            image_path: Path to NIfTI image file.
            mask_path: Path to NIfTI segmentation mask file.
            patient_id: Patient identifier for the record.
            modality: Imaging modality (e.g. "ct", "mri").
            study_date: Study date string (YYYY-MM-DD).  Defaults to today.
            region_names: Optional mapping of label int -> region name string.
            **kwargs: Additional fetch parameters (ignored).

        Returns:
            Dict with keys: regions (feature dicts per label), patient_id,
            modality, study_date, region_names.
        """
        if not study_date:
            study_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")

        logger.info(
            f"Fetching radiomics features: patient={patient_id}, "
            f"modality={modality}, image={image_path}"
        )

        # Extract features for all regions
        regions = self.radiomics_engine.extract_all_regions(image_path, mask_path)

        return {
            "regions": regions,
            "patient_id": patient_id,
            "modality": modality,
            "study_date": study_date,
            "region_names": region_names or {},
        }

    def parse(self, raw_data: Dict[str, Any]) -> List[RadiomicsRecord]:
        """Parse radiomics extraction results into RadiomicsRecord models.

        Creates one record per feature class per region, so that each
        record has a focused text summary suitable for embedding.

        Args:
            raw_data: Output from fetch() containing regions and metadata.

        Returns:
            List of validated RadiomicsRecord instances.
        """
        regions: Dict[int, Dict[str, float]] = raw_data.get("regions", {})
        patient_id: str = raw_data.get("patient_id", "")
        modality: str = raw_data.get("modality", "ct")
        study_date: str = raw_data.get("study_date", "")
        region_names: Dict[int, str] = raw_data.get("region_names", {})

        records: List[RadiomicsRecord] = []

        for label, features in regions.items():
            region_name = region_names.get(label, f"region_{label}")

            # Group features by class prefix
            class_groups: Dict[str, Dict[str, float]] = {}
            for fname, fval in features.items():
                parts = fname.split("_", 1)
                if len(parts) == 2:
                    cls_name = parts[0].lower()
                else:
                    cls_name = "other"
                class_groups.setdefault(cls_name, {})[fname] = fval

            # Create one record per feature class
            for feat_class, class_features in class_groups.items():
                record_id = self._generate_id(
                    patient_id, study_date, region_name, feat_class
                )

                # Extract key scalar values
                mean_intensity = features.get("firstorder_Mean", 0.0)
                volume_mm3 = features.get("shape_VoxelVolume", 0.0)
                sphericity = features.get("shape_Sphericity", 0.0)
                entropy = features.get("firstorder_Entropy", 0.0)

                # Generate text summary
                summary = RadiomicsEngine.generate_feature_summary(
                    features=class_features,
                    region_name=region_name,
                    modality=modality,
                )

                try:
                    record = RadiomicsRecord(
                        id=record_id[:200],
                        patient_id=patient_id[:50],
                        study_date=study_date[:20],
                        region_label=region_name[:100],
                        modality=modality[:10],
                        feature_class=feat_class[:50],
                        feature_count=len(class_features),
                        mean_intensity=round(mean_intensity, 4),
                        volume_mm3=round(volume_mm3, 4),
                        sphericity=round(sphericity, 4),
                        entropy=round(entropy, 4),
                        feature_summary=summary[:5000],
                    )
                    records.append(record)
                except Exception as e:
                    logger.warning(
                        f"Failed to create RadiomicsRecord for "
                        f"label={label}, class={feat_class}: {e}"
                    )
                    continue

        logger.info(f"Parsed {len(records)} radiomics records from {len(regions)} regions")
        return records

    @staticmethod
    def _generate_id(
        patient_id: str,
        study_date: str,
        region_name: str,
        feature_class: str,
    ) -> str:
        """Generate a deterministic unique ID for a radiomics record.

        Args:
            patient_id: Patient identifier.
            study_date: Study date string.
            region_name: Region name string.
            feature_class: Feature class name.

        Returns:
            Deterministic ID string (max 200 chars).
        """
        raw = f"{patient_id}|{study_date}|{region_name}|{feature_class}"
        digest = hashlib.sha256(raw.encode()).hexdigest()[:16]
        prefix = f"rad-{patient_id[:20]}-{region_name[:20]}-{feature_class[:10]}"
        return f"{prefix}-{digest}"

    def run(
        self,
        collection_name: Optional[str] = None,
        batch_size: int = 32,
        **fetch_kwargs,
    ) -> int:
        """Execute the full radiomics ingest pipeline.

        Args:
            collection_name: Target collection (defaults to 'imaging_radiomics').
            batch_size: Batch size for embedding and insertion.
            **fetch_kwargs: Passed to fetch() (image_path, mask_path,
                patient_id, modality, study_date, region_names).

        Returns:
            Total number of records ingested.
        """
        target = collection_name or self.COLLECTION_NAME
        logger.info(f"Starting radiomics ingest pipeline -> {target}")

        raw = self.fetch(**fetch_kwargs)
        logger.info(
            f"Extracted features for {len(raw.get('regions', {}))} regions"
        )

        records = self.parse(raw)
        logger.info(f"Parsed {len(records)} RadiomicsRecord instances")

        count = self.embed_and_store(records, target, batch_size)
        logger.info(f"Ingested {count} radiomics records into {target}")
        return count
