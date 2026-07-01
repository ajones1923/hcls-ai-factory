"""Base MONAI Deploy Application for imaging workflows.

Uses MONAI Deploy App SDK to create portable DICOM-aware containers
that follow clinical deployment standards. Each MAP:
1. Receives DICOM studies via standard DICOM input
2. Selects appropriate series based on configurable rules
3. Preprocesses and runs inference (via existing workflow classes)
4. Generates DICOM SR output with structured findings
5. Optionally generates FHIR R4 DiagnosticReport

Apache 2.0 Licensed.
"""

import json
import tempfile
import time
from abc import abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from loguru import logger

# ── MONAI Deploy App SDK (graceful degradation) ─────────────────────────
try:
    from monai.deploy.core import Application, resource, env, IOType
    from monai.deploy.core import InputContext, OutputContext, ExecutionContext
    from monai.deploy.operators import (
        DICOMDataLoaderOperator,
        DICOMSeriesSelectorOperator,
    )

    MONAI_DEPLOY_AVAILABLE = True
except ImportError:
    MONAI_DEPLOY_AVAILABLE = False

    # Stub classes so the module is importable without monai-deploy-app-sdk
    class Application:  # type: ignore[no-redef]
        """Stub Application when monai-deploy-app-sdk is not installed."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

        def compose(self) -> None:
            pass

    class resource:  # type: ignore[no-redef]
        """Stub decorator."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

        def __call__(self, cls: Any) -> Any:
            return cls

    class env:  # type: ignore[no-redef]
        """Stub decorator."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

        def __call__(self, cls: Any) -> Any:
            return cls

    logger.warning(
        "monai-deploy-app-sdk not installed. MAP classes will operate in "
        "standalone mode. Install with: pip install monai-deploy-app-sdk>=0.6.0"
    )

# ── Medical imaging libraries (graceful degradation) ────────────────────
try:
    import pydicom
    from pydicom.dataset import Dataset as DicomDataset

    PYDICOM_AVAILABLE = True
except ImportError:
    PYDICOM_AVAILABLE = False
    logger.warning("pydicom not installed. DICOM I/O will be unavailable.")

try:
    import highdicom

    HIGHDICOM_AVAILABLE = True
except ImportError:
    HIGHDICOM_AVAILABLE = False
    logger.warning("highdicom not installed. DICOM SR generation will be unavailable.")

try:
    import SimpleITK as sitk

    SIMPLEITK_AVAILABLE = True
except ImportError:
    SIMPLEITK_AVAILABLE = False
    logger.warning("SimpleITK not installed. NIfTI conversion will be unavailable.")


# ═══════════════════════════════════════════════════════════════════════
# Series Selection Rules
# ═══════════════════════════════════════════════════════════════════════


class SeriesSelectionRule:
    """Configurable DICOM series selection rule for MAP input filtering.

    Attributes:
        modality: Required DICOM modality (e.g. 'CT', 'CR', 'MR').
        body_part: Expected BodyPartExamined value (e.g. 'HEAD', 'CHEST').
        description_keywords: Keywords to match in SeriesDescription.
        min_instances: Minimum number of instances (slices) required.
        preferred_slice_thickness_mm: Preferred slice thickness for scoring.
    """

    def __init__(
        self,
        modality: str,
        body_part: str = "",
        description_keywords: Optional[List[str]] = None,
        min_instances: int = 1,
        preferred_slice_thickness_mm: Optional[float] = None,
    ) -> None:
        self.modality = modality.upper()
        self.body_part = body_part.upper()
        self.description_keywords = [
            kw.lower() for kw in (description_keywords or [])
        ]
        self.min_instances = min_instances
        self.preferred_slice_thickness_mm = preferred_slice_thickness_mm

    def matches(self, series_meta: Dict[str, Any]) -> Tuple[bool, float]:
        """Check if a DICOM series matches this rule.

        Returns:
            Tuple of (matches: bool, score: float).  Higher score means
            better match among qualifying series.
        """
        score = 0.0

        # Modality must match
        series_modality = str(series_meta.get("Modality", "")).upper()
        if series_modality != self.modality:
            return False, 0.0
        score += 1.0

        # Body part (optional but scored)
        if self.body_part:
            series_body = str(series_meta.get("BodyPartExamined", "")).upper()
            if series_body == self.body_part:
                score += 1.0

        # Description keyword matching
        series_desc = str(series_meta.get("SeriesDescription", "")).lower()
        for kw in self.description_keywords:
            if kw in series_desc:
                score += 0.5

        # Instance count filtering
        num_instances = int(series_meta.get("NumberOfInstances", 0))
        if num_instances < self.min_instances:
            return False, 0.0
        score += min(num_instances / 100.0, 1.0)  # Reward more slices

        # Slice thickness scoring (prefer thinner slices)
        if self.preferred_slice_thickness_mm is not None:
            thickness = float(series_meta.get("SliceThickness", 0.0))
            if 0 < thickness <= self.preferred_slice_thickness_mm:
                score += 1.0
            elif thickness > 0:
                score += 0.5

        return True, score


# ═══════════════════════════════════════════════════════════════════════
# Base MONAI Deploy Application
# ═══════════════════════════════════════════════════════════════════════


class BaseImagingMAP(Application):
    """Base MONAI Deploy Application for HCLS AI Factory imaging workflows.

    Subclasses must implement:
        - ``get_series_selection_rules()`` — return series filter rules
        - ``get_workflow()`` — return the instantiated BaseImagingWorkflow
        - ``build_dicom_sr()`` — create DICOM SR from WorkflowResult
    """

    APP_NAME: str = "base_imaging_map"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = "Base Imaging MAP"
    TARGET_LATENCY_SEC: float = 120.0

    def __init__(self, mock_mode: bool = False, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.mock_mode = mock_mode
        self._workflow: Optional[Any] = None
        logger.info(
            f"Initialized {self.APP_NAME} v{self.APP_VERSION} "
            f"(mock={mock_mode}, deploy_sdk={MONAI_DEPLOY_AVAILABLE})"
        )

    # ── Abstract interface ──────────────────────────────────────────────

    @abstractmethod
    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return DICOM series selection rules for this application."""
        ...

    @abstractmethod
    def get_workflow(self) -> Any:
        """Return the instantiated BaseImagingWorkflow for this MAP."""
        ...

    @abstractmethod
    def build_dicom_sr(
        self,
        workflow_result: Any,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM Structured Report from workflow result.

        Args:
            workflow_result: WorkflowResult from the imaging workflow.
            source_dicom: Optional source DICOM dataset for patient/study context.

        Returns:
            Serialized DICOM SR bytes, or None if highdicom is unavailable.
        """
        ...

    # ── DICOM input loading ─────────────────────────────────────────────

    def load_dicom_series(self, input_dir: str) -> List[Dict[str, Any]]:
        """Load all DICOM series from input directory.

        Returns a list of dicts with series-level metadata and file paths.
        """
        input_path = Path(input_dir)

        if not input_path.exists():
            logger.warning(f"Input directory does not exist: {input_dir}")
            return []

        if not PYDICOM_AVAILABLE:
            logger.warning("pydicom not available — returning empty series list")
            return []

        series_map: Dict[str, Dict[str, Any]] = {}

        for dcm_file in sorted(input_path.rglob("*.dcm")):
            try:
                ds = pydicom.dcmread(str(dcm_file), stop_before_pixels=True)
                series_uid = str(getattr(ds, "SeriesInstanceUID", "unknown"))

                if series_uid not in series_map:
                    series_map[series_uid] = {
                        "SeriesInstanceUID": series_uid,
                        "StudyInstanceUID": str(
                            getattr(ds, "StudyInstanceUID", "")
                        ),
                        "Modality": str(getattr(ds, "Modality", "")),
                        "BodyPartExamined": str(
                            getattr(ds, "BodyPartExamined", "")
                        ),
                        "SeriesDescription": str(
                            getattr(ds, "SeriesDescription", "")
                        ),
                        "SliceThickness": float(
                            getattr(ds, "SliceThickness", 0.0)
                        ),
                        "PatientID": str(getattr(ds, "PatientID", "")),
                        "PatientName": str(getattr(ds, "PatientName", "")),
                        "NumberOfInstances": 0,
                        "files": [],
                    }

                series_map[series_uid]["NumberOfInstances"] += 1
                series_map[series_uid]["files"].append(str(dcm_file))

            except Exception as e:
                logger.debug(f"Skipping non-DICOM file {dcm_file}: {e}")

        series_list = list(series_map.values())
        logger.info(
            f"Loaded {len(series_list)} DICOM series from {input_dir} "
            f"({sum(s['NumberOfInstances'] for s in series_list)} total instances)"
        )
        return series_list

    def select_best_series(
        self, series_list: List[Dict[str, Any]]
    ) -> Optional[Dict[str, Any]]:
        """Select the best matching DICOM series using configured rules.

        Evaluates each series against all selection rules and returns the
        series with the highest cumulative score.
        """
        rules = self.get_series_selection_rules()

        if not series_list:
            logger.warning("No DICOM series to select from")
            return None

        best_series: Optional[Dict[str, Any]] = None
        best_score = -1.0

        for series in series_list:
            total_score = 0.0
            any_match = False

            for rule in rules:
                matches, score = rule.matches(series)
                if matches:
                    any_match = True
                    total_score += score

            if any_match and total_score > best_score:
                best_score = total_score
                best_series = series

        if best_series is not None:
            logger.info(
                f"Selected series {best_series['SeriesInstanceUID'][:16]}... "
                f"({best_series['Modality']}, "
                f"{best_series['NumberOfInstances']} instances, "
                f"score={best_score:.2f})"
            )
        else:
            logger.warning(
                "No series matched selection rules. "
                f"Available modalities: "
                f"{[s.get('Modality', '?') for s in series_list]}"
            )

        return best_series

    # ── DICOM to NIfTI conversion ───────────────────────────────────────

    def dicom_to_nifti(
        self, dicom_files: List[str], output_dir: Optional[str] = None
    ) -> Optional[str]:
        """Convert a DICOM series to NIfTI for MONAI/workflow inference.

        Uses SimpleITK's DICOM reader to create a properly oriented
        3D volume from the series files.

        Args:
            dicom_files: List of DICOM file paths belonging to one series.
            output_dir: Directory for the NIfTI output. Uses tempdir if None.

        Returns:
            Path to the generated NIfTI file, or None on failure.
        """
        if not SIMPLEITK_AVAILABLE:
            logger.warning("SimpleITK not available — cannot convert DICOM to NIfTI")
            return None

        if not dicom_files:
            logger.warning("No DICOM files provided for NIfTI conversion")
            return None

        try:
            reader = sitk.ImageSeriesReader()
            reader.SetFileNames(dicom_files)
            reader.MetaDataDictionaryArrayUpdateOn()
            reader.LoadPrivateTagsOn()
            image = reader.Execute()

            if output_dir is None:
                output_dir = tempfile.mkdtemp(prefix="map_nifti_")

            output_path = Path(output_dir) / "volume.nii.gz"
            sitk.WriteImage(image, str(output_path))

            logger.info(
                f"Converted {len(dicom_files)} DICOM files to NIfTI: "
                f"{output_path} (size={image.GetSize()}, "
                f"spacing={image.GetSpacing()})"
            )
            return str(output_path)

        except Exception as e:
            logger.error(f"DICOM to NIfTI conversion failed: {e}")
            return None

    # ── Output writing ──────────────────────────────────────────────────

    def write_output(
        self,
        output_dir: str,
        workflow_result: Any,
        dicom_sr_bytes: Optional[bytes] = None,
    ) -> Dict[str, str]:
        """Write workflow outputs to the output directory.

        Generates:
            - result.json: Full WorkflowResult as JSON
            - dicom_sr.dcm: DICOM Structured Report (if available)

        Returns:
            Dict mapping output type to file path.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        outputs: Dict[str, str] = {}

        # JSON result
        json_path = output_path / "result.json"
        result_dict = (
            workflow_result.model_dump()
            if hasattr(workflow_result, "model_dump")
            else vars(workflow_result)
        )
        json_path.write_text(json.dumps(result_dict, indent=2, default=str))
        outputs["json"] = str(json_path)
        logger.info(f"Wrote JSON result: {json_path}")

        # DICOM SR
        if dicom_sr_bytes is not None:
            sr_path = output_path / "dicom_sr.dcm"
            sr_path.write_bytes(dicom_sr_bytes)
            outputs["dicom_sr"] = str(sr_path)
            logger.info(f"Wrote DICOM SR: {sr_path}")

        return outputs

    # ── Main execution pipeline ─────────────────────────────────────────

    def run_standalone(
        self, input_dir: str, output_dir: str
    ) -> Dict[str, Any]:
        """Run the MAP pipeline outside of MONAI Deploy runtime.

        This method provides the same DICOM-in / result-out pipeline
        that the MONAI Deploy operator graph would execute, but can be
        called directly for testing or integration without the SDK.

        Args:
            input_dir: Directory containing DICOM input files.
            output_dir: Directory for pipeline outputs.

        Returns:
            Dict with workflow_result and output file paths.
        """
        start = time.time()
        logger.info(
            f"Running {self.APP_NAME} standalone: "
            f"input={input_dir}, output={output_dir}"
        )

        # 1. Load DICOM series
        series_list = self.load_dicom_series(input_dir)

        # 2. Select best series
        selected = self.select_best_series(series_list)
        source_dicom = None

        # 3. Convert to NIfTI (for 3D workflows)
        nifti_path: Optional[str] = None
        if selected and selected.get("files"):
            nifti_path = self.dicom_to_nifti(selected["files"])
            # Load first DICOM for SR context
            if PYDICOM_AVAILABLE and selected["files"]:
                try:
                    source_dicom = pydicom.dcmread(selected["files"][0])
                except Exception:
                    pass

        # 4. Run workflow
        workflow = self.get_workflow()
        input_path = nifti_path or input_dir
        workflow_result = workflow.run(input_path=input_path)

        # 5. Generate DICOM SR
        dicom_sr_bytes = self.build_dicom_sr(workflow_result, source_dicom)

        # 6. Write outputs
        output_files = self.write_output(output_dir, workflow_result, dicom_sr_bytes)

        elapsed = time.time() - start
        logger.info(
            f"{self.APP_NAME} completed in {elapsed:.1f}s "
            f"(target={self.TARGET_LATENCY_SEC}s, "
            f"severity={workflow_result.severity.value})"
        )

        return {
            "workflow_result": workflow_result,
            "output_files": output_files,
            "elapsed_sec": elapsed,
        }

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph.

        When running under the MONAI Deploy runtime, this method wires
        together the SDK operators (DICOM loader -> series selector ->
        inference -> SR writer).  In standalone mode this is unused.
        """
        if not MONAI_DEPLOY_AVAILABLE:
            logger.info("MONAI Deploy SDK not available — skipping compose()")
            return

        logger.info(f"Composing MONAI Deploy operator graph for {self.APP_NAME}")
        # Operator graph is built by concrete subclasses that can
        # instantiate SDK operators and call self.add_flow().
