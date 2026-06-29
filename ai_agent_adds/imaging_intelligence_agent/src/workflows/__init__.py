"""Reference imaging workflows for the Imaging Intelligence Agent."""

from src.workflows.base import BaseImagingWorkflow
from src.workflows.ct_head_hemorrhage import CTHeadHemorrhageWorkflow
from src.workflows.ct_chest_lung_nodule import CTChestLungNoduleWorkflow
from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow
from src.workflows.mri_brain_ms_lesion import MRIBrainMSLesionWorkflow

WORKFLOW_REGISTRY = {
    "ct_head_hemorrhage": CTHeadHemorrhageWorkflow,
    "ct_chest_lung_nodule": CTChestLungNoduleWorkflow,
    "cxr_rapid_findings": CXRRapidFindingsWorkflow,
    "mri_brain_ms_lesion": MRIBrainMSLesionWorkflow,
}

__all__ = [
    "BaseImagingWorkflow",
    "CTHeadHemorrhageWorkflow",
    "CTChestLungNoduleWorkflow",
    "CXRRapidFindingsWorkflow",
    "MRIBrainMSLesionWorkflow",
    "WORKFLOW_REGISTRY",
]
