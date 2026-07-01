"""Reference imaging workflows for the Imaging Intelligence Agent."""

from src.workflows.base import BaseImagingWorkflow
from src.workflows.ct_head_hemorrhage import CTHeadHemorrhageWorkflow
from src.workflows.ct_chest_lung_nodule import CTChestLungNoduleWorkflow
from src.workflows.ct_coronary_angiography import CTCoronaryAngiographyWorkflow
from src.workflows.cxr_rapid_findings import CXRRapidFindingsWorkflow
from src.workflows.mri_brain_ms_lesion import MRIBrainMSLesionWorkflow
from src.workflows.mri_prostate_pirads import MRIProstateWorkflow
from src.workflows.breast_birads import BreastBIRADSWorkflow
from src.workflows.thyroid_tirads import ThyroidTIRADSWorkflow
from src.workflows.liver_lirads import LiverLIRADSWorkflow

WORKFLOW_REGISTRY = {
    "ct_head_hemorrhage": CTHeadHemorrhageWorkflow,
    "ct_chest_lung_nodule": CTChestLungNoduleWorkflow,
    "ct_coronary_angiography": CTCoronaryAngiographyWorkflow,
    "cxr_rapid_findings": CXRRapidFindingsWorkflow,
    "mri_brain_ms_lesion": MRIBrainMSLesionWorkflow,
    "mri_prostate_pirads": MRIProstateWorkflow,
    "breast_birads": BreastBIRADSWorkflow,
    "thyroid_tirads": ThyroidTIRADSWorkflow,
    "liver_lirads": LiverLIRADSWorkflow,
}

__all__ = [
    "BaseImagingWorkflow",
    "CTHeadHemorrhageWorkflow",
    "CTChestLungNoduleWorkflow",
    "CTCoronaryAngiographyWorkflow",
    "CXRRapidFindingsWorkflow",
    "MRIBrainMSLesionWorkflow",
    "MRIProstateWorkflow",
    "BreastBIRADSWorkflow",
    "ThyroidTIRADSWorkflow",
    "LiverLIRADSWorkflow",
    "WORKFLOW_REGISTRY",
]
