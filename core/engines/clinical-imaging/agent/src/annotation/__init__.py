"""MONAI Label interactive annotation integration.

Provides session-tracked annotation management connecting MONAI Label
with OHIF Viewer, VISTA-3D segmentation, and NVIDIA FLARE federated
learning for the clinical data flywheel.
"""

from .monai_label_config import AnnotationConfig, AnnotationSession, MONAILabelManager

__all__ = [
    "AnnotationConfig",
    "AnnotationSession",
    "MONAILabelManager",
]
