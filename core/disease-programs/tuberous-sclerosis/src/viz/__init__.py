"""
TSC Anatomical Digital Twin — Surface (d) visualization layer (Phase V0).

Authors the engine's projections as OpenUSD scenes for volumetric, time-resolved viewing
in NVIDIA Omniverse. Authoring is a deterministic, CPU-side transform that runs on the
DGX Spark with no GPU and no third-party dependency (the ASCII `.usda` writer); RTX
rendering happens off-box (RunPod / an RTX workstation). The render's crispness equals the
engine's certainty: the uncertainty envelope radii ARE the forecast prediction intervals.
"""
from src.viz.export import author_population, author_scene, scene_for  # noqa: F401
from src.viz.scene_spec import SceneSpec  # noqa: F401
