"""
MDL materials for film-quality RTX rendering (FR-VZ-21, V1 polish). Emits a `/World/Materials`
scope of NVIDIA MDL materials (OmniGlass / OmniPBR — built into Omniverse) and binds them to
the static-colour elements: the uncertainty envelope and threshold render as true glass, the
mosaic variant cells and the recovery halo as emissive glow.

Honest fallbacks: the gprims keep their `displayColor`/`displayOpacity`, so a non-Omniverse
USD viewer still shows colour; and the LESION is deliberately NOT material-bound so its
time-sampled threshold-state colour animation survives. The MDL modules resolve by name in
Omniverse; elsewhere the render falls back to displayColor.
"""
from __future__ import annotations


def _glass(name: str, rgb, ior: float, roughness: float) -> str:
    r, g, b = rgb
    return (
        f'        def Material "{name}"\n        {{\n'
        f'            token outputs:mdl:surface.connect = </World/Materials/{name}/Shader.outputs:out>\n'
        f'            def Shader "Shader"\n            {{\n'
        '                uniform token info:implementationSource = "sourceAsset"\n'
        "                uniform asset info:mdl:sourceAsset = @OmniGlass.mdl@\n"
        '                uniform token info:mdl:sourceAsset:subIdentifier = "OmniGlass"\n'
        f"                color3f inputs:glass_color = ({r}, {g}, {b})\n"
        f"                float inputs:glass_ior = {ior}\n"
        "                bool inputs:thin_walled = 1\n"
        f"                float inputs:frosting_roughness = {roughness}\n"
        "                token outputs:out\n            }\n        }\n"
    )


def _emissive(name: str, rgb, intensity: float, opacity: float) -> str:
    r, g, b = rgb
    op = (
        "                bool inputs:enable_opacity = 1\n"
        f"                float inputs:opacity_constant = {opacity}\n"
    ) if opacity < 1.0 else ""
    return (
        f'        def Material "{name}"\n        {{\n'
        f'            token outputs:mdl:surface.connect = </World/Materials/{name}/Shader.outputs:out>\n'
        f'            def Shader "Shader"\n            {{\n'
        '                uniform token info:implementationSource = "sourceAsset"\n'
        "                uniform asset info:mdl:sourceAsset = @OmniPBR.mdl@\n"
        '                uniform token info:mdl:sourceAsset:subIdentifier = "OmniPBR"\n'
        f"                color3f inputs:diffuse_color_constant = ({r}, {g}, {b})\n"
        "                bool inputs:enable_emission = 1\n"
        f"                color3f inputs:emissive_color = ({r}, {g}, {b})\n"
        f"                float inputs:emissive_intensity = {intensity}\n"
        f"{op}"
        "                token outputs:out\n            }\n        }\n"
    )


# name -> definition. Bound by the writers; the lesion stays displayColor-animated (unbound).
_MATERIALS = (
    _glass("GlassEnvelope", (0.5, 0.69, 1.0), 1.05, 0.12),
    _glass("ThresholdGlass", (0.84, 0.2, 0.2), 1.1, 0.05),
    _emissive("VariantEmissive", (1.0, 0.78, 0.2), 9.0, 1.0),
    _emissive("HaloEmissive", (1.0, 0.78, 0.2), 4.0, 0.25),
)

MATERIAL_NAMES = ("GlassEnvelope", "ThresholdGlass", "VariantEmissive", "HaloEmissive")


def materials_scope() -> str:
    """The `/World/Materials` scope (deterministic, static)."""
    return '    def Scope "Materials"\n    {\n' + "".join(_MATERIALS) + "    }\n\n"


def bind(material: str, indent: str) -> str:
    """A `material:binding` line for a gprim (assumes the Materials scope under /World)."""
    return f"{indent}rel material:binding = </World/Materials/{material}>\n"
