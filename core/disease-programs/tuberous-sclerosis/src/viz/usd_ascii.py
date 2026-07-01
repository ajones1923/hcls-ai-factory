"""
Dependency-free OpenUSD (.usda) writer (FR-VZ-10b). Emits a valid ASCII USD stage from a
SceneSpec with no `pxr`/usd-core dependency, so the DGX Spark authors render-ready scenes
on CPU alone. Output is deterministic (stable ordering, fixed float formatting, no
timestamps) so scenes diff like source. Open the result in NVIDIA Omniverse / USD Composer
(RunPod RTX) for the path-traced view.

Encoding: each lesion is a Sphere whose `radius` and `displayColor` are USD time-samples
keyed by month (= time code); two translucent Sphere shells carry the 50%/90% prediction-
interval radii (the uncertainty envelope); a translucent Sphere marks the clinical
threshold; provenance + the SYNTHETIC watermark are stage/prim customData.
"""
from __future__ import annotations

from src.viz.materials import bind, materials_scope
from src.viz.palette import ANATOMY_RGB, ENVELOPE_RGB, RECOVERY_RGB, THRESHOLD_RGB
from src.viz.scene_spec import (
    AtlasScene, BodyFigure, LesionScene, MosaicScene, PopulationScene, SceneSpec,
)

_CELL_NORMAL_RGB = (0.16, 0.20, 0.30)    # quiet wild-type cell
_CELL_VARIANT_RGB = (1.0, 0.78, 0.20)    # glowing variant-carrying cell (the recovered mosaic)


def _n(x: float) -> str:
    """Deterministic number formatting."""
    if x == int(x):
        return str(int(x))
    return f"{x:.4f}".rstrip("0").rstrip(".")


def _c(rgb) -> str:
    return f"({_n(rgb[0])}, {_n(rgb[1])}, {_n(rgb[2])})"


def _radius_samples(kfs, attr: str) -> str:
    rows = ",\n".join(f"            {_n(k.time_code)}: {_n(getattr(k, attr))}" for k in kfs)
    return "{\n" + rows + "\n        }"


def _color_samples(kfs) -> str:
    rows = ",\n".join(f"            {_n(k.time_code)}: [{_c(k.color)}]" for k in kfs)
    return "{\n" + rows + "\n        }"


def _sphere(name, radius_attr, color_attr, opacity, indent="        ", material=None) -> str:
    mat = bind(material, indent + "    ") if material else ""
    return (
        f'{indent}def Sphere "{name}"\n{indent}{{\n'
        f'{indent}    double radius{radius_attr}\n'
        f'{indent}    color3f[] primvars:displayColor{color_attr}\n'
        f'{indent}    float[] primvars:displayOpacity = [{_n(opacity)}]\n'
        f'{mat}'
        f'{indent}}}\n'
    )


def _lesion_block(le: LesionScene) -> str:
    x, y, z = le.anchor
    out = [
        f'    def Xform "{_id(le.label)}"\n    {{\n',
        f'        double3 xformOp:translate = ({_n(x)}, {_n(y)}, {_n(z)})\n',
        f'        uniform token[] xformOpOrder = ["xformOp:translate"]\n',
        f'        custom string tsc:location = "{le.location}"\n',
        f'        custom string tsc:crossing_grade = "{le.crossing_grade}"\n',
        f'        custom double tsc:threshold_{le.unit.replace("/", "_")} = {_n(le.threshold)}\n',
    ]
    # stylized context organ (very translucent)
    out.append(_sphere("Context", f" = {_n(le.context_radius)}",
                       f" = [{_c(ANATOMY_RGB)}]", 0.05))
    # the lesion: animated radius + threshold-state colour
    out.append(
        f'        def Sphere "Lesion"\n        {{\n'
        f'            double radius.timeSamples = {_radius_samples(le.keyframes, "radius")}\n'
        f'            color3f[] primvars:displayColor.timeSamples = {_color_samples(le.keyframes)}\n'
        f'            float[] primvars:displayOpacity = [1]\n'
        f'            custom string tsc:provenance_section = "trajectory"\n'
        f'        }}\n'
    )
    # uncertainty envelope: 50% (inner, denser) and 90% (outer, softer) — radii ARE the PIs.
    # Glass MDL in Omniverse; displayColor fallback elsewhere.
    out.append(_sphere("Envelope50", f".timeSamples = {_radius_samples(le.keyframes, 'env50_radius')}",
                       f" = [{_c(ENVELOPE_RGB)}]", 0.22, material="GlassEnvelope"))
    out.append(_sphere("Envelope90", f".timeSamples = {_radius_samples(le.keyframes, 'env90_radius')}",
                       f" = [{_c(ENVELOPE_RGB)}]", 0.10, material="GlassEnvelope"))
    # the intervention-threshold membrane (red glass)
    out.append(_sphere("ThresholdShell", f" = {_n(le.threshold_radius)}",
                       f" = [{_c(THRESHOLD_RGB)}]", 0.12, material="ThresholdGlass"))
    out.append("    }\n")
    return "".join(out)


def _id(label: str) -> str:
    s = "".join(ch if ch.isalnum() else "_" for ch in label)
    return s if s[:1].isalpha() or s[:1] == "_" else f"L_{s}"


def write_usda(spec: SceneSpec) -> str:
    header = (
        "#usda 1.0\n(\n"
        '    defaultPrim = "World"\n'
        '    upAxis = "Y"\n'
        "    metersPerUnit = 0.01\n"
        f"    startTimeCode = {_n(spec.start_tc)}\n"
        f"    endTimeCode = {_n(spec.end_tc)}\n"
        f"    timeCodesPerSecond = {_n(spec.time_codes_per_second)}\n"
        f'    doc = "TSC Anatomical Digital Twin — {spec.scene_kind} — patient {spec.patient_id}. '
        'SYNTHETIC demonstration data; atlas-anchored stylized anatomy (not patient imaging); '
        'decision support, clinician review required; not FDA-cleared."\n'
        "    customData = {\n"
        f'        string watermark = "{spec.watermark}"\n'
        f'        string patient_id = "{spec.patient_id}"\n'
        f'        string scene_kind = "{spec.scene_kind}"\n'
        f"        double cm_to_unit_scale = {_n(spec.scale)}\n"
        '        string anatomy = "atlas-anchored stylized (Phase V0); segmented patient anatomy is Phase V3"\n'
        "    }\n)\n\n"
        'def Xform "World"\n{\n'
        '    def Scope "Watermark" (\n'
        f'        customData = {{ string text = "{spec.watermark} — decision support, not FDA-cleared" }}\n'
        "    ) {}\n\n"
        + materials_scope()
    )
    body = "\n".join(_lesion_block(le) for le in spec.lesions)
    prov_rows = "; ".join(f"{p.get('section')}:{p.get('status','?')}" for p in spec.provenance) or "n/a"
    footer = (
        f'\n    def Scope "Provenance" (\n'
        f'        customData = {{ string records = "{prov_rows}"; '
        'string note = "every rendered element traces to a validated projection field" }\n'
        "    ) {}\n}\n"
    )
    return header + body + footer


def write_mosaic_usda(spec: MosaicScene) -> str:
    """Author the mosaic 'powers-of-ten' cell field as a USD PointInstancer: exactly the
    recovered VAF fraction of cells carries the (glowing) variant prototype."""
    positions = ", ".join(f"({_n(c.pos[0])}, {_n(c.pos[1])}, {_n(c.pos[2])})" for c in spec.cells)
    proto = ", ".join("1" if c.variant else "0" for c in spec.cells)
    crit = "; ".join(f"{c['code']}({c['bucket']})" for c in spec.criteria) or "—"
    header = (
        "#usda 1.0\n(\n"
        '    defaultPrim = "World"\n    upAxis = "Y"\n    metersPerUnit = 1\n'
        f'    doc = "TSC Anatomical Digital Twin — mosaic powers-of-ten — patient {spec.patient_id}. '
        f'{spec.gene} {spec.variant} at VAF {_n(spec.vaf)} ({spec.n_variant}/{spec.n_cells} cells carry it) '
        '-> recovered from affected tissue; standard blood testing reports negative. '
        'SYNTHETIC; decision support; not FDA-cleared."\n'
        "    customData = {\n"
        f'        string watermark = "{spec.watermark}"\n'
        f'        string patient_id = "{spec.patient_id}"\n'
        f'        string scene_kind = "mosaic"\n'
        f'        string gene = "{spec.gene}"\n'
        f'        string variant = "{spec.variant}"\n'
        f"        double vaf = {_n(spec.vaf)}\n"
        f'        string classification = "{spec.classification} ({spec.acmg_rule})"\n'
        f'        string acmg_criteria = "{crit}"\n'
        f"        bool recovered_from_tissue = {'true' if spec.recovered else 'false'}\n"
        "    }\n)\n\n"
        'def Xform "World"\n{\n'
        '    def Scope "Watermark" (\n'
        f'        customData = {{ string text = "{spec.watermark} — decision support, not FDA-cleared" }}\n'
        "    ) {}\n\n"
        + materials_scope() +
        f'    def PointInstancer "CellField"\n    {{\n'
        f"        point3f[] positions = [{positions}]\n"
        f"        int[] protoIndices = [{proto}]\n"
        "        rel prototypes = [\n"
        "            </World/CellField/WildType>,\n"
        "            </World/CellField/VariantCarrier>,\n"
        "        ]\n"
        '        def Sphere "WildType"\n        {\n'
        "            double radius = 0.38\n"
        f"            color3f[] primvars:displayColor = [{_c(_CELL_NORMAL_RGB)}]\n"
        "            float[] primvars:displayOpacity = [0.6]\n        }\n"
        '        def Sphere "VariantCarrier"\n        {\n'
        "            double radius = 0.46\n"
        f"            color3f[] primvars:displayColor = [{_c(_CELL_VARIANT_RGB)}]\n"
        "            float[] primvars:displayOpacity = [1]\n"
        + bind("VariantEmissive", "            ") +
        f'            custom string tsc:variant = "{spec.gene} {spec.variant}"\n'
        "        }\n"
        "    }\n"
    )
    annotation = (
        f'\n    def Scope "VariantCall" (\n'
        f'        customData = {{ string gene = "{spec.gene}"; string variant = "{spec.variant}"; '
        f'string classification = "{spec.classification}"; string criteria = "{crit}"; '
        f'double vaf = {_n(spec.vaf)}; string recovery = '
        f'"{spec.n_variant} of {spec.n_cells} cells carry the variant; recovered from affected tissue" }}\n'
        "    ) {}\n}\n"
    )
    return header + annotation


# ── whole-child organ atlas + population (Scene 3) ──────────────────────────
def _figure_block(fig: BodyFigure, indent: str, body_h: float, body_r: float,
                  organ_r: float, name: str) -> str:
    x, y, z = fig.pos
    out = [
        f'{indent}def Xform "{name}"\n{indent}{{\n',
        f'{indent}    double3 xformOp:translate = ({_n(x)}, {_n(y)}, {_n(z)})\n',
        f'{indent}    uniform token[] xformOpOrder = ["xformOp:translate"]\n',
        f'{indent}    custom string tsc:patient_id = "{fig.patient_id}"\n',
        f'{indent}    custom string tsc:classification = "{fig.classification}"\n',
        f'{indent}    custom bool tsc:recovered_mosaic = {"true" if fig.recovered else "false"}\n',
        f'{indent}    custom int tsc:organ_burden = {fig.burden}\n',
    ]
    if fig.featured:
        out.append(f'{indent}    custom string tsc:featured = "{fig.featured}"\n')
    if fig.overdue:
        out.append(f'{indent}    custom string tsc:overdue_surveillance = "{", ".join(fig.overdue)}"\n')
    # body — capsule coloured by ACMG classification
    out.append(
        f'{indent}    def Capsule "Body"\n{indent}    {{\n'
        f'{indent}        uniform token axis = "Y"\n'
        f'{indent}        double height = {_n(body_h)}\n'
        f'{indent}        double radius = {_n(body_r)}\n'
        f'{indent}        color3f[] primvars:displayColor = [{_c(fig.body_rgb)}]\n'
        f'{indent}        float[] primvars:displayOpacity = [0.45]\n{indent}    }}\n'
    )
    # organs — lit if involved, faint otherwise
    for o in fig.organs:
        op = 1.0 if o.involved else 0.07
        rr = organ_r if o.involved else organ_r * 0.7
        ox, oy, oz = o.pos
        out.append(
            f'{indent}    def Sphere "{_id(o.organ)}"\n{indent}    {{\n'
            f'{indent}        double radius = {_n(round(rr, 3))}\n'
            f'{indent}        double3 xformOp:translate = ({_n(ox)}, {_n(oy)}, {_n(oz)})\n'
            f'{indent}        uniform token[] xformOpOrder = ["xformOp:translate"]\n'
            f'{indent}        color3f[] primvars:displayColor = [{_c(o.rgb)}]\n'
            f'{indent}        float[] primvars:displayOpacity = [{_n(op)}]\n{indent}    }}\n'
        )
    # gold halo for a recovered mosaic
    if fig.recovered:
        out.append(
            f'{indent}    def Sphere "RecoveryHalo"\n{indent}    {{\n'
            f'{indent}        double radius = {_n(round(body_h * 0.75, 3))}\n'
            f'{indent}        color3f[] primvars:displayColor = [{_c(RECOVERY_RGB)}]\n'
            f'{indent}        float[] primvars:displayOpacity = [0.12]\n'
            + bind("HaloEmissive", indent + "        ") +
            f'{indent}    }}\n'
        )
    out.append(f'{indent}}}\n')
    return "".join(out)


def write_atlas_usda(spec: AtlasScene) -> str:
    fig = spec.figure
    involved = ", ".join(o.organ for o in fig.organs if o.involved) or "none"
    header = (
        "#usda 1.0\n(\n"
        '    defaultPrim = "World"\n    upAxis = "Y"\n    metersPerUnit = 0.01\n'
        f'    doc = "TSC Anatomical Digital Twin — whole-child atlas — patient {spec.patient_id}. '
        f'Organs involved: {involved}. SYNTHETIC; atlas-anchored (not patient imaging); '
        'decision support; not FDA-cleared."\n'
        "    customData = {\n"
        f'        string watermark = "{spec.watermark}"\n'
        f'        string patient_id = "{spec.patient_id}"\n'
        f'        string scene_kind = "atlas"\n'
        f'        string classification = "{fig.classification}"\n'
        f'        string organs_involved = "{involved}"\n'
        f'        string overdue_surveillance = "{", ".join(fig.overdue) or "none"}"\n'
        "    }\n)\n\n"
        'def Xform "World"\n{\n'
        '    def Scope "Watermark" (\n'
        f'        customData = {{ string text = "{spec.watermark} — decision support, not FDA-cleared" }}\n'
        "    ) {}\n\n"
        + materials_scope()
    )
    body = _figure_block(fig, "    ", body_h=6.0, body_r=2.0, organ_r=0.7, name="Child")
    return header + body + "}\n"


def write_population_usda(spec: PopulationScene) -> str:
    cls = "; ".join(f"{k}:{v}" for k, v in spec.distributions.get("classification", {}).items())
    header = (
        "#usda 1.0\n(\n"
        '    defaultPrim = "World"\n    upAxis = "Y"\n    metersPerUnit = 0.01\n'
        f'    doc = "TSC Anatomical Digital Twin — population ({spec.n_patients} patients) — '
        f'{spec.n_recovered} mosaic variants recovered (gold halos). Body colour = ACMG class; '
        'lit organs = phenome involvement. SYNTHETIC; decision support; not FDA-cleared."\n'
        "    customData = {\n"
        f'        string watermark = "{spec.watermark}"\n'
        f'        string scene_kind = "population"\n'
        f"        int n_patients = {spec.n_patients}\n"
        f"        int n_recovered_mosaics = {spec.n_recovered}\n"
        f'        string classification_distribution = "{cls}"\n'
        "    }\n)\n\n"
        'def Xform "World"\n{\n'
        '    def Scope "Legend" (\n'
        f'        customData = {{ string note = "gold halo = recovered mosaic; body colour = ACMG class"; '
        f'string watermark = "{spec.watermark} — not FDA-cleared" }}\n'
        "    ) {}\n\n"
        + materials_scope() +
        '    def Xform "Cohort"\n    {\n'
    )
    figures = "".join(
        _figure_block(f, "        ", body_h=4.0, body_r=1.2, organ_r=0.45, name=f"P_{_id(f.patient_id)}")
        for f in spec.figures
    )
    return header + figures + "    }\n}\n"
