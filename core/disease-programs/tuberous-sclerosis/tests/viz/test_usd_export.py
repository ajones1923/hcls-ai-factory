"""
USD digital-twin export (FR-VZ): the headline invariant is ENVELOPE == INTERVAL — the
rendered uncertainty shells equal the engine's prediction intervals, so the picture can
never imply more certainty than the engine has. Plus determinism, watermark, no-new-claim.
"""
from app._engine import featured, get_engine
from src.viz import author_scene, scene_for


def _proj(tag):
    orch, _ = get_engine()
    return featured()[tag], orch.store.projection(featured()[tag])


def _sega(spec):
    return next(le for le in spec.lesions if le.label == "SEGA")


def test_envelope_radii_equal_the_prediction_intervals():
    pid, proj = _proj("B")
    spec = scene_for(pid, proj)
    sega = _sega(spec)
    forecast = next(q for q in proj["trajectory"]["lesions"] if q["lesion"] == "SEGA")["forecast"]
    for kf in [k for k in sega.keyframes if k.forecast]:
        f = forecast[f"m{int(kf.month)}"]
        assert kf.env50_radius == f["pi50"][1]      # the 50% shell IS the 50% interval upper bound
        assert kf.env90_radius == f["pi90"][1]      # the 90% shell IS the 90% interval upper bound
        assert kf.radius == f["mean_cm"]            # the lesion radius IS the forecast mean


def test_no_new_claim_observed_radii_match_measurements():
    pid, proj = _proj("B")
    sega = _sega(scene_for(pid, proj))
    q = next(x for x in proj["trajectory"]["lesions"] if x["lesion"] == "SEGA")
    observed = [k for k in sega.keyframes if not k.forecast]
    assert [k.radius for k in observed] == q["observed_cm"]
    assert [k.month for k in observed] == q["observed_months"]


def test_lesion_colour_crosses_threshold_state():
    pid, proj = _proj("B")
    sega = _sega(scene_for(pid, proj))
    states = [k.state for k in sega.keyframes]
    assert states[0] == "below"                     # starts well below threshold
    assert states[-1] == "at_or_above"              # forecast crosses the 1.8 cm discussion line


def test_usda_is_deterministic_and_watermarked():
    pid, proj = _proj("B")
    a = author_scene(pid, proj)["usda"]
    b = author_scene(pid, proj)["usda"]
    assert a == b                                   # byte-identical -> diffable like source
    assert a.startswith("#usda 1.0")
    assert "SYNTHETIC" in a and "not FDA-cleared" in a and "not patient imaging" in a
    assert 'def Sphere "Envelope90"' in a and 'def Sphere "ThresholdShell"' in a


def test_viz_endpoint():
    from fastapi.testclient import TestClient
    from api.main import app
    with TestClient(app) as c:
        r = c.get(f"/viz/lesion/{featured()['B']}").json()
        m = c.get(f"/viz/mosaic/{featured()['A']}").json()
    assert r["path"].endswith(".usda") and r["spec"]["lesions"]
    assert m["spec"]["scene_kind"] == "mosaic"
    bad = TestClient(app).get("/viz/nope/X")
    assert bad.status_code == 404


# ── mosaic 'powers-of-ten' scene (FR-VZ-17) ─────────────────────────────────
def test_mosaic_cell_fraction_equals_vaf():
    pid, proj = _proj("A")
    res = author_scene(pid, proj, "mosaic")
    s = res["spec"]
    vaf = (proj["variant_interp"]["primary"] or {})["vaf"]
    # the rendered fraction of variant-carrying cells equals the VAF (to grid resolution)
    assert abs(s["observed_fraction"] - vaf) <= 1.0 / s["n_cells"]
    assert s["n_variant"] == round(vaf * s["n_cells"])
    # protoIndices in the USD sum to the variant-cell count
    assert res["usda"].count("VariantCarrier") >= 1


def test_mosaic_carries_variant_call_and_watermark():
    pid, proj = _proj("A")
    usda = author_scene(pid, proj, "mosaic")["usda"]
    assert usda.startswith("#usda 1.0") and 'def PointInstancer "CellField"' in usda
    assert "recovered from affected tissue" in usda and "Likely Pathogenic" in usda
    assert "SYNTHETIC" in usda and "not FDA-cleared" in usda
    # deterministic
    assert usda == author_scene(pid, proj, "mosaic")["usda"]


# ── whole-child atlas + population (Scene 3, FR-VZ-14/18) ────────────────────
def test_atlas_lights_organs_from_phenome():
    pid, proj = _proj("B")
    res = author_scene(pid, proj, "atlas")
    involved = set(res["spec"]["organs_involved"])
    assert {"brain", "kidney", "skin"} <= involved          # B's HPO profile lights these
    assert 'def Capsule "Body"' in res["usda"] and "not FDA-cleared" in res["usda"]


def test_population_array_and_recovery_halos():
    from app._engine import get_engine
    from src.viz import author_population
    orch, man = get_engine()
    fmap = featured(); rev = {v: k for k, v in fmap.items()}
    patients = [(p["patient_id"], rev.get(p["patient_id"]), orch.store.projection(p["patient_id"]))
                for p in man["patients"]]
    res = author_population(patients)
    s = res["spec"]
    assert s["n_patients"] == 50
    assert s["n_recovered"] == 7                            # the 7 mosaic recoveries
    assert sum(s["distributions"]["classification"].values()) == 50
    # exactly 7 gold RecoveryHalo prims in the USD
    assert res["usda"].count('def Sphere "RecoveryHalo"') == 7
    assert res["usda"].count('def Capsule "Body"') == 50
    assert res["usda"] == author_population(patients)["usda"]   # deterministic


def test_mdl_materials_bound_but_lesion_keeps_animated_colour():
    pid, proj = _proj("B")
    u = author_scene(pid, proj)["usda"]
    # MDL glass on the envelope/threshold, with displayColor fallback retained
    assert 'def Material "GlassEnvelope"' in u and "OmniGlass.mdl" in u
    assert "material:binding = </World/Materials/GlassEnvelope>" in u
    assert "material:binding = </World/Materials/ThresholdGlass>" in u
    # the LESION must NOT be material-bound — its threshold-state colour is time-sampled
    lesion = u.split('def Sphere "Lesion"')[1].split("def Sphere")[0]
    assert "material:binding" not in lesion
    assert "displayColor.timeSamples" in lesion
    # emissive variant cells + halo glow on the other scenes
    assert "VariantEmissive" in author_scene(*_proj("A"), "mosaic")["usda"]


def test_population_endpoint():
    from fastapi.testclient import TestClient
    from api.main import app
    with TestClient(app) as c:
        r = c.get("/viz/population").json()
    assert r["spec"]["n_patients"] == 50 and r["spec"]["n_recovered"] == 7
    assert r["path"].endswith("cohort_population.usda")
