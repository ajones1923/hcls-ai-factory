#!/usr/bin/env python3
"""Pre-compute the coronary geometry analysis and cache it for the demo.

The GUI walkthrough needs the numbers present and correct, not computed live. Measuring once and
caching keeps the demo instant while keeping every figure REAL: the stenosis values below are
computed from the segmented vessel surface by src/coronary_geometry, not written by hand.

    python3 scripts/precompute_coronary_analysis.py
"""
import json, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from src.coronary_geometry import analyse, rater_agreement
from src.models import CADRADS
from src.workflows.ct_coronary_angiography import cadrads_report

BASE = pathlib.Path(__file__).resolve().parents[1]
# V1P2 is the demo case because its MEASURED grade is CAD-RADS 4A, matching the narrative the
# film and car-01.png already established. (V3P2 sits closer to the film's 72% figure at 68.5 %,
# but that grades CAD-RADS 3 and would break the story.) The 4A grade is the clinically meaningful
# claim and it is measured, not asserted; the exact percentage differs from the film's 72 %.
CASE = "V1P2"
# Representative inputs to the CAD-RADS modifiers, matching demo_cases.json DEMO-003.
AGATSTON_REPR = 385
HIGH_RISK_PLAQUE_REPR = True
MESH = BASE / "data" / "cardiac_ct" / "CoronariesNC6" / CASE
OUT = BASE / "data" / "demo" / "coronary" / "coronary_analysis.json"

def main():
    res = analyse(MESH / "manualGT_rater1.stl")
    agree = rater_agreement(MESH / "manualGT_rater1.stl", MESH / "manualGT_rater2.stl")
    ml = analyse(MESH / "ML.stl")
    payload = res.to_dict()
    payload["case"] = CASE
    payload["rater_agreement"] = agree
    payload["machine_segmentation"] = {
        "max_diameter_stenosis_pct": ml.max_diameter_stenosis_pct,
        "cad_rads": ml.cad_rads,
        "vs_manual_pct_points": round(abs(ml.max_diameter_stenosis_pct
                                          - res.max_diameter_stenosis_pct), 1),
    }
    # The full CAD-RADS 2.0 report string, built by the SAME function the workflow uses, so the
    # rendered panels and the API can never disagree about the modifiers. The plaque-burden and
    # high-risk-plaque inputs are representative for this case (Agatston needs Hounsfield units a
    # surface mesh does not carry), which is why the panel keeps saying "(repr.)".
    payload["cadrads_report"] = cadrads_report(
        CADRADS(payload["cad_rads"]), AGATSTON_REPR, HIGH_RISK_PLAQUE_REPR)
    payload["agatston_representative"] = AGATSTON_REPR
    payload["precomputed"] = True
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"  wrote {OUT.relative_to(BASE)}")
    print(f"  max stenosis {payload['max_diameter_stenosis_pct']}%  {payload['cadrads_report']}")
    print(f"  rater agreement: {agree}")
    print(f"  ML vs manual: {payload['machine_segmentation']['vs_manual_pct_points']} pct points")

    # Re-render the panels from the analysis we just wrote, WITH the reference heart.
    #
    # The renders draw their numbers out of coronary_analysis.json, so regenerating the analysis
    # without regenerating the images leaves the panels quoting a stale stenosis. Doing both here
    # means the JSON and the pictures can never disagree -- and the heart is included by default,
    # rather than depending on someone remembering a flag.
    import subprocess

    def _run(script, *args, label=""):
        cmd = [sys.executable, str(BASE / "scripts" / script), *args]
        print(f"  {label or script}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        for line in (r.stdout or "").splitlines():
            print(f"   {line}")
        if r.returncode != 0:
            print(f"  !! {script} failed ({r.returncode})")
            print((r.stderr or "")[-600:])
        return r.returncode

    _run("render_coronary_mesh.py", "--case", CASE, "--mesh", "manualGT_rater1.stl",
         "--spin", "--frames", "60", label="rendering 3D panels")
    # The cardiac story illustration quotes the same stenosis and marks the same segment, so it is
    # regenerated here too -- otherwise a re-run leaves the fourth panel disagreeing with the
    # other three about both the number and the location of the lesion.
    _run("restyle_cardiac_story.py", label="restyling the cardiac story illustration")

    # Panel D is generated art, so its numbers are pixels and cannot follow the analysis. Compare
    # them here: a silent divergence would put an illustration claiming one stenosis next to three
    # panels measuring another, which is the exact failure this whole pass existed to remove.
    side = OUT.parent / "cardiac_pathway_dark.json"
    if side.exists():
        baked = json.loads(side.read_text())["baked"]
        # cadrads_report carries the P/HRP modifiers. Checking only the bare grade let the panel
        # sit on "CAD-RADS 4A" while every other surface had moved to "4A/P3/HRP" — the grade
        # matched, so the guard stayed silent about a caption that was a modifier behind.
        live = {"max_diameter_stenosis_pct": round(float(payload["max_diameter_stenosis_pct"]), 1),
                "cad_rads": payload["cad_rads"],
                "cadrads_report": payload["cadrads_report"],
                "calcium_score_agatston": baked["calcium_score_agatston"]}
        drift = {k: (baked[k], live[k]) for k in baked if k in live and baked[k] != live[k]}
        if drift:
            print("  !! cardiac_pathway_dark.png is STALE — it has these values painted into it:")
            for k, (was, now) in drift.items():
                print(f"       {k}: image says {was}, analysis says {now}")
            print("     regenerate it with services/google + cardiac_pathway_dark.prompt.txt, "
                  "then update cardiac_pathway_dark.json")
        else:
            print(f"  panel D illustration agrees with the analysis ({baked['max_diameter_stenosis_pct']}%"
                  f" · {baked.get('cadrads_report', 'CAD-RADS ' + baked['cad_rads'])})")

if __name__ == "__main__":
    main()
