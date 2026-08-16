#!/usr/bin/env python3
"""Generate one infographic per subject — the 5th asset of the 85.

House style is fixed in STYLE so all 17 read as one family, matching the existing
services/google/out/hcls-eight-engines.png reference set.

Two rules learned the hard way and encoded here:
  * "no drop shadow, no card, no panel, no frame" — a generated figure with a drop shadow does not
    flatten onto a paper background and ships as a visible grey card edge.
  * the exact label text is supplied, because the model invents captions otherwise and a wrong
    clinical label in a published figure is worse than no figure.

Usage:
    .venv/bin/python scripts/gen_subject_infographics.py            # all 17 (costs money)
    .venv/bin/python scripts/gen_subject_infographics.py A3         # one
    .venv/bin/python scripts/gen_subject_infographics.py --dry      # print prompts, generate nothing
"""
from __future__ import annotations
import json, pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
NB = "/home/adam/projects/services/google/.venv/bin/nanobanana"
OUT = ROOT / "docs" / "assets" / "infographics" / "subjects"

STYLE = (
    "STYLE: clean minimal flat line-art scientific diagram, educational textbook style, on a solid "
    "warm cream background (hex F5F1E8) filling the whole canvas edge to edge. Thin dark navy "
    "(hex 1B2A4A) outlines, white or cream fills, colour accents only in sage green (hex 5E9B8C) "
    "and warm amber (hex E8A44C). Absolutely flat vector artwork: no shading, no gradients, no drop "
    "shadows, no card, no panel, no border or frame around the artwork. Clean sans-serif dark-navy "
    "labels centred beneath each element. No title text, no extra words beyond the labels given."
)

# scene + exact labels, per subject. Written by hand: a wrong label in a clinical figure is worse
# than no figure at all.
SCENES = {
 "E1": ("four stages left to right joined by thin arrows: a DNA double helix; a set of short "
        "fragment lines; a magnifying glass over a column of letters; a single highlighted letter "
        "in a sequence", ["Genome", "Reads", "Aligned", "Variant"]),
 "E2": ("three stages: a question mark in a speech bubble; a stack of document cards with small "
        "dots suggesting a vector index; a document with a highlighted line and a citation marker",
        ["Question", "Evidence", "Cited Answer"]),
 "E3": ("four stages: a ribbon protein with a visible cleft; a single small molecule; a fan of "
        "several different small molecules; a ranked list of three molecules with bars",
        ["Target", "Seed", "Candidates", "Ranked"]),
 # The first attempt let the model write its own report rows and it invented clinical grades
 # ("Severe Blockage", "Significant Narrowing") that are NOT CAD-RADS terminology. A wrong clinical
 # label in a published figure is worse than no figure, so the report card is now explicitly blank.
 "E4": ("three stages: a CT scanner gantry; a cross-sectional scan slice with one highlighted "
        "vessel; a clipboard showing a report card with THREE BLANK UNLABELLED ROWS and empty "
        "checkboxes, containing absolutely no words or text of any kind inside the card",
        ["Scan", "Finding", "Structured Report"]),
 "E5": ("four stages: a tumour cell cluster; a DNA strand with one marked position; a target "
        "symbol; three stacked therapy cards of decreasing size",
        ["Tumour", "Variant", "Target", "Therapy Tiers"]),
 "E6": ("three stages: a stylised heart; a line chart trending upward; a gauge dial with a marked "
        "band", ["Heart", "Signals", "Risk"]),
 "E7": ("three stages: a folded ribbon protein; the same protein with a highlighted surface pocket; "
        "a small molecule sitting in that pocket", ["Structure", "Pocket", "Fit"]),
 "E8": ("three stages: a tube of mixed cells; a scatter plot with distinct coloured clusters; four "
        "small labelled cell icons", ["Sample", "Clusters", "Cell Types"]),
 "A1": ("four stages joined in a cycle: a T-cell; a T-cell with a receptor added on its surface; "
        "several identical T-cells; a T-cell contacting a cancer cell",
        ["T-Cell", "Engineered", "Expanded", "Targeted"]),
 "A2": ("three stages: a blood tube; a marker molecule highlighted among others; two diverging "
        "arrows labelled beneath", ["Sample", "Marker", "Decision"]),
 "A3": ("three stages: two identical pill icons; two different enzyme shapes below them; two "
        "diverging outcome arrows of different lengths", ["Same Dose", "Enzyme Variant", "Different Outcome"]),
 "A4": ("three stages: an immune cell; the same cell with arrows turned toward a body outline; a "
        "calendar with one marked earlier date", ["Immune Cell", "Self-Attack", "Earlier Action"]),
 "A5": ("three stages: a brain outline; a numbered assessment checklist; a clock face with a "
        "highlighted narrow segment", ["Brain", "Assessment", "Time Window"]),
 "A6": ("three stages: a patient outline; a funnel with criteria lines; a folder marked with a "
        "check", ["Patient", "Criteria", "Matched Trial"]),
 "A7": ("four stages: a person outline with several small feature markers; a list of structured "
        "phenotype terms; a ranked list of three conditions; a single highlighted answer",
        ["Features", "Terms", "Differential", "Candidate"]),
 "A8": ("three stages: a scatter plot of coloured clusters; a magnifier over one cluster; a short "
        "written interpretation card", ["Clusters", "Focus", "Interpretation"]),
 "P1": ("four stages around a central infant outline: a DNA strand; a brain with small nodules; a "
        "heart with a small mass; a kidney outline",
        ["Genome", "Neurological", "Cardiac", "Renal"]),
}


def prompt_for(key, scene, labels):
    lab = ", ".join(f"'{x}'" for x in labels)
    return (f"A scientific process diagram showing {scene}. "
            f"Label the elements exactly and only: {lab}. {STYLE}")


def main():
    only = [a.upper() for a in sys.argv[1:] if not a.startswith("--")]
    dry = "--dry" in sys.argv
    OUT.mkdir(parents=True, exist_ok=True)
    keys = [k for k in SCENES if not only or k in only]
    print(f"{len(keys)} infographic(s) -> {OUT.relative_to(ROOT)}")
    ok = 0
    for k in keys:
        scene, labels = SCENES[k]
        p = prompt_for(k, scene, labels)
        dst = OUT / f"{k.lower()}.png"
        if dry:
            print(f"\n--- {k} -> {dst.name}\n{p}\n")
            continue
        r = subprocess.run([NB, p, "--aspect", "16:9", "--resolution", "2K", "-o", str(dst)],
                           capture_output=True, text=True, timeout=900)
        if dst.exists():
            ok += 1
            print(f"  OK   {k}  {dst.name}  ({dst.stat().st_size/1024:.0f} KB)")
        else:
            print(f"  FAIL {k}  {(r.stderr or r.stdout).strip()[:160]}")
    if not dry:
        print(f"\n{ok}/{len(keys)} generated")
    return 0 if dry or ok == len(keys) else 1


if __name__ == "__main__":
    raise SystemExit(main())
