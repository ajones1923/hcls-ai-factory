#!/usr/bin/env python
"""
Author OpenUSD digital-twin scenes from the engine's projections (Spark-side, CPU only).

    venv/bin/python scripts/export_usd.py            # featured A/B/C
    venv/bin/python scripts/export_usd.py --all      # whole cohort

Render the resulting .usda files with RTX in NVIDIA Omniverse / USD Composer on a RunPod
RTX instance (or any RTX workstation). SYNTHETIC, atlas-anchored; not patient imaging.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app._engine import featured, get_engine          # noqa: E402
from src.viz import author_population, author_scene     # noqa: E402


def main() -> None:
    orch, manifest = get_engine()
    fmap = featured(); rev = {v: k for k, v in fmap.items()}
    pids = [p["patient_id"] for p in manifest["patients"]] if "--all" in sys.argv else list(fmap.values())
    print(f"Authoring scenes for {len(pids)} patient(s)…\n")
    n = 0
    for pid in pids:
        proj = orch.store.projection(pid)
        les = author_scene(pid, proj, "lesion_trajectory")
        if les["spec"]["lesions"]:
            n += 1
            names = ", ".join(f"{l['label']} (grade {l['crossing_grade']})" for l in les["spec"]["lesions"])
            print(f"  {pid} [lesion]: {names}\n      -> {les['path']}")
        atl = author_scene(pid, proj, "atlas"); n += 1
        print(f"  {pid} [atlas]: organs {atl['spec']['organs_involved']}\n      -> {atl['path']}")
        if (proj.get("variant_interp") or {}).get("primary", {}).get("recovered"):
            mos = author_scene(pid, proj, "mosaic"); s = mos["spec"]; n += 1
            print(f"  {pid} [mosaic]: {s['gene']} {s['variant']} VAF {s['vaf']} "
                  f"({s['n_variant']}/{s['n_cells']} cells)\n      -> {mos['path']}")
    # the cohort population scene (always, over all 50)
    patients = [(p["patient_id"], rev.get(p["patient_id"]), orch.store.projection(p["patient_id"]))
                for p in manifest["patients"]]
    pop = author_population(patients); ps = pop["spec"]; n += 1
    print(f"\n  POPULATION: {ps['n_patients']} patients, {ps['n_recovered']} mosaic recoveries (gold halos)"
          f"\n      -> {pop['path']}")
    print(f"\n{n} scene(s) written to data/usd/. Open in Omniverse (RunPod RTX) to render.")


if __name__ == "__main__":
    main()
