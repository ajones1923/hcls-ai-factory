#!/usr/bin/env python3
"""Generate synthetic clinical PDF documents for Precision Autoimmune Agent demo.

Creates ~120 realistic clinical PDFs across 5 demo patients, each representing
a multi-year autoimmune diagnostic odyssey. Documents include progress notes,
lab reports, imaging reports, pathology reports, genetic testing reports,
referral letters, medication reconciliations, and dismissal documentation.

Usage:
    python generate_demo_patients.py

Output: ../demo_data/{patient_name}/*.pdf

Author: Adam Jones
Date: March 2026
"""

import os
import sys
import time

# Ensure script directory is on path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Core patient generators
# New disease type patient generators
from patient_emma import generate as gen_emma
from patient_james import generate as gen_james
from patient_karen import generate as gen_karen
from patient_maya import generate as gen_maya
from patient_michael import generate as gen_michael
from patient_sarah import generate as gen_sarah
from patients_345 import generate_david, generate_linda, generate_rachel

# Additional specialist reports and subtle labs
from patients_additional import (
    generate_david_additional,
    generate_linda_additional,
    generate_maya_additional,
    generate_rachel_genetics,
    generate_sarah_additional,
)

# Enhanced timeline documents (gap-filling, post-dx monitoring)
from patients_david_enhanced import generate_david_enhanced

# Dismissal documents (symptoms attributed to wrong diagnoses)
from patients_dismissals import (
    generate_david_dismissals,
    generate_linda_dismissals,
    generate_maya_dismissals,
    generate_sarah_dismissals,
)
from patients_linda_enhanced import generate_linda_enhanced

# Medication reconciliation (polypharmacy revealing diagnostic signal)
from patients_med_lists import (
    generate_david_med_list,
    generate_linda_med_list,
    generate_maya_med_list,
    generate_rachel_med_list,
    generate_sarah_med_list,
)
from patients_rachel_enhanced import generate_rachel_enhanced

# Cross-specialist referral letters
from patients_referrals import (
    generate_david_referrals,
    generate_linda_referrals,
    generate_maya_referrals,
    generate_rachel_referrals,
    generate_sarah_referrals,
)

DEMO_DATA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "demo_data"
)


def main():
    print("=" * 70)
    print("Precision Autoimmune Agent -- Demo Patient Generator")
    print("Generating synthetic clinical PDFs for diagnostic odyssey demo")
    print("=" * 70)
    print()

    start = time.time()

    # ── Patient 1: Sarah Mitchell -- 31F, 3.5-year lupus odyssey ──
    print("[1/9] Sarah Mitchell (SLE / Lupus Nephritis, 3.5 years)...")
    sarah_dir = os.path.join(DEMO_DATA_DIR, "sarah_mitchell")
    gen_sarah(sarah_dir)
    generate_sarah_dismissals(sarah_dir)
    generate_sarah_med_list(sarah_dir)
    generate_sarah_referrals(sarah_dir)
    generate_sarah_additional(sarah_dir)

    # ── Patient 2: Maya Rodriguez -- 25F, 4-year POTS/hEDS/MCAS odyssey ──
    print("[2/9] Maya Rodriguez (POTS / hEDS / MCAS, 4 years)...")
    maya_dir = os.path.join(DEMO_DATA_DIR, "maya_rodriguez")
    gen_maya(maya_dir)
    generate_maya_dismissals(maya_dir)
    generate_maya_med_list(maya_dir)
    generate_maya_referrals(maya_dir)
    generate_maya_additional(maya_dir)

    # ── Patient 3: David Park -- 38M, 6-year ankylosing spondylitis odyssey ──
    print("[3/9] David Park (Ankylosing Spondylitis, 6 years)...")
    david_dir = os.path.join(DEMO_DATA_DIR, "david_park")
    generate_david(david_dir)
    generate_david_enhanced(david_dir)
    generate_david_dismissals(david_dir)
    generate_david_med_list(david_dir)
    generate_david_referrals(david_dir)
    generate_david_additional(david_dir)

    # ── Patient 4: Linda Chen -- 47F, 4-year Sjogren's odyssey ──
    print("[4/9] Linda Chen (Sjogren's Syndrome, 4 years)...")
    linda_dir = os.path.join(DEMO_DATA_DIR, "linda_chen")
    generate_linda(linda_dir)
    generate_linda_enhanced(linda_dir)
    generate_linda_dismissals(linda_dir)
    generate_linda_med_list(linda_dir)
    generate_linda_referrals(linda_dir)
    generate_linda_additional(linda_dir)

    # ── Patient 5: Rachel Thompson -- 32F, 3-year MCTD odyssey ──
    print("[5/9] Rachel Thompson (Mixed Connective Tissue Disease)...")
    rachel_dir = os.path.join(DEMO_DATA_DIR, "rachel_thompson")
    generate_rachel(rachel_dir)
    generate_rachel_enhanced(rachel_dir)
    generate_rachel_med_list(rachel_dir)
    generate_rachel_referrals(rachel_dir)
    generate_rachel_genetics(rachel_dir)

    # ── Patient 6: Emma Williams -- 34F, Multiple Sclerosis ──
    print("[6/9] Emma Williams (Multiple Sclerosis, 2 years)...")
    emma_dir = os.path.join(DEMO_DATA_DIR, "emma_williams")
    gen_emma(emma_dir)

    # ── Patient 7: James Cooper -- 19M, T1D + Celiac overlap ──
    print("[7/9] James Cooper (Type 1 Diabetes + Celiac Disease, 18 months)...")
    james_dir = os.path.join(DEMO_DATA_DIR, "james_cooper")
    gen_james(james_dir)

    # ── Patient 8: Karen Foster -- 48F, Systemic Sclerosis ──
    print("[8/9] Karen Foster (Systemic Sclerosis / Scleroderma, 3 years)...")
    karen_dir = os.path.join(DEMO_DATA_DIR, "karen_foster")
    gen_karen(karen_dir)

    # ── Patient 9: Michael Torres -- 41M, Graves' Disease ──
    print("[9/9] Michael Torres (Graves' Disease, 2 years)...")
    michael_dir = os.path.join(DEMO_DATA_DIR, "michael_torres")
    gen_michael(michael_dir)

    elapsed = time.time() - start

    # Count documents per patient and total
    print()
    print("=" * 70)
    patient_counts = {}
    total = 0
    for patient_dir in sorted(os.listdir(DEMO_DATA_DIR)):
        pdir = os.path.join(DEMO_DATA_DIR, patient_dir)
        if os.path.isdir(pdir):
            count = len([f for f in os.listdir(pdir) if f.endswith(".pdf")])
            patient_counts[patient_dir] = count
            total += count

    print(f"COMPLETE: {total} clinical PDFs generated in {elapsed:.1f}s")
    print(f"Output: {os.path.abspath(DEMO_DATA_DIR)}")
    print()
    print("Document Counts:")
    for name, count in patient_counts.items():
        print(f"  {name}: {count} PDFs")
    print()
    print("Patient Summary:")
    print("  Sarah Mitchell  -- SLE with Class IV lupus nephritis")
    print("                     3.5-year odyssey: PCP, UC, ER, Derm, Rheum, Nephrology")
    print("                     Includes: ER dismissal, urgent care misdiagnosis, referral chain")
    print("  Maya Rodriguez  -- POTS / hEDS / MCAS / Small Fiber Neuropathy")
    print("                     4-year odyssey: PCP, Psych, ER, Cardio, GI, Derm, Neuro, Allergy")
    print("                     Includes: Psychiatry GAD/panic misdiagnosis, rosacea misdiagnosis")
    print("  David Park      -- Ankylosing Spondylitis")
    print("                     6-year odyssey: PCP, PT, Ophth, Sleep, Ortho, Rheum, GI")
    print("                     Includes: Sleep study dead-end, ophthalmology warning letters")
    print("  Linda Chen      -- Primary Sjogren's Syndrome")
    print("                     4-year odyssey: Ophth, Dentistry, PCP, Rheum, Pulm")
    print("                     Includes: Perimenopause misattribution, dismissed ANA, provider letters")
    print("  Rachel Thompson -- Mixed Connective Tissue Disease")
    print("                     3-year journey: PCP, Rheum, Cardio, Pulm, Genetics")
    print("                     Includes: HLA/genetic risk panel, flare/remission cycle, PAH screening")
    print("  Emma Williams   -- Multiple Sclerosis (RRMS)")
    print("                     2-year journey: PCP, Ophth, Neurology, Radiology")
    print("                     Includes: Optic neuritis, MRI brain/spine, LP, Ocrelizumab")
    print("  James Cooper    -- Type 1 Diabetes + Celiac Disease (overlap)")
    print("                     18-month journey: ER (DKA), Endocrine, GI, Genetics, Dietitian")
    print("                     Includes: DKA admission, autoantibodies, duodenal biopsy, HLA typing")
    print("  Karen Foster    -- Systemic Sclerosis (Scleroderma)")
    print("                     3-year journey: PCP, Rheum, Derm, Pulm, Cardiology")
    print("                     Includes: Anti-Scl-70, ILD on HRCT, PFTs, PAH screening, mRSS")
    print("  Michael Torres  -- Graves' Disease")
    print("                     2-year journey: PCP, Endocrine, Cardio, Ophth, Nuclear Medicine")
    print("                     Includes: TRAb/TSI, thyroid scan, RAI therapy, TED evaluation")
    print()
    print("Document Types:")
    print("  Progress Notes  |  Lab Reports  |  Imaging  |  Pathology  |  Genetics")
    print("  Referral Letters  |  Medication Reconciliations  |  Dismissal Documentation")
    print("=" * 70)


if __name__ == "__main__":
    main()
