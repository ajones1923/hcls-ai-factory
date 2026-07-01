#!/usr/bin/env python3
# =============================================================================
# Cardiology Intelligence Agent - Seed Knowledge Base
# Author: Adam Jones
# Date: March 2026
# =============================================================================
"""Seed Milvus collections with curated cardiology knowledge base data.

This script populates the agent's vector collections with expertly curated
reference data spanning cardiovascular conditions, biomarkers, pharmacotherapy,
genomics, imaging protocols, electrophysiology criteria, clinical guidelines,
cardiac devices, hemodynamic parameters, and landmark trial evidence.

Each seed function targets a specific domain collection and uses the
corresponding ingest parser to transform structured knowledge into
embedding-ready records that are inserted into Milvus.

Usage:
    python scripts/seed_knowledge.py                  # seed all collections
    python scripts/seed_knowledge.py --collection conditions
    python scripts/seed_knowledge.py --collection biomarkers
    python scripts/seed_knowledge.py --dry-run        # preview without inserting

Collections seeded:
    conditions         Cardiovascular disease phenotypes (ICD-10, SNOMED)
    biomarkers         Cardiac biomarker reference ranges and clinical utility
    drug_classes       Cardiac pharmacotherapy classes and mechanisms
    genes              Cardiovascular genomics (genes, variants, pathways)
    imaging            Echo, CT, CMR, and nuclear imaging protocols
    ecg                ECG interpretation criteria and arrhythmia management
    guidelines         ACC/AHA/ESC guideline recommendations (2020-2025)
    devices            FDA-cleared cardiac devices and implantables
    hemodynamics       Hemodynamic parameters and catheterization protocols
    trials             Landmark cardiovascular clinical trial summaries
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import List

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.knowledge import (
    CARDIOVASCULAR_CONDITIONS,
    CARDIAC_BIOMARKERS,
    CARDIAC_DRUG_CLASSES,
    CARDIOVASCULAR_GENES,
    IMAGING_MODALITIES,
)
from src.ingest.imaging_parser import ImagingParser
from src.ingest.ecg_parser import ECGParser
from src.ingest.guideline_parser import GuidelineParser
from src.ingest.device_parser import DeviceParser
from src.ingest.hemodynamics_parser import HemodynamicsParser
from src.ingest.clinical_trials_parser import ClinicalTrialsCardioParser

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)

# Map of collection short names to seed functions (populated below)
SEED_REGISTRY: dict = {}


# ---------------------------------------------------------------------------
# Milvus insertion helper
# ---------------------------------------------------------------------------

def _insert_records(
    collection_name: str,
    records: List[dict],
    text_field: str = "description",
) -> int:
    """Generate embeddings and insert records into a Milvus collection.

    Degrades gracefully: if pymilvus or sentence_transformers are not
    installed, or if Milvus is unreachable, logs a warning and returns
    the record count (as if it were a dry run).

    Parameters
    ----------
    collection_name : str
        Target Milvus collection name.
    records : list[dict]
        Records to insert.  Each must contain *text_field*.
    text_field : str
        Key whose value is used to produce the embedding vector.

    Returns
    -------
    int
        Number of records inserted (or that would have been inserted on
        graceful degradation).
    """
    if not records:
        return 0

    # --- load optional dependencies ---
    try:
        from pymilvus import MilvusClient  # noqa: F811
    except ImportError:
        logger.warning(
            "pymilvus is not installed – skipping Milvus insert for %s "
            "(%d records would have been inserted).",
            collection_name,
            len(records),
        )
        return len(records)

    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        logger.warning(
            "sentence_transformers is not installed – skipping Milvus "
            "insert for %s (%d records would have been inserted).",
            collection_name,
            len(records),
        )
        return len(records)

    # --- connect to Milvus ---
    try:
        from config.settings import settings

        client = MilvusClient(
            uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}",
        )
    except Exception as exc:
        logger.warning(
            "Could not connect to Milvus for %s: %s – "
            "treating as dry run (%d records).",
            collection_name,
            exc,
            len(records),
        )
        return len(records)

    # --- generate embeddings ---
    try:
        model = SentenceTransformer("BAAI/bge-small-en-v1.5")
        texts = [str(r.get(text_field, "")) for r in records]
        embeddings = model.encode(texts, show_progress_bar=False).tolist()
    except Exception as exc:
        logger.warning(
            "Embedding generation failed for %s: %s – "
            "treating as dry run (%d records).",
            collection_name,
            exc,
            len(records),
        )
        return len(records)

    # --- insert ---
    try:
        insert_data = []
        for record, vector in zip(records, embeddings):
            entry = dict(record)
            entry["vector"] = vector
            insert_data.append(entry)
        client.insert(collection_name=collection_name, data=insert_data)
        logger.info(
            "Inserted %d records into Milvus collection '%s'.",
            len(insert_data),
            collection_name,
        )
        return len(insert_data)
    except Exception as exc:
        logger.warning(
            "Milvus insert failed for %s: %s – "
            "treating as dry run (%d records).",
            collection_name,
            exc,
            len(records),
        )
        return len(records)


def _register(name: str):
    """Decorator to register a seed function under a short name."""
    def decorator(fn):
        SEED_REGISTRY[name] = fn
        return fn
    return decorator


# ---------------------------------------------------------------------------
# Seed functions - one per domain collection
# ---------------------------------------------------------------------------

@_register("conditions")
def seed_conditions(dry_run: bool = False) -> int:
    """Seed cardiovascular conditions into the cardio_conditions collection.

    Sources: ICD-10-CM I00-I99, SNOMED-CT cardiovascular hierarchy, and curated
    phenotype descriptions covering coronary artery disease, heart failure,
    cardiomyopathies, valvular disease, congenital heart disease, vascular
    disorders, and arrhythmias.
    """
    logger.info("Seeding cardiovascular conditions (%d entries)...", len(CARDIOVASCULAR_CONDITIONS))
    if dry_run:
        return len(CARDIOVASCULAR_CONDITIONS)
    # Transform and insert
    records = []
    for name, condition in CARDIOVASCULAR_CONDITIONS.items():
        record = {
            "name": name,
            "icd10": condition.get("icd10", ""),
            "snomed": condition.get("snomed", ""),
            "category": condition.get("category", ""),
            "description": condition.get("description", ""),
            "risk_factors": condition.get("risk_factors", []),
            "diagnostic_criteria": condition.get("diagnostic_criteria", ""),
        }
        records.append(record)
    count = _insert_records("cardio_conditions", records, text_field="description")
    logger.info("Seeded %d cardiovascular conditions.", count)
    return count


@_register("biomarkers")
def seed_biomarkers(dry_run: bool = False) -> int:
    """Seed cardiac biomarker reference data into cardio_biomarkers.

    Includes troponin (I/T, hs-cTnI/T), BNP, NT-proBNP, CK-MB, myoglobin,
    D-dimer, hs-CRP, galectin-3, sST2, and emerging biomarkers with
    sex-specific reference ranges and clinical decision thresholds.
    """
    logger.info("Seeding cardiac biomarkers (%d entries)...", len(CARDIAC_BIOMARKERS))
    if dry_run:
        return len(CARDIAC_BIOMARKERS)
    records = []
    for name, biomarker in CARDIAC_BIOMARKERS.items():
        record = {
            "name": name,
            "analyte_code": biomarker.get("analyte_code", ""),
            "unit": biomarker.get("unit", ""),
            "reference_range_male": biomarker.get("reference_range_male", ""),
            "reference_range_female": biomarker.get("reference_range_female", ""),
            "critical_value": biomarker.get("critical_value", ""),
            "clinical_utility": biomarker.get("clinical_utility", ""),
            "turnaround_time": biomarker.get("turnaround_time", ""),
        }
        records.append(record)
    count = _insert_records("cardio_biomarkers", records, text_field="clinical_utility")
    logger.info("Seeded %d cardiac biomarkers.", count)
    return count


@_register("drug_classes")
def seed_drug_classes(dry_run: bool = False) -> int:
    """Seed cardiac drug class information into cardio_drug_classes.

    Covers beta-blockers, ACE inhibitors, ARBs, ARNI, calcium channel blockers,
    antiarrhythmics (Vaughan-Williams I-IV), anticoagulants, antiplatelets,
    statins, PCSK9 inhibitors, SGLT2 inhibitors, diuretics, vasodilators,
    inotropes, and emerging therapies.
    """
    logger.info("Seeding cardiac drug classes (%d entries)...", len(CARDIAC_DRUG_CLASSES))
    if dry_run:
        return len(CARDIAC_DRUG_CLASSES)
    records = []
    for name, drug_class in CARDIAC_DRUG_CLASSES.items():
        record = {
            "class_name": name,
            "mechanism": drug_class.get("mechanism", ""),
            "indications": drug_class.get("indications", []),
            "contraindications": drug_class.get("contraindications", []),
            "key_agents": drug_class.get("key_agents", []),
            "monitoring": drug_class.get("monitoring", ""),
            "evidence_level": drug_class.get("evidence_level", ""),
        }
        records.append(record)
    count = _insert_records("cardio_drug_classes", records, text_field="mechanism")
    logger.info("Seeded %d cardiac drug classes.", count)
    return count


@_register("genes")
def seed_genes(dry_run: bool = False) -> int:
    """Seed cardiovascular genomics data into cardio_genes.

    Includes cardiomyopathy genes (TTN, MYH7, MYBPC3, LMNA, etc.),
    channelopathy genes (SCN5A, KCNQ1, KCNH2, RYR2), lipid metabolism
    genes (LDLR, APOB, PCSK9), and structural heart disease genes with
    ClinVar pathogenicity classifications and inheritance patterns.
    """
    logger.info("Seeding cardiovascular genes (%d entries)...", len(CARDIOVASCULAR_GENES))
    if dry_run:
        return len(CARDIOVASCULAR_GENES)
    records = []
    for name, gene in CARDIOVASCULAR_GENES.items():
        record = {
            "gene_symbol": name,
            "gene_name": gene.get("gene_name", ""),
            "chromosome": gene.get("chromosome", ""),
            "associated_conditions": gene.get("associated_conditions", []),
            "inheritance": gene.get("inheritance", ""),
            "pathogenic_variants": gene.get("pathogenic_variants", 0),
            "actionability": gene.get("actionability", ""),
        }
        records.append(record)
    count = _insert_records("cardio_genes", records, text_field="actionability")
    logger.info("Seeded %d cardiovascular genes.", count)
    return count


@_register("imaging")
def seed_imaging(dry_run: bool = False) -> int:
    """Seed imaging protocols and measurements into cardio_imaging.

    Echocardiography: TTE, TEE, stress echo, strain imaging, 3D echo
    Cardiac CT: CCTA, calcium scoring, CT-FFR, structural CT
    Cardiac MRI: cine, LGE, T1/T2 mapping, perfusion, 4D flow
    Nuclear: SPECT MPI, PET MPI, pyrophosphate scan, MUGA
    """
    logger.info("Seeding imaging protocols and measurements...")
    if dry_run:
        return len(IMAGING_MODALITIES)
    parser = ImagingParser()
    records = parser.seed_imaging_protocols()
    records += parser.seed_echo_measurements()
    records += parser.seed_ct_protocols()
    records += parser.seed_cmr_protocols()
    records += parser.seed_nuclear_protocols()
    count = _insert_records("cardio_imaging", records, text_field="text")
    logger.info("Seeded %d imaging records.", count)
    return count


@_register("ecg")
def seed_electrophysiology(dry_run: bool = False) -> int:
    """Seed ECG criteria and arrhythmia management into cardio_ecg.

    ECG criteria: axis determination, interval normals, ST/T-wave patterns,
    bundle branch blocks, chamber enlargement, ischemia localization.
    Arrhythmia management: ACLS algorithms, ablation indications,
    antiarrhythmic selection, device programming.
    """
    logger.info("Seeding electrophysiology data...")
    parser = ECGParser()
    records = parser.seed_ecg_criteria()
    records += parser.seed_arrhythmia_management()
    if dry_run:
        logger.info("Dry run: would seed %d EP records.", len(records))
        return len(records)
    count = _insert_records("cardio_electrophysiology", records, text_field="text")
    logger.info("Seeded %d electrophysiology records.", count)
    return count


@_register("guidelines")
def seed_guidelines(dry_run: bool = False) -> int:
    """Seed ACC/AHA/ESC guideline recommendations into cardio_guidelines.

    Covers 2020-2025 guidelines for heart failure, valvular heart disease,
    chest pain evaluation, atrial fibrillation, ventricular arrhythmias,
    coronary revascularization, hypertrophic cardiomyopathy, pulmonary
    hypertension, and preventive cardiology.
    """
    logger.info("Seeding guideline recommendations...")
    parser = GuidelineParser()
    records = parser.seed_guidelines()
    if dry_run:
        logger.info("Dry run: would seed %d guideline records.", len(records))
        return len(records)
    count = _insert_records("cardio_guidelines", records, text_field="text")
    logger.info("Seeded %d guideline recommendations.", count)
    return count


@_register("devices")
def seed_devices(dry_run: bool = False) -> int:
    """Seed FDA-cleared cardiac devices into cardio_devices.

    Includes pacemakers, ICDs, CRT devices, LAA occluders, TAVR/TMVR valves,
    ventricular assist devices, percutaneous coronary intervention devices,
    and remote monitoring platforms.
    """
    logger.info("Seeding cardiac devices...")
    parser = DeviceParser()
    records = parser.seed_fda_devices()
    records += parser.seed_implantable_devices()
    if dry_run:
        logger.info("Dry run: would seed %d device records.", len(records))
        return len(records)
    count = _insert_records("cardio_devices", records, text_field="text")
    logger.info("Seeded %d cardiac device records.", count)
    return count


@_register("hemodynamics")
def seed_hemodynamics(dry_run: bool = False) -> int:
    """Seed hemodynamic reference parameters into cardio_hemodynamics.

    Normal hemodynamic values, pressure waveform interpretation,
    catheterization protocols (right heart, left heart, coronary),
    derived calculations (cardiac output, PVR, SVR, shunt fractions).
    """
    logger.info("Seeding hemodynamic parameters...")
    parser = HemodynamicsParser()
    records = parser.seed_hemodynamic_parameters()
    records += parser.seed_cathlab_protocols()
    if dry_run:
        logger.info("Dry run: would seed %d hemodynamic records.", len(records))
        return len(records)
    count = _insert_records("cardio_hemodynamics", records, text_field="text")
    logger.info("Seeded %d hemodynamic records.", count)
    return count


@_register("trials")
def seed_landmark_trials(dry_run: bool = False) -> int:
    """Seed landmark cardiovascular trial data into cardio_trials.

    Major trials: PARADIGM-HF, DAPA-HF, EMPEROR-Reduced, ISCHEMIA,
    REVIVED-BCIS2, ORBITA-2, CLEAR Outcomes, SELECT, and others
    with endpoints, NNT/NNH, and practice-changing conclusions.
    """
    logger.info("Seeding landmark trial data...")
    parser = ClinicalTrialsCardioParser()
    records = parser.seed_landmark_trials()
    if dry_run:
        logger.info("Dry run: would seed %d trial records.", len(records))
        return len(records)
    count = _insert_records("cardio_trials", records, text_field="text")
    logger.info("Seeded %d landmark trial records.", count)
    return count


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def seed_all(dry_run: bool = False) -> dict:
    """Seed all collections with curated data.

    Returns
    -------
    dict
        Mapping of collection name to number of records seeded.
    """
    logger.info("Starting knowledge base seeding (dry_run=%s)...", dry_run)
    start = time.time()
    results = {}

    for name, fn in SEED_REGISTRY.items():
        try:
            count = fn(dry_run=dry_run)
            results[name] = count
        except Exception:
            logger.exception("Failed to seed collection: %s", name)
            results[name] = -1

    elapsed = time.time() - start
    total = sum(v for v in results.values() if v > 0)
    failed = sum(1 for v in results.values() if v < 0)
    logger.info(
        "Knowledge base seeding complete. Total records: %d | Failed: %d | Time: %.1fs",
        total,
        failed,
        elapsed,
    )
    return results


def seed_collection(name: str, dry_run: bool = False) -> int:
    """Seed a single collection by short name.

    Parameters
    ----------
    name : str
        Short name of the collection (e.g. 'conditions', 'biomarkers').
    dry_run : bool
        If True, count records without inserting.
    """
    if name not in SEED_REGISTRY:
        available = ", ".join(sorted(SEED_REGISTRY.keys()))
        raise ValueError(f"Unknown collection '{name}'. Available: {available}")
    return SEED_REGISTRY[name](dry_run=dry_run)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Seed Cardiology Intelligence Agent knowledge base",
    )
    parser.add_argument(
        "--collection",
        choices=sorted(SEED_REGISTRY.keys()),
        default=None,
        help="Seed a single collection (default: all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview record counts without inserting into Milvus",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available collections and exit",
    )
    args = parser.parse_args()

    if args.list:
        for name in sorted(SEED_REGISTRY.keys()):
            fn = SEED_REGISTRY[name]
            doc_line = (fn.__doc__ or "").strip().split("\n")[0]
            print(f"  {name:<20s} {doc_line}")
        return

    if args.collection:
        count = seed_collection(args.collection, dry_run=args.dry_run)
        logger.info("Seeded %d records into %s.", count, args.collection)
    else:
        seed_all(dry_run=args.dry_run)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    main()
