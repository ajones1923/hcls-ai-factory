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
    python scripts/seed_knowledge.py --collection imaging
    python scripts/seed_knowledge.py --dry-run        # preview without inserting
    python scripts/seed_knowledge.py --force          # seed a non-empty collection (DUPLICATES)

Seed targets, and the Milvus collection each writes to:
    imaging        -> cardio_imaging            Echo, CT, CMR, and nuclear protocols
    ecg            -> cardio_electrophysiology  ECG criteria and arrhythmia management
    guidelines     -> cardio_guidelines         ACC/AHA/ESC recommendations (2020-2025)
    devices        -> cardio_devices            FDA-cleared devices and implantables
    hemodynamics   -> cardio_hemodynamics       Hemodynamic parameters and cath protocols
    trials         -> cardio_trials             Landmark cardiovascular trial summaries

    conditions     -> cardio_conditions         } NOT DECLARED COLLECTIONS. These four
    biomarkers     -> cardio_biomarkers         } targets are absent from
    drug_classes   -> cardio_drug_classes       } src/vector_collections.ALL_COLLECTIONS,
    genes          -> cardio_genes              } which is the authority on what exists.

The four undeclared targets hold 162 curated records with nowhere to go: setup_collections.py
will never create them, so seeding them logs "no such collection" and writes nothing. Whether
they should become collections (which changes every cardiology workflow's coverage denominator)
or map onto the declared twelve is a design decision -- see docs/build/EMPTY_COLLECTIONS.md.

Every collection above uses `auto_id=True`, so there is no key to upsert against and a second
run APPENDS. Seeding a non-empty collection is therefore refused unless --force.
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

from hcls_common.ingest_ledger import IngestEvent
from hcls_common.ingest_ledger import append as ledger_append
from hcls_common.ingest_ledger import model_conflict
from hcls_common.ingest_persist import DEFAULT_EMBED_MODEL as EMBED_MODEL
from hcls_common.ingest_persist import persist_records

from src.ingest.base import IngestRecord
from src.knowledge import (
    CARDIOVASCULAR_CONDITIONS,
    CARDIAC_BIOMARKERS,
    CARDIAC_DRUG_CLASSES,
    CARDIOVASCULAR_GENES,
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
#: Seeding an auto_id collection appends, so a second run duplicates the corpus. Set by
#: --force when that is what the operator actually wants.
FORCE = False

#: What actually happened per collection: name -> ("written" | "skipped" | "failed", rows).
#: The run summary is computed from this rather than from return values, because a seed
#: function returning 0 can mean "already seeded, correctly left alone" or "wrote nothing it
#: was asked to write", and reporting those as the same thing is how this script came to
#: print "Failed: 0" over a run that stored nothing.
OUTCOMES: dict = {}

#: Repo-relative, so a run from anywhere lands in the one committed ledger.
_LEDGER = Path(__file__).resolve().parents[4] / "docs" / "build" / "corpus_provenance.jsonl"


def _milvus_uri() -> str:
    try:
        from config.settings import settings
        return f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}"
    except Exception:
        import os
        return (f"http://{os.getenv('MILVUS_HOST', 'localhost')}:"
                f"{os.getenv('MILVUS_PORT', '19530')}")


def _existing_rows(collection_name: str) -> int:
    """Rows already in the collection, or -1 if that cannot be determined."""
    try:
        from pymilvus import MilvusClient
        client = MilvusClient(uri=_milvus_uri())
        stats = client.get_collection_stats(collection_name=collection_name)
        return int(stats.get("row_count", 0))
    except Exception:
        return -1


def _insert_records(
    collection_name: str,
    records: List[IngestRecord],
    text_field: str = "description",   # retained for call-site compatibility; unused
) -> int:
    """Embed and write records into a Milvus collection. Returns rows ACTUALLY written.

    This delegates to :func:`hcls_common.ingest_persist.persist_records` rather than doing
    its own insert. The reason is the history: this function used to flatten records itself,
    hard-code the vector column as ``vector`` when every cardio_* schema calls it
    ``embedding``, and -- in all five of its failure paths -- ``return len(records)``. A run
    therefore printed "Total records: 392 | Failed: 0" while writing nothing at all, and the
    cardiology corpus sat at 2% coverage looking like a collection nobody had got round to
    seeding rather than one whose seeder was broken.

    ``persist_records`` exists precisely to project a record onto whatever a collection
    declares, and each of its rules came from a real rejected insert: UTF-8 byte truncation
    for VARCHAR widths, int64 vs float, BOOL defaults, lists only into ARRAY columns, and
    auto_id primary keys (which all twelve cardio_* collections use). Re-deriving any of that
    here is how the divergence happened the first time.

    *text_field* is ignored: the target schema decides which column holds text, and the
    embedding comes from each record's own ``text``.
    """
    if not records:
        return 0

    # persist_records groups on `collection_name`, which IngestRecord does not have, so the
    # default below is what every row actually goes to. If a parser declares a different
    # target, that intent is being silently overridden -- say so.
    declared = {getattr(r, "collection", collection_name) for r in records}
    if declared - {collection_name}:
        logger.warning(
            "%s: records declare a different target (%s); writing all of them to %s.",
            collection_name, ", ".join(sorted(declared - {collection_name})), collection_name,
        )

    existing = _existing_rows(collection_name)
    if existing > 0 and not FORCE:
        logger.warning(
            "%s already holds %d rows -- SKIPPING (%d records not written). These collections "
            "use auto_id, so there is no key to upsert against and seeding again would "
            "duplicate the curated set. Re-run with --force if that is intended.",
            collection_name, existing, len(records),
        )
        OUTCOMES[collection_name] = ("skipped", existing)
        return 0

    # Refuse before writing, not after: two 384-dim models in one collection insert equally
    # happily and nothing raises, so the only moment this can be caught is here.
    conflict = model_conflict(collection_name, EMBED_MODEL, _LEDGER)
    if conflict and not FORCE:
        logger.error("%s -- NOTHING was written (%d records). Use --force to override.",
                     conflict, len(records))
        OUTCOMES[collection_name] = ("failed", 0)
        return 0

    written = persist_records(records, collection_name, milvus_uri=_milvus_uri(),
                              embed_model=EMBED_MODEL)
    ledger_append(
        IngestEvent(collection=collection_name, source="curated (cardiology seed_knowledge.py)",
                    query="seed", embed_model=EMBED_MODEL,
                    fetched=len(records), written=written,
                    skipped=len(records) - written),
        _LEDGER,
    )
    if written < len(records):
        logger.error(
            "%s: %d of %d records were NOT written. The log above says why; this run did "
            "not seed what it set out to.",
            collection_name, len(records) - written, len(records),
        )
    OUTCOMES[collection_name] = ("written" if written else "failed", written)
    return written


def _as_records(collection: str, source: str, text_key: str,
                rows: List[dict]) -> List[IngestRecord]:
    """Wrap curated dicts as IngestRecords, embedding on *text_key*.

    The parser-backed seed functions already return IngestRecords; these four built bare
    dicts, which meant two record shapes flowing into one writer. One shape is enough.
    """
    out = []
    for row in rows:
        out.append(IngestRecord(
            text=str(row.get(text_key, "")),
            metadata=dict(row),
            collection=collection,
            source=source,
        ))
    return out


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
    count = _insert_records(
        "cardio_conditions",
        _as_records("cardio_conditions", "curated (ICD-10-CM / SNOMED-CT)",
                    "description", records))
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
    count = _insert_records(
        "cardio_biomarkers",
        _as_records("cardio_biomarkers", "curated", "clinical_utility", records))
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
    count = _insert_records(
        "cardio_drug_classes",
        _as_records("cardio_drug_classes", "curated", "mechanism", records))
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
    count = _insert_records(
        "cardio_genes",
        _as_records("cardio_genes", "curated (ClinVar-derived)", "actionability", records))
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
    parser = ImagingParser()
    records = parser.seed_imaging_protocols()
    records += parser.seed_echo_measurements()
    records += parser.seed_ct_protocols()
    records += parser.seed_cmr_protocols()
    records += parser.seed_nuclear_protocols()
    if dry_run:
        # This used to return len(IMAGING_MODALITIES) -- the number of modalities, not of
        # records -- so the dry run reported 15 where the real run would write 54.
        logger.info("Dry run: would seed %d imaging records.", len(records))
        return len(records)
    count = _insert_records("cardio_imaging", records)
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
    count = _insert_records("cardio_electrophysiology", records)
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
    count = _insert_records("cardio_guidelines", records)
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
    count = _insert_records("cardio_devices", records)
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
    count = _insert_records("cardio_hemodynamics", records)
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
    count = _insert_records("cardio_trials", records)
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
    intended = {}

    for name, fn in SEED_REGISTRY.items():
        # Count first, so a zero written can be told from an empty source.
        try:
            intended[name] = fn(dry_run=True)
        except Exception:
            intended[name] = -1
        if dry_run:
            results[name] = intended[name]
            continue
        try:
            results[name] = fn(dry_run=False)
        except Exception:
            logger.exception("Failed to seed collection: %s", name)
            results[name] = -1

    elapsed = time.time() - start
    total = sum(v for v in results.values() if v > 0)

    # "Failed" used to count only the seed functions that RAISED, which is why a run that
    # stored nothing still printed "Failed: 0". Attribute it per collection instead, and
    # keep "skipped because already seeded" separate from "stored nothing it was given" --
    # collapsing those two is the same mistake one level up.
    failed = sorted(c for c, (what, _) in OUTCOMES.items() if what == "failed")
    skipped = sorted(c for c, (what, _) in OUTCOMES.items() if what == "skipped")
    logger.info(
        "Knowledge base seeding complete. Total records: %d | Failed: %d | Skipped: %d | "
        "Time: %.1fs",
        total, len(failed), len(skipped), elapsed,
    )
    if skipped and not dry_run:
        logger.info("Left alone (already seeded): %s", ", ".join(skipped))
    if failed and not dry_run:
        logger.error(
            "%d collection(s) stored NOTHING of what they were given: %s. "
            "Scroll up for the per-collection reason. A target absent from "
            "src/vector_collections.ALL_COLLECTIONS will never exist, however often this runs.",
            len(failed), ", ".join(failed),
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
        "--force",
        action="store_true",
        help="Seed a collection that already holds rows (these use auto_id, so this "
             "DUPLICATES rather than replaces)",
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
        return 0

    global FORCE
    FORCE = args.force

    if args.collection:
        count = seed_collection(args.collection, dry_run=args.dry_run)
        logger.info("Seeded %d records into %s.", count, args.collection)
        return 0 if (count > 0 or args.dry_run) else 1

    seed_all(dry_run=args.dry_run)
    # Exit status has to carry the outcome: a seeding run is normally invoked from a shell
    # script or a Makefile that checks nothing but $?, and this one returned 0 while writing
    # nothing for as long as it has existed. A skip is not a failure -- a re-run over an
    # already-seeded corpus is the expected case and must stay green.
    failed = [c for c, (what, _) in OUTCOMES.items() if what == "failed"]
    return 1 if failed and not args.dry_run else 0


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    raise SystemExit(main())
