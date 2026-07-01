#!/usr/bin/env python3
"""
Ingest VCF variants into Milvus vector database.

Usage:
    python scripts/ingest_vcf.py [--vcf PATH] [--limit N] [--drop-existing]
"""
import argparse
import sys
from pathlib import Path
from tqdm import tqdm
from loguru import logger

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config.settings import settings
from src.vcf_parser import VCFParser
from src.annotator import VariantAnnotator, ClinVarAnnotator, AlphaMissenseAnnotator
from src.embedder import EvidenceEmbedder, CachedEmbedder
from src.milvus_client import MilvusClient


def parse_args():
    parser = argparse.ArgumentParser(description="Ingest VCF into Milvus")
    parser.add_argument(
        "--vcf",
        type=Path,
        default=settings.VCF_INPUT_PATH,
        help="Path to VCF file"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=settings.MAX_VARIANTS_TO_INGEST,
        help="Maximum variants to ingest (None = all)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=settings.INGESTION_BATCH_SIZE,
        help="Batch size for ingestion"
    )
    parser.add_argument(
        "--drop-existing",
        action="store_true",
        help="Drop existing collection before ingesting"
    )
    parser.add_argument(
        "--no-annotation",
        action="store_true",
        help="Skip annotation (faster but less info)"
    )
    parser.add_argument(
        "--clinvar",
        type=Path,
        default=Path("data/annotations/clinvar_variant_summary.txt.gz"),
        help="Path to ClinVar variant_summary.txt.gz file"
    )
    parser.add_argument(
        "--alphamissense",
        type=Path,
        default=Path("data/annotations/AlphaMissense_hg38.tsv.gz"),
        help="Path to AlphaMissense predictions file"
    )
    parser.add_argument(
        "--use-vep",
        action="store_true",
        help="Use VEP API instead of ClinVar (slower but more comprehensive)"
    )
    parser.add_argument(
        "--include-ref-calls",
        action="store_true",
        help="Include homozygous reference (0/0) calls (not recommended)"
    )
    parser.add_argument(
        "--annotated-only",
        action="store_true",
        help="Only ingest variants with ClinVar annotations (best for demos)"
    )
    parser.add_argument(
        "--use-cache",
        action="store_true",
        help="Use embedding cache for resumable ingestion"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Configure logging
    logger.remove()
    logger.add(sys.stderr, level="INFO")

    logger.info("=" * 60)
    logger.info("VCF to Milvus Ingestion")
    logger.info("=" * 60)

    # Validate VCF path
    if not args.vcf.exists():
        logger.error(f"VCF file not found: {args.vcf}")
        sys.exit(1)

    logger.info(f"VCF file: {args.vcf}")
    logger.info(f"Batch size: {args.batch_size}")
    logger.info(f"Max variants: {args.limit or 'all'}")

    # Initialize components
    logger.info("Initializing components...")

    # VCF Parser
    parser = VCFParser(
        vcf_path=args.vcf,
        min_qual=settings.MIN_VARIANT_QUAL,
        include_chromosomes=settings.INCLUDE_CHROMOSOMES,
        max_variants=args.limit,
        exclude_ref_calls=not args.include_ref_calls,
    )

    # Annotators - ClinVar for clinical significance, AlphaMissense for AI predictions
    clinvar_annotator = None
    alphamissense_annotator = None

    if not args.no_annotation:
        # ClinVar annotator
        if args.use_vep:
            logger.info("Initializing VEP annotator (slow, uses REST API)...")
            clinvar_annotator = VariantAnnotator(
                use_vep=True,
                species=settings.VEP_SPECIES,
                assembly=settings.VEP_ASSEMBLY,
            )
        elif args.clinvar.exists():
            logger.info(f"Initializing ClinVar annotator from: {args.clinvar}")
            clinvar_annotator = ClinVarAnnotator(
                clinvar_file=args.clinvar,
                assembly="GRCh38",
            )
            logger.info("Loading ClinVar database (this may take a minute)...")
            count = clinvar_annotator.load()
            stats = clinvar_annotator.get_stats()
            logger.info(f"Loaded {count:,} ClinVar variants")
            logger.info(f"Significance breakdown: {stats.get('by_significance', {})}")
        else:
            logger.warning(f"ClinVar file not found: {args.clinvar}")
            logger.info("ClinVar annotation will be skipped")

        # AlphaMissense annotator
        if args.alphamissense.exists():
            logger.info(f"Initializing AlphaMissense annotator from: {args.alphamissense}")
            alphamissense_annotator = AlphaMissenseAnnotator(
                alphamissense_file=args.alphamissense,
            )
            logger.info("Loading AlphaMissense database (this may take 2-3 minutes)...")
            logger.info("AlphaMissense: 71M+ predictions, requires ~8GB RAM")
            count = alphamissense_annotator.load()
            stats = alphamissense_annotator.get_stats()
            logger.info(f"Loaded {count:,} AlphaMissense predictions")
            logger.info(f"Class breakdown: {stats.get('by_class', {})}")
        else:
            logger.warning(f"AlphaMissense file not found: {args.alphamissense}")
            logger.info("AlphaMissense annotation will be skipped")
    else:
        logger.info("Skipping annotation (--no-annotation specified)")

    # Embedder
    if args.use_cache:
        embedder = CachedEmbedder(
            cache_path=settings.CACHE_DIR / "embeddings",
            model_name=settings.EMBEDDING_MODEL,
            batch_size=settings.EMBEDDING_BATCH_SIZE,
        )
    else:
        embedder = EvidenceEmbedder(
            model_name=settings.EMBEDDING_MODEL,
            batch_size=settings.EMBEDDING_BATCH_SIZE,
        )

    # Milvus Client
    milvus = MilvusClient(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
        collection_name=settings.MILVUS_COLLECTION,
        embedding_dim=embedder.dimension,
    )

    # Connect to Milvus
    logger.info(f"Connecting to Milvus at {settings.MILVUS_HOST}:{settings.MILVUS_PORT}...")
    try:
        milvus.connect()
    except Exception as e:
        logger.error(f"Failed to connect to Milvus: {e}")
        logger.error("Make sure Milvus is running: docker-compose up -d")
        sys.exit(1)

    # Create collection
    milvus.create_collection(drop_existing=args.drop_existing)

    # Count variants for progress bar
    logger.info("Counting variants...")
    # For large files, skip counting to save time
    total_variants = args.limit if args.limit else None

    # Process variants in batches
    logger.info("Processing variants...")
    batch = []
    total_ingested = 0
    total_clinvar_annotated = 0
    total_alphamissense_annotated = 0

    pbar = tqdm(
        parser.parse(),
        total=total_variants,
        desc="Processing variants",
        unit=" variants"
    )

    for variant in pbar:
        # Annotate with ClinVar
        if clinvar_annotator:
            try:
                variant = clinvar_annotator.annotate(variant)
                if variant.gene:
                    total_clinvar_annotated += 1
            except Exception as e:
                logger.debug(f"ClinVar annotation failed for {variant.variant_id}: {e}")

        # Annotate with AlphaMissense
        if alphamissense_annotator:
            try:
                variant = alphamissense_annotator.annotate(variant)
                if variant.am_pathogenicity is not None:
                    total_alphamissense_annotated += 1
            except Exception as e:
                logger.debug(f"AlphaMissense annotation failed for {variant.variant_id}: {e}")

        # Skip unannotated variants if --annotated-only is set
        if args.annotated_only and not variant.gene:
            continue

        batch.append(variant)

        # Process batch
        if len(batch) >= args.batch_size:
            # Generate embeddings
            embeddings = embedder.embed_evidence_batch(batch)

            # Insert into Milvus
            inserted = milvus.insert(batch, embeddings)
            total_ingested += inserted

            pbar.set_postfix({
                "ingested": total_ingested,
                "clinvar": total_clinvar_annotated,
                "alphamissense": total_alphamissense_annotated
            })

            batch = []

    # Process remaining batch
    if batch:
        embeddings = embedder.embed_evidence_batch(batch)
        inserted = milvus.insert(batch, embeddings)
        total_ingested += inserted

    # Flush and save cache
    logger.info("Flushing data to disk...")
    milvus.flush()

    if args.use_cache and hasattr(embedder, 'flush_cache'):
        embedder.flush_cache()

    # Print summary
    logger.info("=" * 60)
    logger.info("Ingestion Complete!")
    logger.info("=" * 60)
    logger.info(f"Total variants ingested: {total_ingested:,}")
    logger.info(f"Variants with ClinVar annotation: {total_clinvar_annotated:,}")
    logger.info(f"Variants with AlphaMissense prediction: {total_alphamissense_annotated:,}")

    stats = milvus.get_stats()
    logger.info(f"Collection size: {stats['num_entities']:,}")

    milvus.disconnect()
    logger.info("Done!")


if __name__ == "__main__":
    main()
