"""
Milvus Client - Vector database operations for genomic evidence.
"""
from typing import List, Optional, Dict, Any
import re
import numpy as np
from loguru import logger

from pymilvus import (
    connections,
    utility,
    Collection,
    CollectionSchema,
    FieldSchema,
    DataType,
)

from .vcf_parser import VariantEvidence


class MilvusClient:
    """
    Client for storing and retrieving genomic evidence from Milvus.
    """

    # Input validation patterns for filter expression safety
    _GENE_PATTERN = re.compile(r'^[A-Za-z0-9_\-]+$')
    _CHROM_PATTERN = re.compile(r'^(chr)?[0-9XYM]+$')

    def __init__(
        self,
        host: str = "localhost",
        port: int = 19530,
        collection_name: str = "genomic_evidence",
        embedding_dim: int = 384,
    ):
        self.host = host
        self.port = port
        self.collection_name = collection_name
        self.embedding_dim = embedding_dim
        self._collection: Optional[Collection] = None

    @staticmethod
    def _sanitize_gene(gene: str) -> str:
        """Validate gene name to prevent filter expression injection."""
        if not gene or not MilvusClient._GENE_PATTERN.match(gene):
            raise ValueError(f"Invalid gene name: {gene!r}")
        return gene

    @staticmethod
    def _sanitize_chrom(chrom: str) -> str:
        """Validate chromosome name to prevent filter expression injection."""
        if not chrom or not MilvusClient._CHROM_PATTERN.match(chrom):
            raise ValueError(f"Invalid chromosome: {chrom!r}")
        return chrom

    def connect(self) -> None:
        """Connect to Milvus server."""
        logger.info(f"Connecting to Milvus at {self.host}:{self.port}")
        connections.connect(
            alias="default",
            host=self.host,
            port=self.port,
        )
        logger.info("Connected to Milvus")

    def disconnect(self) -> None:
        """Disconnect from Milvus server."""
        connections.disconnect("default")
        logger.info("Disconnected from Milvus")

    def collection_exists(self) -> bool:
        """Check if collection exists."""
        return utility.has_collection(self.collection_name)

    def create_collection(self, drop_existing: bool = False) -> Collection:
        """Create the genomic evidence collection."""
        if drop_existing and self.collection_exists():
            logger.warning(f"Dropping existing collection: {self.collection_name}")
            utility.drop_collection(self.collection_name)

        if self.collection_exists():
            logger.info(f"Collection {self.collection_name} already exists")
            self._collection = Collection(self.collection_name)
            return self._collection

        logger.info(f"Creating collection: {self.collection_name}")

        # Define schema
        fields = [
            FieldSchema(
                name="id",
                dtype=DataType.VARCHAR,
                is_primary=True,
                max_length=200,
                description="Variant ID (chr_pos_ref_alt)"
            ),
            FieldSchema(
                name="embedding",
                dtype=DataType.FLOAT_VECTOR,
                dim=self.embedding_dim,
                description="Text embedding vector"
            ),
            FieldSchema(
                name="chrom",
                dtype=DataType.VARCHAR,
                max_length=10,
                description="Chromosome"
            ),
            FieldSchema(
                name="pos",
                dtype=DataType.INT64,
                description="Position"
            ),
            FieldSchema(
                name="ref",
                dtype=DataType.VARCHAR,
                max_length=500,
                description="Reference allele"
            ),
            FieldSchema(
                name="alt",
                dtype=DataType.VARCHAR,
                max_length=500,
                description="Alternate allele"
            ),
            FieldSchema(
                name="qual",
                dtype=DataType.FLOAT,
                description="Quality score"
            ),
            FieldSchema(
                name="gene",
                dtype=DataType.VARCHAR,
                max_length=50,
                description="Gene symbol"
            ),
            FieldSchema(
                name="consequence",
                dtype=DataType.VARCHAR,
                max_length=100,
                description="Variant consequence"
            ),
            FieldSchema(
                name="impact",
                dtype=DataType.VARCHAR,
                max_length=20,
                description="Impact level"
            ),
            FieldSchema(
                name="genotype",
                dtype=DataType.VARCHAR,
                max_length=10,
                description="Sample genotype"
            ),
            FieldSchema(
                name="text_summary",
                dtype=DataType.VARCHAR,
                max_length=2000,
                description="Human-readable summary"
            ),
            FieldSchema(
                name="clinical_significance",
                dtype=DataType.VARCHAR,
                max_length=200,
                description="Clinical significance"
            ),
            FieldSchema(
                name="rsid",
                dtype=DataType.VARCHAR,
                max_length=20,
                description="dbSNP rsID"
            ),
            FieldSchema(
                name="disease_associations",
                dtype=DataType.VARCHAR,
                max_length=500,
                description="Associated diseases/phenotypes"
            ),
            FieldSchema(
                name="am_pathogenicity",
                dtype=DataType.FLOAT,
                description="AlphaMissense pathogenicity score (0-1)"
            ),
            FieldSchema(
                name="am_class",
                dtype=DataType.VARCHAR,
                max_length=30,
                description="AlphaMissense classification"
            ),
        ]

        schema = CollectionSchema(
            fields=fields,
            description="Genomic variant evidence for RAG"
        )

        self._collection = Collection(
            name=self.collection_name,
            schema=schema,
        )

        # Create index for vector search
        logger.info("Creating vector index...")
        index_params = {
            "metric_type": "COSINE",
            "index_type": "IVF_FLAT",
            "params": {"nlist": 1024},
        }
        self._collection.create_index(
            field_name="embedding",
            index_params=index_params,
        )

        logger.info(f"Collection {self.collection_name} created with index")
        return self._collection

    def get_collection(self) -> Collection:
        """Get or create collection."""
        if self._collection is None:
            if self.collection_exists():
                self._collection = Collection(self.collection_name)
            else:
                self._collection = self.create_collection()
        return self._collection

    def insert(
        self,
        evidence_list: List[VariantEvidence],
        embeddings: np.ndarray,
    ) -> int:
        """Insert evidence with embeddings into collection."""
        collection = self.get_collection()

        # Prepare data for insertion
        data = [
            [e.variant_id for e in evidence_list],  # id
            embeddings.tolist(),  # embedding
            [e.chrom for e in evidence_list],  # chrom
            [e.pos for e in evidence_list],  # pos
            [e.ref for e in evidence_list],  # ref
            [e.alt for e in evidence_list],  # alt
            [e.qual for e in evidence_list],  # qual
            [e.gene or "" for e in evidence_list],  # gene
            [e.consequence or "" for e in evidence_list],  # consequence
            [e.impact or "" for e in evidence_list],  # impact
            [e.genotype or "" for e in evidence_list],  # genotype
            [e.text_summary[:2000] for e in evidence_list],  # text_summary
            [e.clinical_significance or "" for e in evidence_list],  # clinical_significance
            [e.rsid or "" for e in evidence_list],  # rsid
            [getattr(e, 'disease_associations', '') or "" for e in evidence_list],  # disease_associations
            [e.am_pathogenicity if e.am_pathogenicity is not None else -1.0 for e in evidence_list],  # am_pathogenicity
            [e.am_class or "" for e in evidence_list],  # am_class
        ]

        result = collection.insert(data)
        return result.insert_count

    def flush(self) -> None:
        """Flush data to disk."""
        collection = self.get_collection()
        collection.flush()
        logger.info("Collection flushed to disk")

    def load(self) -> None:
        """Load collection into memory for searching."""
        collection = self.get_collection()
        collection.load()
        logger.info("Collection loaded into memory")

    def search(
        self,
        query_embedding: np.ndarray,
        top_k: int = 10,
        score_threshold: float = 0.0,
        filter_expr: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Search for similar evidence.

        Args:
            query_embedding: Query vector
            top_k: Number of results to return
            score_threshold: Minimum similarity score (0-1 for cosine)
            filter_expr: Milvus filter expression (e.g., "chrom == 'chr7'")

        Returns:
            List of matching evidence with scores
        """
        collection = self.get_collection()

        search_params = {
            "metric_type": "COSINE",
            "params": {"nprobe": 16},
        }

        # Fields to return
        output_fields = [
            "chrom", "pos", "ref", "alt", "qual", "gene",
            "consequence", "impact", "genotype", "text_summary",
            "clinical_significance", "rsid", "am_pathogenicity", "am_class"
        ]

        results = collection.search(
            data=[query_embedding.tolist()],
            anns_field="embedding",
            param=search_params,
            limit=top_k,
            expr=filter_expr,
            output_fields=output_fields,
        )

        # Parse results
        evidence_results = []
        for hits in results:
            for hit in hits:
                score = hit.score  # Cosine similarity (0-1)
                if score < score_threshold:
                    continue

                # Get AlphaMissense score, convert -1.0 sentinel back to None
                am_score = hit.entity.get("am_pathogenicity")
                if am_score == -1.0:
                    am_score = None

                evidence_results.append({
                    "id": hit.id,
                    "score": score,
                    "chrom": hit.entity.get("chrom"),
                    "pos": hit.entity.get("pos"),
                    "ref": hit.entity.get("ref"),
                    "alt": hit.entity.get("alt"),
                    "qual": hit.entity.get("qual"),
                    "gene": hit.entity.get("gene"),
                    "consequence": hit.entity.get("consequence"),
                    "impact": hit.entity.get("impact"),
                    "genotype": hit.entity.get("genotype"),
                    "text_summary": hit.entity.get("text_summary"),
                    "clinical_significance": hit.entity.get("clinical_significance"),
                    "rsid": hit.entity.get("rsid"),
                    "am_pathogenicity": am_score,
                    "am_class": hit.entity.get("am_class") or None,
                })

        return evidence_results

    def search_by_gene(
        self,
        gene: str,
        top_k: int = 100,
    ) -> List[Dict[str, Any]]:
        """Search for all variants in a specific gene."""
        gene = self._sanitize_gene(gene)
        collection = self.get_collection()

        results = collection.query(
            expr=f'gene == "{gene}"',
            output_fields=[
                "chrom", "pos", "ref", "alt", "qual", "gene",
                "consequence", "impact", "genotype", "text_summary",
                "clinical_significance", "rsid", "am_pathogenicity", "am_class"
            ],
            limit=top_k,
        )

        return results

    def search_by_region(
        self,
        chrom: str,
        start: int,
        end: int,
        top_k: int = 100,
    ) -> List[Dict[str, Any]]:
        """Search for variants in a genomic region."""
        chrom = self._sanitize_chrom(chrom)
        start = int(start)
        end = int(end)
        if start < 0 or end < 0 or end < start:
            raise ValueError(f"Invalid region: start={start}, end={end}")
        collection = self.get_collection()

        results = collection.query(
            expr=f'chrom == "{chrom}" and pos >= {start} and pos <= {end}',
            output_fields=[
                "chrom", "pos", "ref", "alt", "qual", "gene",
                "consequence", "impact", "genotype", "text_summary",
                "clinical_significance", "rsid", "am_pathogenicity", "am_class"
            ],
            limit=top_k,
        )

        return results

    def get_stats(self) -> Dict[str, Any]:
        """Get collection statistics."""
        collection = self.get_collection()
        return {
            "name": self.collection_name,
            "num_entities": collection.num_entities,
            "schema": str(collection.schema),
        }

    def drop_collection(self) -> None:
        """Drop the collection."""
        if self.collection_exists():
            utility.drop_collection(self.collection_name)
            self._collection = None
            logger.info(f"Collection {self.collection_name} dropped")
