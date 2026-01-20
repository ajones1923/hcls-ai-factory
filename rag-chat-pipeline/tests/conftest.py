"""
Shared pytest fixtures for RAG Chat Pipeline tests.
"""
import pytest
import tempfile
import json
from pathlib import Path
from typing import List
from unittest.mock import Mock, MagicMock
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.vcf_parser import VariantEvidence


@pytest.fixture
def sample_variant() -> VariantEvidence:
    """Create a sample variant for testing."""
    return VariantEvidence(
        variant_id="chr9_123456_A_G",
        chrom="chr9",
        pos=123456,
        ref="A",
        alt="G",
        qual=99.0,
        filter_status="PASS",
        genotype="0/1",
        depth=50,
        gene="VCP",
        consequence="missense_variant",
        impact="HIGH",
        hgvs_c="c.123A>G",
        hgvs_p="p.Arg41Gly",
        clinical_significance="Pathogenic",
        rsid="rs12345",
        disease_associations="Frontotemporal dementia",
        source="GIAB WGS VCF (GRCh38)",
        tags=["giab", "wgs", "hg002", "grch38", "disease:Frontotemporal dementia"],
    )


@pytest.fixture
def sample_variants() -> List[VariantEvidence]:
    """Create a list of sample variants for testing."""
    return [
        VariantEvidence(
            variant_id="chr9_123456_A_G",
            chrom="chr9",
            pos=123456,
            ref="A",
            alt="G",
            qual=99.0,
            filter_status="PASS",
            genotype="0/1",
            depth=50,
            gene="VCP",
            consequence="missense_variant",
            impact="HIGH",
        ),
        VariantEvidence(
            variant_id="chr7_55259515_C_T",
            chrom="chr7",
            pos=55259515,
            ref="C",
            alt="T",
            qual=85.0,
            filter_status="PASS",
            genotype="0/1",
            depth=45,
            gene="EGFR",
            consequence="missense_variant",
            impact="MODERATE",
        ),
        VariantEvidence(
            variant_id="chr10_89692904_G_A",
            chrom="chr10",
            pos=89692904,
            ref="G",
            alt="A",
            qual=75.0,
            filter_status="PASS",
            genotype="1/1",
            depth=40,
            gene="CYP2C19",
            consequence="synonymous_variant",
            impact="LOW",
            clinical_significance="Benign",
        ),
    ]


@pytest.fixture
def temp_vcf_file():
    """Create a temporary VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr7,length=159345973>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr9\t123456\trs12345\tA\tG\t99\tPASS\t.\tGT:DP\t0/1:50
chr9\t123457\t.\tC\tT\t85\tPASS\t.\tGT:DP\t0/1:45
chr7\t55259515\t.\tC\tT\t75\tPASS\t.\tGT:DP\t1/1:40
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
        f.write(vcf_content)
        f.flush()
        yield Path(f.name)

    # Cleanup
    Path(f.name).unlink(missing_ok=True)


@pytest.fixture
def temp_hypotheses_dir():
    """Create a temporary directory for hypothesis storage."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def mock_embedder():
    """Create a mock embedder for testing."""
    embedder = Mock()
    embedder.embed_query.return_value = np.random.randn(384).astype(np.float32)
    embedder.embed_batch.return_value = np.random.randn(10, 384).astype(np.float32)
    embedder.dimension = 384
    return embedder


@pytest.fixture
def mock_milvus_client():
    """Create a mock Milvus client for testing."""
    client = Mock()
    client.connect.return_value = None
    client.collection_exists.return_value = True
    client.get_collection.return_value = Mock()
    client.search.return_value = [
        {
            "id": "chr9_123456_A_G",
            "score": 0.95,
            "gene": "VCP",
            "chrom": "chr9",
            "pos": 123456,
            "text_summary": "Variant at chr9:123456 A>G in gene VCP (missense_variant) with HIGH impact",
        }
    ]
    client.get_stats.return_value = {"num_entities": 1000, "name": "genomic_evidence"}
    return client


@pytest.fixture
def mock_llm_client():
    """Create a mock LLM client for testing."""
    client = Mock()
    client.generate.return_value = "This is a mock LLM response about genomic variants."

    def mock_stream(*args, **kwargs):
        tokens = ["This ", "is ", "a ", "streaming ", "response."]
        for token in tokens:
            yield token

    client.generate_stream.side_effect = mock_stream
    return client


@pytest.fixture
def sample_evidence_results():
    """Sample evidence results from Milvus search."""
    return [
        {
            "id": "chr9_123456_A_G",
            "score": 0.95,
            "gene": "VCP",
            "chrom": "chr9",
            "pos": 123456,
            "ref": "A",
            "alt": "G",
            "consequence": "missense_variant",
            "impact": "HIGH",
            "clinical_significance": "Pathogenic",
            "text_summary": "Variant at chr9:123456 A>G in gene VCP (missense_variant) with HIGH impact",
        },
        {
            "id": "chr9_123789_G_A",
            "score": 0.82,
            "gene": "C9orf72",
            "chrom": "chr9",
            "pos": 123789,
            "ref": "G",
            "alt": "A",
            "consequence": "intronic_variant",
            "impact": "MODIFIER",
            "clinical_significance": None,
            "text_summary": "Variant at chr9:123789 G>A in gene C9orf72 (intronic_variant)",
        },
    ]
