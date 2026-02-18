"""Tests for Milvus client input sanitization and query safety."""
import re
import sys
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.milvus_client import MilvusClient


class TestSanitizeGene:
    """Test gene name sanitization."""

    @pytest.mark.parametrize("gene", [
        'BRCA1', 'TP53', 'C9orf72', 'BCR-ABL1', 'HLA_A', 'VCP', 'EGFR',
    ])
    def test_valid_gene_names(self, gene):
        assert MilvusClient._sanitize_gene(gene) == gene

    @pytest.mark.parametrize("gene", [
        '"; DROP',
        '" OR 1==1',
        'BRCA1"',
        'gene" or "1"=="1',
        'BRCA1; rm -rf',
        'gene name with spaces',
        'gene\ttab',
    ])
    def test_injection_attempts_rejected(self, gene):
        with pytest.raises(ValueError, match="Invalid gene name"):
            MilvusClient._sanitize_gene(gene)

    def test_empty_gene_rejected(self):
        with pytest.raises(ValueError):
            MilvusClient._sanitize_gene('')

    def test_none_gene_rejected(self):
        with pytest.raises(ValueError):
            MilvusClient._sanitize_gene(None)


class TestSanitizeChrom:
    """Test chromosome name sanitization."""

    @pytest.mark.parametrize("chrom", [
        'chr1', 'chr22', 'chrX', 'chrY', 'chrM', '1', '22', 'X', 'Y', 'M',
    ])
    def test_valid_chromosomes(self, chrom):
        assert MilvusClient._sanitize_chrom(chrom) == chrom

    @pytest.mark.parametrize("chrom", [
        'chr1"',
        'chr1"; DROP',
        'chr1 OR 1==1',
        'chrZ',
        'chromosome1',
    ])
    def test_injection_attempts_rejected(self, chrom):
        with pytest.raises(ValueError, match="Invalid chromosome"):
            MilvusClient._sanitize_chrom(chrom)

    def test_empty_chrom_rejected(self):
        with pytest.raises(ValueError):
            MilvusClient._sanitize_chrom('')

    def test_none_chrom_rejected(self):
        with pytest.raises(ValueError):
            MilvusClient._sanitize_chrom(None)


class TestSearchByGeneWithSanitization:
    """Test that search_by_gene uses sanitization."""

    def test_normal_gene_constructs_safe_filter(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client._collection.query.return_value = []
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        client.search_by_gene('BRCA1')

        call_args = client._collection.query.call_args
        assert call_args[1]['expr'] == 'gene == "BRCA1"'

    def test_injection_attempt_raises_before_query(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        with pytest.raises(ValueError):
            client.search_by_gene('" OR gene != "')

        # Verify the query was never executed
        client._collection.query.assert_not_called()


class TestSearchByRegionWithSanitization:
    """Test that search_by_region uses sanitization."""

    def test_normal_region_constructs_safe_filter(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client._collection.query.return_value = []
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        client.search_by_region('chr7', 1000, 2000)

        call_args = client._collection.query.call_args
        assert call_args[1]['expr'] == 'chrom == "chr7" and pos >= 1000 and pos <= 2000'

    def test_injection_in_chrom_raises(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        with pytest.raises(ValueError):
            client.search_by_region('chr7" or chrom != "', 0, 999999999)

        client._collection.query.assert_not_called()

    def test_negative_start_raises(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        with pytest.raises(ValueError):
            client.search_by_region('chr1', -1, 1000)

    def test_end_before_start_raises(self):
        client = MilvusClient.__new__(MilvusClient)
        client._collection = Mock()
        client.collection_name = "test"
        client.get_collection = Mock(return_value=client._collection)

        with pytest.raises(ValueError):
            client.search_by_region('chr1', 2000, 1000)
