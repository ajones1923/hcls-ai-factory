"""Tests for hcls_common.security — input validation and sanitization."""

import pytest

from hcls_common.security import (
    sanitize_search_query,
    validate_milvus_filter,
    sanitize_gene_name,
    sanitize_chromosome,
    validate_patient_id,
    validate_smiles,
    validate_pdb_id,
    rate_limit_check,
    add_security_headers,
    sanitize_filename,
    _rate_windows,
)


class TestSanitizeSearchQuery:
    def test_basic_query_passthrough(self):
        assert sanitize_search_query("BRCA1 pathogenic") == "BRCA1 pathogenic"

    def test_strips_whitespace(self):
        assert sanitize_search_query("  BRCA1  ") == "BRCA1"

    def test_removes_shell_metacharacters(self):
        result = sanitize_search_query("BRCA1; rm -rf /")
        assert ";" not in result
        assert "|" not in sanitize_search_query("gene | cat /etc/passwd")

    def test_removes_path_traversal(self):
        result = sanitize_search_query("../../../etc/passwd")
        assert ".." not in result

    def test_truncates_to_max_length(self):
        long_query = "A" * 600
        result = sanitize_search_query(long_query, max_length=500)
        assert len(result) <= 500

    def test_empty_query_raises(self):
        with pytest.raises(ValueError, match="empty"):
            sanitize_search_query("")

    def test_none_query_raises(self):
        with pytest.raises(ValueError):
            sanitize_search_query(None)

    def test_only_metacharacters_raises(self):
        with pytest.raises(ValueError, match="invalid characters"):
            sanitize_search_query(";&|`$\\")


class TestValidateMilvusFilter:
    def test_valid_gene_filter(self):
        expr = 'gene == "BRCA1"'
        assert validate_milvus_filter(expr) == expr

    def test_valid_compound_filter(self):
        expr = 'gene == "TP53" and impact == "HIGH"'
        assert validate_milvus_filter(expr) == expr

    def test_valid_numeric_filter(self):
        expr = "pos >= 100000 and pos <= 200000"
        assert validate_milvus_filter(expr) == expr

    def test_rejects_injection_attempt(self):
        with pytest.raises(ValueError, match="disallowed"):
            validate_milvus_filter('gene == "BRCA1"; DROP COLLECTION')

    def test_empty_filter_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_milvus_filter("")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_milvus_filter("   ")


class TestSanitizeGeneName:
    def test_valid_gene(self):
        assert sanitize_gene_name("BRCA1") == "BRCA1"

    def test_valid_gene_with_hyphen(self):
        assert sanitize_gene_name("HLA-A") == "HLA-A"

    def test_valid_gene_with_underscore(self):
        assert sanitize_gene_name("TP53_v1") == "TP53_v1"

    def test_strips_whitespace(self):
        assert sanitize_gene_name("  BRCA2  ") == "BRCA2"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            sanitize_gene_name("")

    def test_too_long_raises(self):
        with pytest.raises(ValueError, match="too long"):
            sanitize_gene_name("A" * 51)

    def test_special_chars_raise(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_gene_name("BRCA1; DROP TABLE")


class TestSanitizeChromosome:
    def test_numeric_chromosome(self):
        assert sanitize_chromosome("1") == "1"
        assert sanitize_chromosome("22") == "22"

    def test_chr_prefix(self):
        assert sanitize_chromosome("chr1") == "chr1"
        assert sanitize_chromosome("chrX") == "chrX"

    def test_sex_chromosomes(self):
        assert sanitize_chromosome("X") == "X"
        assert sanitize_chromosome("Y") == "Y"

    def test_mitochondrial(self):
        assert sanitize_chromosome("MT") == "MT"
        assert sanitize_chromosome("M") == "M"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            sanitize_chromosome("")

    def test_invalid_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_chromosome("chromosome_one")


class TestValidatePatientId:
    def test_valid_id(self):
        valid, err = validate_patient_id("HG002")
        assert valid is True
        assert err is None

    def test_valid_id_with_underscore(self):
        valid, err = validate_patient_id("PATIENT_001")
        assert valid is True

    def test_empty_id(self):
        valid, err = validate_patient_id("")
        assert valid is False
        assert "required" in err

    def test_too_long(self):
        valid, err = validate_patient_id("A" * 51)
        assert valid is False
        assert "too long" in err

    def test_invalid_start_char(self):
        valid, err = validate_patient_id("-invalid")
        assert valid is False


class TestValidateSmiles:
    def test_valid_simple_smiles(self):
        valid, err = validate_smiles("CCO")  # ethanol
        assert valid is True
        assert err is None

    def test_empty_smiles(self):
        valid, err = validate_smiles("")
        assert valid is False
        assert "empty" in err

    def test_too_long_smiles(self):
        valid, err = validate_smiles("C" * 5001)
        assert valid is False
        assert "too long" in err

    def test_shell_injection_rejected(self):
        valid, err = validate_smiles("CC; rm -rf /")
        assert valid is False
        assert "disallowed" in err


class TestValidatePdbId:
    def test_valid_pdb_id(self):
        valid, err = validate_pdb_id("1A2B")
        assert valid is True
        assert err is None

    def test_empty_pdb_id(self):
        valid, err = validate_pdb_id("")
        assert valid is False

    def test_too_short(self):
        valid, err = validate_pdb_id("1AB")
        assert valid is False

    def test_too_long(self):
        valid, err = validate_pdb_id("1A2B3")
        assert valid is False

    def test_strips_whitespace_and_uppercases(self):
        valid, err = validate_pdb_id("  1a2b  ")
        assert valid is True


class TestRateLimitCheck:
    def test_allows_under_limit(self):
        # Use a unique key to avoid collisions with other tests
        key = "test_rate_limit_under"
        _rate_windows.pop(key, None)
        assert rate_limit_check(key, max_requests=5, window_seconds=60) is True

    def test_blocks_over_limit(self):
        key = "test_rate_limit_over"
        _rate_windows.pop(key, None)
        for _ in range(3):
            rate_limit_check(key, max_requests=3, window_seconds=60)
        assert rate_limit_check(key, max_requests=3, window_seconds=60) is False

    def test_window_expiry(self):
        key = "test_rate_limit_expiry"
        _rate_windows.pop(key, None)
        # Fill up window with a tiny window_seconds
        for _ in range(3):
            rate_limit_check(key, max_requests=3, window_seconds=0)
        # Expired entries should be purged
        assert rate_limit_check(key, max_requests=3, window_seconds=0) is True


class TestAddSecurityHeaders:
    def test_adds_headers_to_mock_response(self):
        class MockResponse:
            def __init__(self):
                self.headers = {}

        resp = MockResponse()
        result = add_security_headers(resp)
        assert result is resp
        assert resp.headers["X-Frame-Options"] == "DENY"
        assert resp.headers["X-Content-Type-Options"] == "nosniff"
        assert "Content-Security-Policy" in resp.headers
        assert "Permissions-Policy" in resp.headers


class TestSanitizeFilename:
    def test_basic_filename(self):
        assert sanitize_filename("report.pdf") == "report.pdf"

    def test_strips_directory_components(self):
        assert sanitize_filename("/etc/passwd") == "passwd"
        assert sanitize_filename("../../secret.txt") == "secret.txt"

    def test_replaces_unsafe_characters(self):
        result = sanitize_filename("my file (1).txt")
        assert " " not in result
        assert "(" not in result

    def test_truncates_long_filename(self):
        result = sanitize_filename("A" * 300, max_length=255)
        assert len(result) <= 255

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            sanitize_filename("")

    def test_all_unsafe_raises(self):
        with pytest.raises(ValueError, match="unsafe"):
            sanitize_filename("   ")
