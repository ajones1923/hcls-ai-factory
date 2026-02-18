"""Tests for genomics pipeline input validation."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))

from validation import (
    VALID_LOG_TYPES,
    VALID_STEPS,
    sanitize_path,
    validate_config_key,
    validate_config_value,
    validate_fastq_file,
    validate_log_type,
    validate_patient_id,
    validate_reference_file,
    validate_step_name,
)


class TestValidateStepName:
    """Tests for validate_step_name()."""

    @pytest.mark.parametrize("step", ['check', 'login', 'download', 'reference', 'test', 'full'])
    def test_valid_steps(self, step):
        valid, err = validate_step_name(step)
        assert valid is True
        assert err is None

    def test_invalid_step(self):
        valid, err = validate_step_name('exploit')
        assert valid is False
        assert 'Invalid step' in err

    def test_empty_step(self):
        valid, err = validate_step_name('')
        assert valid is False
        assert 'required' in err.lower()

    def test_none_step(self):
        valid, err = validate_step_name(None)
        assert valid is False

    def test_injection_attempt(self):
        valid, _ = validate_step_name('check; rm -rf /')
        assert valid is False

    def test_case_sensitive(self):
        valid, _ = validate_step_name('CHECK')
        assert valid is False

    def test_valid_steps_constant(self):
        assert {'check', 'login', 'download', 'reference', 'test', 'full'} == VALID_STEPS


class TestValidateLogType:
    """Tests for validate_log_type()."""

    @pytest.mark.parametrize("log_type", [
        'chr20_fq2bam', 'chr20_deepvariant', 'genome_fq2bam', 'genome_deepvariant'
    ])
    def test_valid_log_types(self, log_type):
        valid, err = validate_log_type(log_type)
        assert valid is True
        assert err is None

    def test_invalid_log_type(self):
        valid, err = validate_log_type('../../etc/passwd')
        assert valid is False

    def test_empty_log_type(self):
        valid, err = validate_log_type('')
        assert valid is False
        assert 'required' in err.lower()

    def test_none_log_type(self):
        valid, _ = validate_log_type(None)
        assert valid is False

    def test_valid_log_types_constant(self):
        assert {'chr20_fq2bam', 'chr20_deepvariant', 'genome_fq2bam', 'genome_deepvariant'} == VALID_LOG_TYPES


class TestValidateConfigKey:
    """Tests for validate_config_key()."""

    @pytest.mark.parametrize("key", ['API_KEY', 'MAX_RETRIES', 'GPU_COUNT', 'A1', 'NGC_API_KEY'])
    def test_valid_keys(self, key):
        valid, _ = validate_config_key(key)
        assert valid is True

    @pytest.mark.parametrize("key", ['lowercase', '1STARTS_NUM', 'has space', 'semi;colon', '_LEADING'])
    def test_invalid_keys(self, key):
        valid, _ = validate_config_key(key)
        assert valid is False

    def test_empty_key(self):
        valid, err = validate_config_key('')
        assert valid is False
        assert 'required' in err.lower()

    def test_none_key(self):
        valid, _ = validate_config_key(None)
        assert valid is False


class TestValidateConfigValue:
    """Tests for validate_config_value()."""

    @pytest.mark.parametrize("value", ['normal_value_123', '/path/to/file', 'true', '42', 'hello world'])
    def test_valid_values(self, value):
        valid, _ = validate_config_value(value)
        assert valid is True

    @pytest.mark.parametrize("value,desc", [
        ('val;ue', 'semicolon'),
        ('val&ue', 'ampersand'),
        ('val|ue', 'pipe'),
        ('val`cmd`', 'backtick'),
        ('val$HOME', 'dollar'),
        ('../etc/passwd', 'path traversal'),
    ])
    def test_dangerous_characters(self, value, desc):
        valid, err = validate_config_value(value)
        assert valid is False, f"Should reject {desc}"
        assert 'invalid characters' in err.lower()

    def test_none_value(self):
        valid, err = validate_config_value(None)
        assert valid is False
        assert 'required' in err.lower()


class TestSanitizePath:
    """Tests for sanitize_path()."""

    def test_valid_relative_path(self, tmp_path):
        result = sanitize_path('subdir/file.txt', tmp_path)
        assert result is not None
        assert str(result).startswith(str(tmp_path.resolve()))

    def test_path_traversal_blocked(self, tmp_path):
        result = sanitize_path('../../etc/passwd', tmp_path)
        assert result is None

    def test_simple_filename(self, tmp_path):
        result = sanitize_path('safe.txt', tmp_path)
        assert result is not None
        assert result.name == 'safe.txt'

    def test_absolute_traversal_blocked(self, tmp_path):
        result = sanitize_path('/etc/passwd', tmp_path)
        # /etc/passwd resolved won't be under tmp_path
        assert result is None

    def test_dot_dot_in_middle(self, tmp_path):
        # Create subdirectory so traversal can be tested
        (tmp_path / 'subdir').mkdir()
        result = sanitize_path('subdir/../../../etc/passwd', tmp_path)
        assert result is None


class TestValidatePatientId:
    """Tests for validate_patient_id()."""

    @pytest.mark.parametrize("pid", ['Patient001', 'HG002', 'sample-1_abc', 'A', 'test123'])
    def test_valid_ids(self, pid):
        valid, _ = validate_patient_id(pid)
        assert valid is True

    def test_empty_id(self):
        valid, err = validate_patient_id('')
        assert valid is False
        assert 'required' in err.lower()

    def test_none_id(self):
        valid, _ = validate_patient_id(None)
        assert valid is False

    def test_too_long_id(self):
        valid, err = validate_patient_id('x' * 51)
        assert valid is False
        assert 'too long' in err.lower()

    def test_exactly_max_length(self):
        valid, _ = validate_patient_id('x' * 50)
        assert valid is True

    def test_invalid_chars(self):
        valid, _ = validate_patient_id('patient;drop')
        assert valid is False

    def test_starts_with_hyphen(self):
        valid, _ = validate_patient_id('-invalid')
        assert valid is False

    def test_starts_with_underscore(self):
        valid, _ = validate_patient_id('_invalid')
        assert valid is False


class TestValidateFastqFile:
    """Tests for validate_fastq_file()."""

    def test_valid_fastq(self, tmp_path):
        f = tmp_path / 'sample.fastq'
        f.write_text('@SEQ\nACGT\n+\nIIII\n')
        valid, _ = validate_fastq_file(f)
        assert valid is True

    def test_valid_fq_extension(self, tmp_path):
        f = tmp_path / 'sample.fq'
        f.write_text('@SEQ\nACGT\n+\nIIII\n')
        valid, _ = validate_fastq_file(f)
        assert valid is True

    def test_valid_fastq_gz(self, tmp_path):
        f = tmp_path / 'sample.fastq.gz'
        f.write_bytes(b'\x1f\x8b' + b'\x00' * 10)
        valid, _ = validate_fastq_file(f)
        assert valid is True

    def test_missing_file(self, tmp_path):
        valid, err = validate_fastq_file(tmp_path / 'missing.fastq')
        assert valid is False
        assert 'not found' in err.lower()

    def test_wrong_extension(self, tmp_path):
        f = tmp_path / 'sample.bam'
        f.write_text('data')
        valid, err = validate_fastq_file(f)
        assert valid is False
        assert 'extension' in err.lower()

    def test_empty_file(self, tmp_path):
        f = tmp_path / 'empty.fastq'
        f.touch()
        valid, err = validate_fastq_file(f)
        assert valid is False
        assert 'empty' in err.lower()


class TestValidateReferenceFile:
    """Tests for validate_reference_file()."""

    def test_missing_reference(self, tmp_path):
        valid, err = validate_reference_file(tmp_path / 'ref.fa')
        assert valid is False
        assert 'not found' in err.lower()

    def test_wrong_extension(self, tmp_path):
        f = tmp_path / 'ref.bam'
        f.write_text('data')
        valid, err = validate_reference_file(f)
        assert valid is False
        assert 'extension' in err.lower()

    def test_missing_indices(self, tmp_path):
        f = tmp_path / 'ref.fa'
        f.write_text('>chr1\nACGT\n')
        valid, err = validate_reference_file(f)
        assert valid is False
        assert 'Missing index' in err

    def test_valid_reference_with_all_indices(self, tmp_path):
        f = tmp_path / 'ref.fa'
        f.write_text('>chr1\nACGT\n')
        for ext in ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']:
            (tmp_path / ('ref.fa' + ext)).touch()
        valid, _ = validate_reference_file(f)
        assert valid is True

    def test_partial_indices(self, tmp_path):
        f = tmp_path / 'ref.fa'
        f.write_text('>chr1\nACGT\n')
        # Only create some index files
        (tmp_path / 'ref.fa.fai').touch()
        (tmp_path / 'ref.fa.amb').touch()
        valid, err = validate_reference_file(f)
        assert valid is False
        assert 'Missing index' in err

    def test_valid_fasta_extension(self, tmp_path):
        f = tmp_path / 'ref.fasta'
        f.write_text('>chr1\nACGT\n')
        for ext in ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']:
            (tmp_path / ('ref.fasta' + ext)).touch()
        valid, _ = validate_reference_file(f)
        assert valid is True
