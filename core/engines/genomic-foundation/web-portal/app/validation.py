"""
Input validation utilities for Genomics Pipeline Web Portal.
"""
import re
from pathlib import Path
from typing import Optional, Tuple


# Valid step names for the pipeline
VALID_STEPS = {'check', 'login', 'download', 'reference', 'test', 'full'}

# Valid log types
VALID_LOG_TYPES = {'chr20_fq2bam', 'chr20_deepvariant', 'genome_fq2bam', 'genome_deepvariant'}

# Pattern for safe file paths (alphanumeric, underscore, hyphen, dot)
SAFE_PATH_PATTERN = re.compile(r'^[a-zA-Z0-9_\-\.]+$')


def validate_step_name(step: str) -> Tuple[bool, Optional[str]]:
    """
    Validate pipeline step name.

    Args:
        step: Step name to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not step:
        return False, "Step name is required"

    if step not in VALID_STEPS:
        return False, f"Invalid step. Valid steps: {', '.join(VALID_STEPS)}"

    return True, None


def validate_log_type(log_type: str) -> Tuple[bool, Optional[str]]:
    """
    Validate log type name.

    Args:
        log_type: Log type to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not log_type:
        return False, "Log type is required"

    if log_type not in VALID_LOG_TYPES:
        return False, f"Invalid log type. Valid types: {', '.join(VALID_LOG_TYPES)}"

    return True, None


def validate_config_key(key: str) -> Tuple[bool, Optional[str]]:
    """
    Validate configuration key name.

    Args:
        key: Configuration key to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not key:
        return False, "Key is required"

    # Config keys should be uppercase alphanumeric with underscores
    if not re.match(r'^[A-Z][A-Z0-9_]*$', key):
        return False, "Invalid key format. Must be uppercase with underscores only."

    return True, None


def validate_config_value(value: str) -> Tuple[bool, Optional[str]]:
    """
    Validate configuration value.

    Args:
        value: Configuration value to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if value is None:
        return False, "Value is required"

    # Check for potentially dangerous characters
    dangerous_patterns = [
        r'[;&|`$]',  # Shell metacharacters
        r'\.\.',     # Path traversal
    ]

    for pattern in dangerous_patterns:
        if re.search(pattern, str(value)):
            return False, "Value contains invalid characters"

    return True, None


def sanitize_path(path: str, base_dir: Path) -> Optional[Path]:
    """
    Sanitize and validate a file path.

    Args:
        path: Path to sanitize
        base_dir: Base directory that path must be within

    Returns:
        Sanitized Path if valid, None if invalid
    """
    try:
        # Resolve to absolute path
        resolved = (base_dir / path).resolve()

        # Ensure path is within base_dir (prevent path traversal)
        if not str(resolved).startswith(str(base_dir.resolve())):
            return None

        return resolved
    except Exception:
        return None


def validate_patient_id(patient_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate patient/sample ID.

    Args:
        patient_id: Patient ID to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not patient_id:
        return False, "Patient ID is required"

    # Patient IDs should be alphanumeric with optional underscores/hyphens
    if not re.match(r'^[A-Za-z0-9][A-Za-z0-9_\-]*$', patient_id):
        return False, "Invalid patient ID format"

    if len(patient_id) > 50:
        return False, "Patient ID too long (max 50 characters)"

    return True, None


def validate_fastq_file(filepath: Path) -> Tuple[bool, Optional[str]]:
    """
    Validate FASTQ file.

    Args:
        filepath: Path to FASTQ file

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not filepath.exists():
        return False, f"File not found: {filepath}"

    # Check extension
    valid_extensions = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}
    suffix = ''.join(filepath.suffixes)
    if suffix not in valid_extensions:
        return False, f"Invalid file extension: {suffix}. Expected: {', '.join(valid_extensions)}"

    # Check file size (should be non-empty)
    if filepath.stat().st_size == 0:
        return False, f"File is empty: {filepath}"

    return True, None


def validate_reference_file(filepath: Path) -> Tuple[bool, Optional[str]]:
    """
    Validate reference genome file.

    Args:
        filepath: Path to reference file

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not filepath.exists():
        return False, f"Reference file not found: {filepath}"

    # Check extension
    valid_extensions = {'.fa', '.fasta', '.fa.gz', '.fasta.gz'}
    suffix = ''.join(filepath.suffixes)
    if suffix not in valid_extensions:
        return False, f"Invalid file extension: {suffix}"

    # Check for index files
    index_extensions = ['.fai', '.amb', '.ann', '.bwt', '.pac', '.sa']
    missing_indices = []
    for ext in index_extensions:
        index_file = filepath.parent / (filepath.name + ext)
        if not index_file.exists():
            # Try without .gz if applicable
            base_name = filepath.name.replace('.gz', '')
            index_file_alt = filepath.parent / (base_name + ext)
            if not index_file_alt.exists():
                missing_indices.append(ext)

    if missing_indices:
        return False, f"Missing index files: {', '.join(missing_indices)}"

    return True, None
