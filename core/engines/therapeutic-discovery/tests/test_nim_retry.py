"""Tests for NIM client retry logic and SMILES validation."""
import pytest
from pathlib import Path
from unittest.mock import patch, Mock, MagicMock
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.nim_clients import MolMIMClient, DiffDockClient, NIMServiceConfig
import requests


@pytest.fixture
def molmim_client():
    """Create a MolMIMClient with fast retry settings."""
    config = NIMServiceConfig(port=8001, max_retries=3, timeout=5)
    client = MolMIMClient(config)
    client.health_checked = True  # Skip health check
    return client


@pytest.fixture
def diffdock_client():
    """Create a DiffDockClient with fast retry settings."""
    config = NIMServiceConfig(port=8002, max_retries=3, timeout=5)
    client = DiffDockClient(config)
    client.health_checked = True
    return client


class TestMolMIMRetry:
    """Tests for MolMIM client retry behavior."""

    @patch('src.nim_clients.requests.post')
    def test_succeeds_first_attempt(self, mock_post, molmim_client):
        """Test successful generation on first attempt."""
        mock_post.return_value = Mock(
            status_code=200,
            json=Mock(return_value={'molecules': [{'smiles': 'CCO', 'score': 0.9}]}),
            raise_for_status=Mock(),
        )
        result = molmim_client.generate('CCO', num_molecules=1)
        assert len(result) == 1
        assert result[0]['smiles'] == 'CCO'
        assert mock_post.call_count == 1

    @patch('src.nim_clients.time.sleep')
    @patch('src.nim_clients.requests.post')
    def test_retries_on_failure(self, mock_post, mock_sleep, molmim_client):
        """Test retry after transient failures."""
        # First 2 calls fail, 3rd succeeds
        mock_post.side_effect = [
            requests.ConnectionError("Connection refused"),
            requests.ConnectionError("Connection refused"),
            Mock(
                status_code=200,
                json=Mock(return_value={'molecules': [{'smiles': 'CCO', 'score': 0.9}]}),
                raise_for_status=Mock(),
            ),
        ]
        result = molmim_client.generate('CCO', num_molecules=1)
        assert len(result) == 1
        assert mock_post.call_count == 3
        # Verify exponential backoff: sleep(2), sleep(4)
        assert mock_sleep.call_count == 2
        mock_sleep.assert_any_call(2)
        mock_sleep.assert_any_call(4)

    @patch('src.nim_clients.time.sleep')
    @patch('src.nim_clients.requests.post')
    def test_raises_after_max_retries(self, mock_post, mock_sleep, molmim_client):
        """Test RuntimeError after exhausting all retries."""
        mock_post.side_effect = requests.ConnectionError("Connection refused")
        with pytest.raises(RuntimeError, match="MolMIM generation failed after 3 attempts"):
            molmim_client.generate('CCO', num_molecules=1)
        assert mock_post.call_count == 3

    @patch('src.nim_clients.time.sleep')
    @patch('src.nim_clients.requests.post')
    def test_exponential_backoff_timing(self, mock_post, mock_sleep, molmim_client):
        """Test that backoff follows 2^attempt pattern."""
        mock_post.side_effect = requests.ConnectionError("timeout")
        with pytest.raises(RuntimeError):
            molmim_client.generate('CCO')
        # For 3 retries: sleep(2^1=2), sleep(2^2=4)
        calls = [c.args[0] for c in mock_sleep.call_args_list]
        assert calls == [2, 4]

    @patch('src.nim_clients.requests.post')
    def test_http_error_triggers_retry(self, mock_post, molmim_client):
        """Test that HTTP 500 errors trigger retry."""
        error_response = Mock(status_code=500)
        error_response.raise_for_status.side_effect = requests.HTTPError("500 Server Error")
        success_response = Mock(
            status_code=200,
            json=Mock(return_value={'molecules': []}),
            raise_for_status=Mock(),
        )
        mock_post.side_effect = [error_response, success_response]

        with patch('src.nim_clients.time.sleep'):
            result = molmim_client.generate('CCO')
        assert mock_post.call_count == 2


class TestDiffDockRetry:
    """Tests for DiffDock client retry behavior."""

    @patch('src.nim_clients.requests.post')
    def test_succeeds_first_attempt(self, mock_post, diffdock_client):
        """Test successful docking on first attempt."""
        mock_post.return_value = Mock(
            status_code=200,
            json=Mock(return_value={'poses': [{'docking_score': -8.5, 'confidence': 0.8}]}),
            raise_for_status=Mock(),
        )
        result = diffdock_client.dock('ATOM...', 'CCO')
        assert len(result) == 1
        assert result[0]['docking_score'] == -8.5
        assert mock_post.call_count == 1

    @patch('src.nim_clients.time.sleep')
    @patch('src.nim_clients.requests.post')
    def test_retries_on_failure(self, mock_post, mock_sleep, diffdock_client):
        """Test retry after transient failures."""
        mock_post.side_effect = [
            requests.ConnectionError("Connection refused"),
            Mock(
                status_code=200,
                json=Mock(return_value={'poses': [{'docking_score': -7.0, 'confidence': 0.7}]}),
                raise_for_status=Mock(),
            ),
        ]
        result = diffdock_client.dock('ATOM...', 'CCO')
        assert len(result) == 1
        assert mock_post.call_count == 2

    @patch('src.nim_clients.time.sleep')
    @patch('src.nim_clients.requests.post')
    def test_raises_after_max_retries(self, mock_post, mock_sleep, diffdock_client):
        """Test RuntimeError after exhausting all retries."""
        mock_post.side_effect = requests.ConnectionError("Connection refused")
        with pytest.raises(RuntimeError, match="DiffDock docking failed after 3 attempts"):
            diffdock_client.dock('ATOM...', 'CCO')
        assert mock_post.call_count == 3


class TestSMILESValidation:
    """Tests for SMILES validation before NIM calls."""

    def test_molmim_rejects_invalid_smiles(self, molmim_client):
        """Test that invalid SMILES raises ValueError."""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not installed")

        with pytest.raises(ValueError, match="Invalid SMILES"):
            molmim_client.generate('not_a_valid_smiles_string!!!')

    def test_diffdock_rejects_invalid_smiles(self, diffdock_client):
        """Test that invalid SMILES raises ValueError."""
        try:
            from rdkit import Chem
        except ImportError:
            pytest.skip("RDKit not installed")

        with pytest.raises(ValueError, match="Invalid SMILES"):
            diffdock_client.dock('ATOM...', 'not_a_valid_smiles!!!')

    @patch('src.nim_clients.requests.post')
    def test_molmim_accepts_valid_smiles(self, mock_post, molmim_client):
        """Test that valid SMILES passes validation."""
        mock_post.return_value = Mock(
            status_code=200,
            json=Mock(return_value={'molecules': [{'smiles': 'CCO', 'score': 0.9}]}),
            raise_for_status=Mock(),
        )
        # CB-5083 seed compound — should not raise
        result = molmim_client.generate(
            'CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5'
        )
        assert mock_post.call_count == 1
