"""
Shared pytest fixtures for Genomics Pipeline Web Portal tests.
"""
import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
import sys

# Add the app directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


@pytest.fixture
def app():
    """Create Flask test app."""
    from server import app as flask_app
    flask_app.config['TESTING'] = True
    return flask_app


@pytest.fixture
def client(app):
    """Create test client."""
    return app.test_client()


@pytest.fixture
def mock_subprocess():
    """Mock subprocess for testing script execution."""
    with patch('server.subprocess') as mock:
        mock.run.return_value = Mock(returncode=0, stdout='', stderr='')
        mock.Popen.return_value = Mock(
            returncode=0,
            stdout=iter(['']),
            wait=Mock()
        )
        yield mock


@pytest.fixture
def mock_pynvml():
    """Mock pynvml for GPU testing without actual GPU."""
    mock_nvml = Mock()
    mock_nvml.nvmlInit.return_value = None
    mock_nvml.nvmlDeviceGetCount.return_value = 1
    mock_nvml.nvmlDeviceGetHandleByIndex.return_value = Mock()
    mock_nvml.nvmlDeviceGetName.return_value = "NVIDIA GB10 Test"
    mock_nvml.nvmlDeviceGetUtilizationRates.return_value = Mock(gpu=50, memory=30)
    mock_nvml.nvmlDeviceGetMemoryInfo.return_value = Mock(
        used=8 * 1024**3,
        total=128 * 1024**3
    )
    mock_nvml.nvmlDeviceGetTemperature.return_value = 45
    mock_nvml.nvmlDeviceGetPowerUsage.return_value = 150000  # 150W
    mock_nvml.nvmlDeviceGetComputeRunningProcesses.return_value = []
    mock_nvml.nvmlShutdown.return_value = None
    return mock_nvml


@pytest.fixture
def temp_data_dir():
    """Create temporary data directory structure."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create directory structure
        (tmpdir / 'input').mkdir()
        (tmpdir / 'ref').mkdir()
        (tmpdir / 'output' / 'logs').mkdir(parents=True)

        yield tmpdir


@pytest.fixture
def mock_psutil():
    """Mock psutil for system metrics."""
    with patch('server.psutil') as mock:
        def cpu_percent_side_effect(interval=None, percpu=False):
            if percpu:
                return [25.5, 20.0, 30.0, 15.0]
            return 25.5
        mock.cpu_percent.side_effect = cpu_percent_side_effect
        mock.cpu_freq.return_value = Mock(current=3200)
        mock.cpu_count.return_value = 128
        mock.virtual_memory.return_value = Mock(
            percent=45.0,
            used=60 * 1024**3,
            total=128 * 1024**3,
            available=68 * 1024**3
        )
        mock.swap_memory.return_value = Mock(percent=5.0)
        yield mock
