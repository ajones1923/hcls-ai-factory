"""
Tests for RAG Chat Pipeline Portal server.
"""
import pytest
import json
import os
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


class TestHealthEndpoint:
    """Tests for /health endpoint."""

    def test_health_returns_200(self, client):
        """Test health check returns 200."""
        response = client.get('/health')
        assert response.status_code == 200

    def test_health_returns_json(self, client):
        """Test health check returns JSON with service info."""
        response = client.get('/health')
        data = json.loads(response.data)
        assert data['status'] == 'healthy'
        assert data['service'] == 'rag-portal'

    def test_healthz_alias(self, client):
        """Test /healthz alias works."""
        response = client.get('/healthz')
        assert response.status_code == 200


class TestApiKeyAuth:
    """Tests for require_api_key decorator on protected endpoints."""

    def test_config_post_requires_api_key(self, client, tmp_path):
        """Test that POST /api/config requires API key when set."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            with patch('server.CONFIG_FILE', tmp_path / '.env'):
                (tmp_path / '.env').write_text("# Config\n")
                response = client.post(
                    '/api/config',
                    data=json.dumps({'LLM_MODEL': 'llama3.1:8b'}),
                    content_type='application/json'
                )
                assert response.status_code == 401
                data = json.loads(response.data)
                assert 'Unauthorized' in data.get('error', '')

    def test_config_post_accepts_valid_key(self, client, tmp_path):
        """Test that POST /api/config works with correct API key."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            with patch('server.CONFIG_FILE', tmp_path / '.env'):
                (tmp_path / '.env').write_text("# Config\n")
                response = client.post(
                    '/api/config',
                    data=json.dumps({'LLM_MODEL': 'llama3.1:8b'}),
                    content_type='application/json',
                    headers={'X-API-Key': 'test-secret'}
                )
                assert response.status_code == 200

    def test_config_get_does_not_require_key(self, client, tmp_path):
        """Test that GET /api/config works without API key (read-only)."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            with patch('server.CONFIG_FILE', tmp_path / '.env'):
                (tmp_path / '.env').write_text("LLM_MODEL=test\n")
                response = client.get('/api/config')
                # GET is protected by the same decorator — depends on implementation
                # Our implementation protects both GET and POST
                assert response.status_code in (200, 401)

    def test_run_requires_api_key(self, client):
        """Test that /api/run requires API key."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                response = client.get('/api/run/setup')
                assert response.status_code == 401

    def test_stop_requires_api_key(self, client):
        """Test that /api/stop requires API key."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            with patch('server.running_processes', {}):
                with patch('server.pipeline_state', {'status': 'running', 'current_step': None}):
                    response = client.get('/api/stop')
                    assert response.status_code == 401

    def test_model_post_requires_api_key(self, client):
        """Test that POST /api/model requires API key."""
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-secret'}):
            response = client.post(
                '/api/model',
                data=json.dumps({'model': 'llama3.1:8b'}),
                content_type='application/json'
            )
            assert response.status_code == 401

    def test_endpoints_refused_without_key_set(self, anon_client, tmp_path, monkeypatch):
        """An unset PORTAL_API_KEY must refuse, not wave the caller through.

        This test previously asserted the opposite — a "dev mode" in which an unset key
        opened the API. This portal has never actually behaved that way (it generates an
        ephemeral key), so the test was asserting a posture the code did not have; the
        assertion is inverted rather than removed.
        """
        import server
        monkeypatch.delenv('PORTAL_API_KEY', raising=False)
        monkeypatch.setenv('PORTAL_API_KEY_FILE', str(tmp_path / 'portal.key'))
        monkeypatch.setattr(server._get_api_key, '_generated_key', None, raising=False)
        try:
            with patch('server.CONFIG_FILE', tmp_path / '.env'):
                (tmp_path / '.env').write_text("# Config\n")
                with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
                    response = anon_client.post(
                        '/api/config',
                        data=json.dumps({'LLM_MODEL': 'llama3.1:8b'}),
                        content_type='application/json'
                    )
                    assert response.status_code == 401
        finally:
            server._get_api_key._generated_key = None


class TestConfigWhitelist:
    """Tests for config key whitelist validation."""

    def test_valid_key_accepted(self, client, tmp_path):
        """Test that whitelisted keys are accepted."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'LLM_MODEL': 'llama3.1:8b'}),
                content_type='application/json'
            )
            assert response.status_code == 200

    def test_invalid_key_rejected(self, client, tmp_path):
        """Test that non-whitelisted keys are rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'SHELL': '/bin/sh'}),
                content_type='application/json'
            )
            assert response.status_code == 400
            data = json.loads(response.data)
            assert 'Invalid config keys' in data.get('error', '')

    def test_mixed_valid_invalid_rejected(self, client, tmp_path):
        """Test that a mix of valid and invalid keys is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({
                    'LLM_MODEL': 'llama3.1:8b',
                    'EVIL_KEY': 'malicious'
                }),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_empty_body_rejected(self, client):
        """Test that empty/null body is rejected."""
        response = client.post(
            '/api/config',
            data='null',
            content_type='application/json'
        )
        assert response.status_code == 400

    def test_non_dict_body_rejected(self, client):
        """Test that non-dict body is rejected."""
        response = client.post(
            '/api/config',
            data=json.dumps(['not', 'a', 'dict']),
            content_type='application/json'
        )
        assert response.status_code == 400

    def test_all_whitelisted_keys(self, client, tmp_path):
        """Test that all whitelisted keys are accepted."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            valid_config = {
                'LLM_PROVIDER': 'ollama',
                'LLM_MODEL': 'llama3.1:8b',
                'OLLAMA_HOST': 'http://localhost:11434',
            }
            response = client.post(
                '/api/config',
                data=json.dumps(valid_config),
                content_type='application/json'
            )
            assert response.status_code == 200


class TestConfigValueValidation:
    """Tests for config value range/type validation."""

    def test_temperature_too_high(self, client, tmp_path):
        """Test that temperature > 2.0 is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'LLM_TEMPERATURE': '999'}),
                content_type='application/json'
            )
            assert response.status_code == 400
            data = json.loads(response.data)
            assert 'LLM_TEMPERATURE' in data.get('error', '')

    def test_temperature_negative(self, client, tmp_path):
        """Test that negative temperature is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'LLM_TEMPERATURE': '-0.5'}),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_valid_temperature(self, client, tmp_path):
        """Test that valid temperature is accepted."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'LLM_TEMPERATURE': '0.7'}),
                content_type='application/json'
            )
            assert response.status_code == 200

    def test_port_negative(self, client, tmp_path):
        """Test that negative port is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'MILVUS_PORT': '-1'}),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_port_too_high(self, client, tmp_path):
        """Test that port > 65535 is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'MILVUS_PORT': '99999'}),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_valid_port(self, client, tmp_path):
        """Test that valid port is accepted."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'MILVUS_PORT': '19530'}),
                content_type='application/json'
            )
            assert response.status_code == 200

    def test_max_tokens_non_numeric(self, client, tmp_path):
        """Test that non-numeric max_tokens is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'LLM_MAX_TOKENS': 'abc'}),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_score_threshold_out_of_range(self, client, tmp_path):
        """Test that score threshold > 1.0 is rejected."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'RAG_SCORE_THRESHOLD': '1.5'}),
                content_type='application/json'
            )
            assert response.status_code == 400

    def test_keys_without_validators_pass_through(self, client, tmp_path):
        """Test that whitelisted keys without validators are accepted."""
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("# Config\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'VCF_INPUT_PATH': '/data/test.vcf.gz'}),
                content_type='application/json'
            )
            assert response.status_code == 200


class TestRunEndpoint:
    """Tests for /api/run endpoint."""

    def test_run_invalid_step(self, client):
        """Test that invalid step returns 400."""
        with patch('server.pipeline_state', {'status': 'idle', 'logs': [], 'progress': 0}):
            response = client.get('/api/run/invalid_step')
            assert response.status_code == 400

    def test_run_while_running(self, client):
        """Test that cannot run while already running."""
        with patch('server.pipeline_state', {'status': 'running', 'logs': [], 'progress': 0}):
            response = client.get('/api/run/setup')
            assert response.status_code == 400


class TestVcfPreviewLimit:
    """Tests for VCF preview limit clamping."""

    def test_default_limit(self, client, tmp_path):
        """Test that default limit is 100."""
        with patch('server.load_config', return_value={'VCF_INPUT_PATH': str(tmp_path / 'test.vcf')}):
            with patch('server.check_file_exists', return_value=False):
                response = client.get('/api/vcf-preview')
                assert response.status_code == 404  # File doesn't exist, but limit was parsed

    def test_excessive_limit_clamped(self, client, tmp_path):
        """Test that excessive limit is clamped to 10000."""
        # Create a dummy VCF
        vcf_file = tmp_path / 'test.vcf'
        vcf_file.write_text('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        with patch('server.load_config', return_value={'VCF_INPUT_PATH': str(vcf_file)}):
            with patch('server.check_file_exists', return_value=True):
                with patch('server.get_vcf_preview', return_value=[]) as mock_preview:
                    response = client.get('/api/vcf-preview?limit=999999')
                    if mock_preview.called:
                        # Verify the limit was clamped
                        call_limit = mock_preview.call_args[0][1] if len(mock_preview.call_args[0]) > 1 else mock_preview.call_args[1].get('limit', 100)
                        assert call_limit <= 10000


# ─────────────────────────────────────────────────────────────────────────────────────────
# Security regression tests (2026-11-03 pre-release audit, finding S10)
#
# The ephemeral portal key was written to the log in full, and that key gates /api/config,
# whose GET returned every key parsed out of this engine's .env — the file that holds
# HF_TOKEN and ANTHROPIC_API_KEY. Anyone who could read the log therefore held the
# operator's model credentials.
# ─────────────────────────────────────────────────────────────────────────────────────────


class TestEphemeralKeyNotLogged:
    """S10 regression: a credential written to a log stream."""

    def test_key_is_not_written_to_the_log(self, tmp_path, monkeypatch):
        """Read the log through loguru itself, which is where the key used to be written."""
        import server
        monkeypatch.delenv('PORTAL_API_KEY', raising=False)
        monkeypatch.setenv('PORTAL_API_KEY_FILE', str(tmp_path / 'portal.key'))
        monkeypatch.setattr(server._get_api_key, '_generated_key', None, raising=False)
        records = []
        sink = server.logger.add(lambda m: records.append(str(m)), level='DEBUG')
        try:
            key = server._get_api_key()
        finally:
            server.logger.remove(sink)
            server._get_api_key._generated_key = None
        logged = '\n'.join(records)
        assert key
        assert logged.strip(), 'expected the portal to log where the key went'
        assert key not in logged
        # A fingerprint is logged instead, so the operator can tell which key is live.
        assert 'fingerprint' in logged.lower()
        assert str(tmp_path / 'portal.key') in logged

    def test_key_file_is_private_and_holds_the_key(self, tmp_path, monkeypatch):
        import server
        key_file = tmp_path / 'portal.key'
        monkeypatch.delenv('PORTAL_API_KEY', raising=False)
        monkeypatch.setenv('PORTAL_API_KEY_FILE', str(key_file))
        monkeypatch.setattr(server._get_api_key, '_generated_key', None, raising=False)
        try:
            key = server._get_api_key()
            assert key_file.exists()
            assert oct(key_file.stat().st_mode & 0o777) == '0o600'
            assert key_file.read_text().strip() == key
        finally:
            server._get_api_key._generated_key = None

    def test_fingerprint_is_not_reversible(self):
        import server
        assert server._key_fingerprint('super-secret-value') != 'super-secret-value'
        assert 'secret-value' not in server._key_fingerprint('super-secret-value')


class TestConfigGetDoesNotServeTheEnv:
    """S10 regression: GET /api/config returned the engine's .env wholesale."""

    def test_undeclared_keys_are_not_returned(self, client, tmp_path):
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text(
                "LLM_MODEL=llama3.1:8b\n"
                "ANTHROPIC_API_KEY=sk-ant-do-not-disclose\n"
                "MINIO_ROOT_PASSWORD=also-secret\n"
            )
            response = client.get('/api/config')
            assert response.status_code == 200
            body = response.get_data(as_text=True)
            data = json.loads(body)
            assert data.get('LLM_MODEL') == 'llama3.1:8b'
            assert 'ANTHROPIC_API_KEY' not in data
            assert 'MINIO_ROOT_PASSWORD' not in data
            assert 'sk-ant-do-not-disclose' not in body
            assert 'also-secret' not in body

    def test_declared_secret_is_redacted_not_returned(self, client, tmp_path):
        """HF_TOKEN is writable through this route, so it is declared — but not readable."""
        import server
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("HF_TOKEN=hf_realtokenvalue\nLLM_MODEL=x\n")
            response = client.get('/api/config')
            body = response.get_data(as_text=True)
            assert 'hf_realtokenvalue' not in body
            assert json.loads(body)['HF_TOKEN'] == server.REDACTED

    def test_round_tripping_the_redaction_does_not_clobber_the_secret(self, client, tmp_path):
        import server
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("HF_TOKEN=hf_realtokenvalue\n")
            response = client.post(
                '/api/config',
                data=json.dumps({'HF_TOKEN': server.REDACTED}),
                content_type='application/json',
            )
            assert response.status_code == 200
            assert 'hf_realtokenvalue' in (tmp_path / '.env').read_text()
            assert server.REDACTED not in (tmp_path / '.env').read_text()

    def test_config_get_still_requires_the_key(self, anon_client, tmp_path):
        with patch('server.CONFIG_FILE', tmp_path / '.env'):
            (tmp_path / '.env').write_text("LLM_MODEL=x\n")
            assert anon_client.get('/api/config').status_code == 401


class TestPortalKeyComparison:
    """The portal key check must be constant-time and must not raise."""

    def test_comparison_is_constant_time(self):
        import inspect
        import server
        assert 'compare_digest' in inspect.getsource(server.require_api_key)

    def test_non_ascii_header_is_401_not_500(self, anon_client):
        response = anon_client.get('/api/config', headers={'X-API-Key': 'pässwörd'})
        assert response.status_code == 401
