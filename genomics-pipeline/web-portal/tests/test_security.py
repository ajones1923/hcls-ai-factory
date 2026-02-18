"""Tests for genomics pipeline security utilities."""
import os
import sys
import time
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


@pytest.fixture
def app():
    """Create a minimal Flask app for testing security decorators."""
    from flask import Flask
    app = Flask(__name__)
    app.config['TESTING'] = True
    return app


class TestSecurityHeaders:
    """Tests for add_security_headers()."""

    def test_x_frame_options(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert response.headers['X-Frame-Options'] == 'DENY'

    def test_x_content_type_options(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert response.headers['X-Content-Type-Options'] == 'nosniff'

    def test_xss_protection(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert response.headers['X-XSS-Protection'] == '1; mode=block'

    def test_referrer_policy(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert response.headers['Referrer-Policy'] == 'strict-origin-when-cross-origin'

    def test_csp_header_present(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert 'Content-Security-Policy' in response.headers
            assert "default-src 'self'" in response.headers['Content-Security-Policy']

    def test_permissions_policy(self, app):
        from security import add_security_headers
        with app.test_request_context():
            response = app.make_response('OK')
            response = add_security_headers(response)
            assert 'Permissions-Policy' in response.headers
            assert 'geolocation=()' in response.headers['Permissions-Policy']


class TestCSRF:
    """Tests for CSRF token generation and verification."""

    def test_generate_token_returns_string(self):
        from security import generate_csrf_token
        token = generate_csrf_token()
        assert isinstance(token, str)
        assert len(token) > 20

    def test_generate_token_unique(self):
        from security import generate_csrf_token
        t1 = generate_csrf_token()
        t2 = generate_csrf_token()
        assert t1 != t2

    def test_verify_valid_token(self):
        from security import generate_csrf_token, verify_csrf_token
        token = generate_csrf_token()
        assert verify_csrf_token(token, token) is True

    def test_verify_invalid_token(self):
        from security import verify_csrf_token
        assert verify_csrf_token('wrong', 'correct') is False

    def test_verify_empty_token(self):
        from security import verify_csrf_token
        assert verify_csrf_token('', 'stored') is False

    def test_verify_empty_stored(self):
        from security import verify_csrf_token
        assert verify_csrf_token('token', '') is False

    def test_verify_none_token(self):
        from security import verify_csrf_token
        assert verify_csrf_token(None, 'stored') is False

    def test_verify_none_stored(self):
        from security import verify_csrf_token
        assert verify_csrf_token('token', None) is False


class TestSimpleAuthenticator:
    """Tests for SimpleAuthenticator."""

    def test_auth_disabled_by_default(self):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {}, clear=True):
            auth = SimpleAuthenticator()
            assert auth.is_enabled() is False

    def test_auth_enabled_with_key(self):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'test-key-123'}):
            auth = SimpleAuthenticator()
            assert auth.is_enabled() is True

    def test_require_auth_passes_when_disabled(self, app):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {}, clear=True):
            auth = SimpleAuthenticator()

            @app.route('/test-open')
            @auth.require_auth
            def test_view():
                return 'OK'

            with app.test_client() as client:
                resp = client.get('/test-open')
                assert resp.status_code == 200

    def test_require_auth_rejects_missing_token(self, app):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'secret'}):
            auth = SimpleAuthenticator()

            @app.route('/test-protected')
            @auth.require_auth
            def protected():
                return 'OK'

            with app.test_client() as client:
                resp = client.get('/test-protected')
                assert resp.status_code == 401

    def test_require_auth_accepts_valid_bearer(self, app):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'secret'}):
            auth = SimpleAuthenticator()

            @app.route('/test-valid')
            @auth.require_auth
            def valid_view():
                return 'OK'

            with app.test_client() as client:
                resp = client.get('/test-valid', headers={'Authorization': 'Bearer secret'})
                assert resp.status_code == 200

    def test_require_auth_rejects_wrong_bearer(self, app):
        from security import SimpleAuthenticator
        with patch.dict(os.environ, {'PORTAL_API_KEY': 'secret'}):
            auth = SimpleAuthenticator()

            @app.route('/test-wrong')
            @auth.require_auth
            def wrong_view():
                return 'OK'

            with app.test_client() as client:
                resp = client.get('/test-wrong', headers={'Authorization': 'Bearer wrong-key'})
                assert resp.status_code == 401


class TestRequireLocalAccess:
    """Tests for require_local_access decorator."""

    def test_allows_localhost(self, app):
        from flask import jsonify
        from security import require_local_access

        @app.route('/test-local')
        @require_local_access
        def local_view():
            return jsonify({'ok': True})

        with app.test_client() as client:
            # Flask test client uses 127.0.0.1 by default
            resp = client.get('/test-local')
            assert resp.status_code == 200


class TestRateLimiter:
    """Tests for RateLimiter."""

    def test_allows_under_limit(self):
        from security import RateLimiter
        limiter = RateLimiter(max_requests=5, window_seconds=60)
        for _ in range(5):
            assert limiter.is_allowed('127.0.0.1') is True

    def test_blocks_over_limit(self):
        from security import RateLimiter
        limiter = RateLimiter(max_requests=3, window_seconds=60)
        for _ in range(3):
            limiter.is_allowed('10.0.0.1')
        assert limiter.is_allowed('10.0.0.1') is False

    def test_different_ips_independent(self):
        from security import RateLimiter
        limiter = RateLimiter(max_requests=1, window_seconds=60)
        assert limiter.is_allowed('10.0.0.1') is True
        assert limiter.is_allowed('10.0.0.2') is True
        assert limiter.is_allowed('10.0.0.1') is False

    def test_window_expiry(self):
        from security import RateLimiter
        limiter = RateLimiter(max_requests=1, window_seconds=1)
        assert limiter.is_allowed('10.0.0.1') is True
        assert limiter.is_allowed('10.0.0.1') is False
        time.sleep(1.1)
        assert limiter.is_allowed('10.0.0.1') is True

    def test_exact_limit(self):
        from security import RateLimiter
        limiter = RateLimiter(max_requests=2, window_seconds=60)
        assert limiter.is_allowed('10.0.0.1') is True
        assert limiter.is_allowed('10.0.0.1') is True
        assert limiter.is_allowed('10.0.0.1') is False
