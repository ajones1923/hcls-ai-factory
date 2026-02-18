"""
Security utilities for Genomics Pipeline Web Portal.

Provides:
- Security headers
- Rate limiting helpers
- Authentication helpers (optional)
"""
import hashlib
import hmac
import os
import secrets
from functools import wraps

from flask import Response, jsonify, request


def add_security_headers(response: Response) -> Response:
    """
    Add security headers to response.

    Args:
        response: Flask response object

    Returns:
        Response with security headers added
    """
    # Prevent clickjacking
    response.headers['X-Frame-Options'] = 'DENY'

    # Prevent MIME type sniffing
    response.headers['X-Content-Type-Options'] = 'nosniff'

    # Enable XSS filter in browsers
    response.headers['X-XSS-Protection'] = '1; mode=block'

    # Referrer policy
    response.headers['Referrer-Policy'] = 'strict-origin-when-cross-origin'

    # Content Security Policy (adjust as needed)
    response.headers['Content-Security-Policy'] = (
        "default-src 'self'; "
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://cdn.jsdelivr.net; "
        "style-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; "
        "font-src 'self' https://cdn.jsdelivr.net; "
        "img-src 'self' data:; "
        "connect-src 'self'"
    )

    # Permissions policy
    response.headers['Permissions-Policy'] = (
        'geolocation=(), microphone=(), camera=()'
    )

    return response


def generate_csrf_token() -> str:
    """Generate a CSRF token."""
    return secrets.token_urlsafe(32)


def verify_csrf_token(token: str, stored_token: str) -> bool:
    """Verify a CSRF token using constant-time comparison."""
    if not token or not stored_token:
        return False
    return hmac.compare_digest(token, stored_token)


class SimpleAuthenticator:
    """
    Simple API key authentication for the portal.

    This is optional - the portal can run without authentication
    for local/demo use, but authentication can be enabled for
    production deployments.
    """

    def __init__(self):
        self.enabled = False
        self.api_key = None
        self._load_config()

    def _load_config(self):
        """Load authentication configuration from environment."""
        api_key = os.environ.get('PORTAL_API_KEY')
        if api_key:
            self.enabled = True
            self.api_key = api_key

    def require_auth(self, f):
        """Decorator to require authentication for an endpoint."""
        @wraps(f)
        def decorated(*args, **kwargs):
            if not self.enabled:
                return f(*args, **kwargs)

            auth_header = request.headers.get('Authorization')
            if not auth_header:
                return jsonify({'error': 'Authentication required'}), 401

            # Check Bearer token
            if auth_header.startswith('Bearer '):
                token = auth_header[7:]
                if hmac.compare_digest(token, self.api_key):
                    return f(*args, **kwargs)

            return jsonify({'error': 'Invalid authentication'}), 401

        return decorated

    def is_enabled(self) -> bool:
        """Check if authentication is enabled."""
        return self.enabled


# Singleton authenticator
authenticator = SimpleAuthenticator()


def require_local_access(f):
    """
    Decorator to restrict access to localhost only.

    This is useful for sensitive operations like running
    pipeline steps or modifying configuration.
    """
    @wraps(f)
    def decorated(*args, **kwargs):
        # Allow localhost access
        client_ip = request.remote_addr
        allowed_ips = {'127.0.0.1', '::1', 'localhost'}

        if client_ip not in allowed_ips:
            # Check if running in Docker (allow Docker network)
            docker_network = os.environ.get('ALLOW_DOCKER_NETWORK', 'false').lower() == 'true'
            if docker_network and client_ip.startswith('172.'):
                return f(*args, **kwargs)

            return jsonify({'error': 'Access denied - localhost only'}), 403

        return f(*args, **kwargs)

    return decorated


class RateLimiter:
    """
    Simple in-memory rate limiter.

    For production use, consider using Flask-Limiter with Redis.
    """

    def __init__(self, max_requests: int = 100, window_seconds: int = 60):
        self.max_requests = max_requests
        self.window_seconds = window_seconds
        self.requests = {}  # ip -> [(timestamp, count), ...]

    def is_allowed(self, client_ip: str) -> bool:
        """Check if request is allowed."""
        import time

        current_time = time.time()
        window_start = current_time - self.window_seconds

        # Clean old entries
        if client_ip in self.requests:
            self.requests[client_ip] = [
                (ts, count) for ts, count in self.requests[client_ip]
                if ts > window_start
            ]
        else:
            self.requests[client_ip] = []

        # Count requests in window
        total = sum(count for _, count in self.requests[client_ip])

        if total >= self.max_requests:
            return False

        # Record this request
        self.requests[client_ip].append((current_time, 1))
        return True

    def limit(self, f):
        """Decorator to apply rate limiting."""
        @wraps(f)
        def decorated(*args, **kwargs):
            client_ip = request.remote_addr

            if not self.is_allowed(client_ip):
                return jsonify({
                    'error': 'Rate limit exceeded',
                    'retry_after': self.window_seconds
                }), 429

            return f(*args, **kwargs)

        return decorated


# Default rate limiter
rate_limiter = RateLimiter()
