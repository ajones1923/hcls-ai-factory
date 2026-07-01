#!/bin/bash
# Start Caddy HTTPS reverse proxy
cd "$(dirname "$0")"
caddy start --config Caddyfile
