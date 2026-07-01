#!/usr/bin/env bash
# Install the dependency-free pre-commit guard into this clone.
#   ./scripts/install-hooks.sh
# For the full framework-based hooks, also run:  pip install pre-commit && pre-commit install
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
hook="$root/.git/hooks/pre-commit"
ln -sf ../../scripts/pre-commit-hook.sh "$hook"
chmod +x "$root/scripts/pre-commit-hook.sh"
echo "Installed pre-commit guard -> $hook"
