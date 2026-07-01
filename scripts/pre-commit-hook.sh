#!/usr/bin/env bash
# Dependency-free pre-commit guard for the HCLS AI Factory.
# Blocks three classes of accidental commit: oversized files, secrets, and
# proprietary/vendor (VAST / AI-OS / DASE) material that must stay out of the
# neutral open-source repo. Installed by scripts/install-hooks.sh.
set -u
MAXBYTES=$((5 * 1024 * 1024))   # 5 MB
fail=0

# Staged, added/copied/modified files only.
mapfile -t files < <(git diff --cached --name-only --diff-filter=ACM)

for f in "${files[@]}"; do
  [ -f "$f" ] || continue

  # 1) Oversized file (data / weights / archives slipping past .gitignore).
  sz=$(wc -c < "$f" 2>/dev/null || echo 0)
  if [ "$sz" -gt "$MAXBYTES" ]; then
    echo "BLOCK: $f is $((sz/1024/1024)) MB (> 5 MB). Data/weights belong outside git."
    fail=1
  fi

  # Only scan text files for the pattern checks.
  git grep -Iq . -- "$f" 2>/dev/null || continue

  # 2) Secrets.
  if grep -InE 'sk-ant-[a-zA-Z0-9_-]{20}|AKIA[0-9A-Z]{16}|-----BEGIN [A-Z]+ PRIVATE KEY|gh[pousr]_[a-zA-Z0-9]{36}' "$f" \
       | grep -vqE '\.example|placeholder|your[_-]|xxxx|<|EXAMPLE'; then
    echo "BLOCK: $f looks like it contains a secret/API key."
    fail=1
  fi

  # 3) Proprietary/vendor terms (keep the OSS repo neutral). \bVAST\b is
  #    word-bounded so drug names like SIMVASTATIN are never matched.
  #    The guard files below legitimately contain these terms (they define the guard).
  case "$f" in
    scripts/pre-commit-hook.sh|.gitignore|.pre-commit-config.yaml) continue ;;
  esac
  if grep -InE '\bVAST\b|\bAI[ _-]OS\b|DataStore|DataBase|DataEngine|DataSpace|InsightEngine|AgentEngine|\bDASE\b' "$f" >/dev/null 2>&1; then
    echo "BLOCK: $f contains VAST/AI-OS/DASE branding — this repo must stay neutral."
    echo "       (If this is a false positive, commit with --no-verify.)"
    fail=1
  fi
done

if [ "$fail" -ne 0 ]; then
  echo ""
  echo "Pre-commit checks failed. Fix the above, or bypass with:  git commit --no-verify"
  exit 1
fi
exit 0
