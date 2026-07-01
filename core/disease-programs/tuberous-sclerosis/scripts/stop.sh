#!/usr/bin/env bash
# Stop the host-run TSC Intelligence Engine services started by serve.sh.
cd "$(dirname "$0")/.."
for s in api briefing dashboard alerts; do
  f="data/state/$s.pid"
  if [ -f "$f" ] && kill "$(cat "$f")" 2>/dev/null; then
    echo "stopped $s (pid $(cat "$f"))"
  fi
  rm -f "$f"
done
