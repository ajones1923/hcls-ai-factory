#!/usr/bin/env bash
# Sync the repo's canonical skills into .claude/skills/ so Claude Code auto-discovers them.
#
# Uses real copies (not symlinks) for maximum discovery compatibility across Claude Code builds.
# The repo `skills/` tree is the single source of truth; re-run this after adding or editing a
# skill. `.claude/` is gitignored, so the synced copies stay local and never publish.
#
#   ./scripts/sync-skills.sh
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
src="$root/skills"
dest="$root/.claude/skills"
mkdir -p "$dest"

# Prune previously-synced skills (any dest subdir containing a SKILL.md) so deletions propagate.
if [ -d "$dest" ]; then
  find "$dest" -maxdepth 2 -name SKILL.md -printf '%h\0' 2>/dev/null | xargs -0 -r rm -rf
fi

# Copy every skill dir (a directory containing SKILL.md) from skills/ into a flat namespace.
count=0
while IFS= read -r skillmd; do
  d="$(dirname "$skillmd")"
  name="$(basename "$d")"
  rm -rf "$dest/$name"
  cp -R "$d" "$dest/$name"
  count=$((count + 1))
done < <(find "$src" -name SKILL.md)

echo "Synced $count skills -> $dest"
