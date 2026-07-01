#!/usr/bin/env bash
# Fetch the official Human Phenotype Ontology release and build the compact index the
# Phenome Mapper validates against. Idempotent; safe to re-run. The engine degrades
# gracefully (terms pass through "unverified") if this is never run, but running it once
# gives every emitted HPO term real ontology grounding.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/ref
if [ ! -s data/ref/hp.obo ]; then
  echo "downloading HPO release (hp.obo)…"
  curl -sL --max-time 120 https://purl.obolibrary.org/obo/hp.obo -o data/ref/hp.obo
fi
echo "building HPO index…"
venv/bin/python -m src.utils.hpo 2>/dev/null || python3 -m src.utils.hpo
