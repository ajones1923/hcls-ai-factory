"""
HPO ontology validator (PRD §3 FR-PM-3 — "every emitted HPO term is verified against
the ontology release, not trusted from the model").

Loads the official Human Phenotype Ontology release (data/ref/hp.obo, ~200k lines) once
into a compact in-memory index, then validates/normalizes the (hpo_id, label) pairs the
Phenome Mapper extracts. This converts a model-guessed code into an ontology-grounded
assertion: the ID is confirmed to exist, the label is replaced with the *official* term,
obsolete/alt IDs are remapped to their primary, and a term named correctly but coded
wrongly is recovered by name lookup.

Falls back gracefully: if the release file is absent (no network at build time), the
validator reports `status="unverified"` and passes terms through unchanged, so the engine
still runs — it simply cannot claim ontology grounding. Build the index once with
`python -m src.utils.hpo` (or it is built lazily on first use).
"""
from __future__ import annotations

import json
import re
import threading
from pathlib import Path

_OBO = Path(__file__).resolve().parents[2] / "data" / "ref" / "hp.obo"
_INDEX = Path(__file__).resolve().parents[2] / "data" / "ref" / "hpo_index.json"
_ID_RE = re.compile(r"^HP:\d{7}$")
_SYN_RE = re.compile(r'^synonym: "(.+?)" (EXACT|NARROW)\b')


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip().lower())


def build_index(obo: Path = _OBO, out: Path = _INDEX) -> dict:
    """Parse hp.obo into {id_to_label, name_to_id, alt_to_primary} and cache as JSON."""
    id_to_label: dict[str, str] = {}
    name_to_id: dict[str, str] = {}
    alt_to_primary: dict[str, str] = {}
    cur: str | None = None
    name: str | None = None
    obsolete = False
    replaced: str | None = None
    syns: list[str] = []

    def flush() -> None:
        nonlocal cur, name, obsolete, replaced, syns
        if cur and name and not obsolete:
            id_to_label[cur] = name
            name_to_id.setdefault(_norm(name), cur)
            for s in syns:
                name_to_id.setdefault(_norm(s), cur)
        if cur and obsolete and replaced:
            alt_to_primary[cur] = replaced
        cur, name, obsolete, replaced, syns = None, None, False, None, []

    for line in obo.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line == "[Term]":
            flush()
            continue
        if line.startswith("id: HP:"):
            cur = line[4:].strip()
        elif line.startswith("name:"):
            name = line[5:].strip()
        elif line.startswith("alt_id: HP:") and cur:
            alt_to_primary[line[8:].strip()] = cur
        elif line.startswith("is_obsolete: true"):
            obsolete = True
        elif line.startswith("replaced_by: HP:"):
            replaced = line[13:].strip()
        elif line.startswith("synonym:"):
            m = _SYN_RE.match(line)
            if m:
                syns.append(m.group(1))
    flush()

    index = {"id_to_label": id_to_label, "name_to_id": name_to_id,
             "alt_to_primary": alt_to_primary, "n_terms": len(id_to_label)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(index), encoding="utf-8")
    return index


class HPOOntology:
    """Lazy-loaded, thread-safe singleton wrapper around the cached HPO index."""

    _lock = threading.Lock()
    _instance: "HPOOntology | None" = None

    def __init__(self) -> None:
        self.available = False
        self._id_to_label: dict[str, str] = {}
        self._name_to_id: dict[str, str] = {}
        self._alt: dict[str, str] = {}
        try:
            if not _INDEX.exists() and _OBO.exists():
                build_index()
            if _INDEX.exists():
                idx = json.loads(_INDEX.read_text(encoding="utf-8"))
                self._id_to_label = idx["id_to_label"]
                self._name_to_id = idx["name_to_id"]
                self._alt = idx["alt_to_primary"]
                self.available = len(self._id_to_label) > 1000
        except Exception:
            self.available = False

    @property
    def n_terms(self) -> int:
        return len(self._id_to_label)

    def label_for(self, hpo_id: str) -> str | None:
        return self._id_to_label.get(hpo_id)

    def resolve(self, hpo_id: str | None, label: str | None) -> dict:
        """Return {hpo_id, label, status} grounding a (possibly model-guessed) pair.

        status: verified | relabeled | remapped | recovered | unknown | unverified
        """
        hid = (hpo_id or "").strip()
        lbl = (label or "").strip()
        if not self.available:
            return {"hpo_id": hid, "label": lbl, "status": "unverified"}

        # follow alt/obsolete -> primary
        primary = self._alt.get(hid, hid)
        if _ID_RE.match(primary) and primary in self._id_to_label:
            official = self._id_to_label[primary]
            remapped = primary != hid
            relabeled = _norm(official) != _norm(lbl)
            status = "remapped" if remapped else ("relabeled" if relabeled else "verified")
            return {"hpo_id": primary, "label": official, "status": status}

        # ID is bad/unknown — try to recover from the label
        by_name = self._name_to_id.get(_norm(lbl))
        if by_name:
            return {"hpo_id": by_name, "label": self._id_to_label[by_name], "status": "recovered"}

        return {"hpo_id": hid, "label": lbl, "status": "unknown"}

    @classmethod
    def get(cls) -> "HPOOntology":
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance


def get_ontology() -> HPOOntology:
    return HPOOntology.get()


if __name__ == "__main__":
    idx = build_index()
    print(f"built HPO index: {idx['n_terms']} terms, "
          f"{len(idx['name_to_id'])} names/synonyms, {len(idx['alt_to_primary'])} alt/obsolete ids")
