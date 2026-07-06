"""
Queryable variant store (E1 Genomic Foundation) — load VCF into DuckDB for cohort/agent queries without a
warehouse. Columnar, SQL-queryable by region/filter/sample, with a real Ts/Tv QC metric.

E1 F1 adds read-backed **VAF / somatic-mosaicism** signal: per-call allelic depths (AD/DP)
are parsed from the FORMAT/sample columns and a low-VAF mosaic-*candidate* flag surfaces the
diagnostic trap standard germline pipelines drop (mosaic TSC1/TSC2). This is *evidence for a
clinician*, never a classification.
"""
import gzip
import re
from pathlib import Path

_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

# QC expected ranges (whole-genome, PASS biallelic) — inline in the report so a caller sees the bar
_QC_EXPECT = {"ts_tv": "~2.0-2.1", "het_hom": "~1.3-2.0"}
_VERDICT_RANK = {"fail": 0, "warn": 1, "pass": 2}


class VariantStore:
    def __init__(self, db_path: str = ":memory:"):
        import duckdb
        self.con = duckdb.connect(db_path)
        self.con.execute(
            "CREATE TABLE IF NOT EXISTS variants("
            "chrom VARCHAR, pos BIGINT, id VARCHAR, ref VARCHAR, alt VARCHAR,"
            "qual DOUBLE, filter VARCHAR, sample VARCHAR,"
            # E1 F1: read-backed genotype + allelic depths + VAF (null for FORMAT-less VCFs)
            "gt VARCHAR, dp BIGINT, ad_ref BIGINT, ad_alt BIGINT, vaf DOUBLE)")

    @staticmethod
    def _parse_format(fmt: str, sample_field: str):
        """Parse a VCF FORMAT + first sample column -> (gt, dp, ad_ref, ad_alt, vaf).

        VAF = ad_alt / dp, and is None when DP or the ALT allelic depth is absent. Multiallelic
        ALT uses the first ALT's AD index (a documented limitation; SV/multiallelic is P2). Any
        malformed field degrades to None rather than raising — a FORMAT-less VCF still loads.
        """
        if not fmt or not sample_field:
            return (None, None, None, None, None)
        d = dict(zip(fmt.split(":"), sample_field.split(":")))

        def _int(x):
            try:
                return int(x) if x not in (None, ".", "") else None
            except (ValueError, TypeError):
                return None

        gt = d.get("GT") or None
        dp = _int(d.get("DP"))
        ad_ref = ad_alt = None
        ad = d.get("AD")
        if ad and ad not in (".", ""):
            parts = ad.split(",")
            ad_ref = _int(parts[0]) if parts else None
            ad_alt = _int(parts[1]) if len(parts) > 1 else None
        vaf = round(ad_alt / dp, 6) if (ad_alt is not None and dp and dp > 0) else None
        return (gt, dp, ad_ref, ad_alt, vaf)

    def load_vcf(self, path, sample: str = "HG002", limit: int | None = None) -> int:
        opener = gzip.open if str(path).endswith(".gz") else open
        batch, n = [], 0
        cols = "(?,?,?,?,?,?,?,?,?,?,?,?,?)"   # 13 columns (8 core + gt/dp/ad_ref/ad_alt/vaf)
        with opener(path, "rt") as f:
            for line in f:
                if not line or line[0] == "#":
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 7:
                    continue
                try:
                    qual = None if p[5] in (".", "") else float(p[5])
                except ValueError:
                    qual = None
                gt = dp = ad_ref = ad_alt = vaf = None
                if len(p) >= 10:                       # FORMAT (p[8]) + first sample (p[9])
                    gt, dp, ad_ref, ad_alt, vaf = self._parse_format(p[8], p[9])
                batch.append((p[0], int(p[1]), p[2], p[3], p[4], qual, p[6], sample,
                              gt, dp, ad_ref, ad_alt, vaf))
                n += 1
                if len(batch) >= 20000:
                    self.con.executemany(f"INSERT INTO variants VALUES {cols}", batch); batch = []
                if limit and n >= limit:
                    break
        if batch:
            self.con.executemany(f"INSERT INTO variants VALUES {cols}", batch)
        return n

    def count(self) -> int:
        return self.con.execute("SELECT count(*) FROM variants").fetchone()[0]

    def query_region(self, chrom: str, start: int, end: int) -> list[dict]:
        cur = self.con.execute(
            "SELECT chrom,pos,ref,alt,qual,filter FROM variants "
            "WHERE chrom=? AND pos BETWEEN ? AND ? ORDER BY pos", (chrom, start, end))
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    def n_pass(self) -> int:
        return self.con.execute("SELECT count(*) FROM variants WHERE filter='PASS'").fetchone()[0]

    def pass_rate(self) -> float:
        total = self.count() or 1
        return round(self.n_pass() / total, 4)

    def ts_tv(self) -> float:
        """Transition/transversion ratio over PASS biallelic SNVs (a real genomics QC metric).
        Computed over PASS calls only -- gVCF reference (RefCall) sites are excluded."""
        rows = self.con.execute(
            "SELECT ref,alt FROM variants WHERE filter='PASS' AND length(ref)=1 AND length(alt)=1 "
            "AND ref IN ('A','C','G','T') AND alt IN ('A','C','G','T')").fetchall()
        ts = sum(1 for r, a in rows if (r, a) in _TRANSITIONS)
        tv = len(rows) - ts
        return round(ts / tv, 3) if tv else 0.0

    # ---- E1 F1: VAF / somatic-mosaicism (evidence, never a classification) ---- #
    def mosaic_candidates(self, vaf_lo: float = 0.02, vaf_hi: float = 0.35,
                          min_dp: int = 30) -> list[dict]:
        """PASS calls with ``vaf`` in ``[vaf_lo, vaf_hi]`` and ``dp >= min_dp`` — low-VAF mosaic
        *candidates* the Rare Disease agent can act on (the mosaic TSC1/TSC2 trap standard pipelines
        filter as noise). Never drops a PASS call; carries the FILTER through. Empty (no rows) when
        the source VCF had no AD/DP — see ``stats()['mosaicism']['reason']``.
        """
        cur = self.con.execute(
            "SELECT chrom,pos,ref,alt,filter,dp,ad_ref,ad_alt,vaf FROM variants "
            "WHERE filter='PASS' AND vaf IS NOT NULL AND vaf BETWEEN ? AND ? AND dp >= ? "
            "ORDER BY vaf", (vaf_lo, vaf_hi, min_dp))
        cols = [c[0] for c in cur.description]
        out = []
        for r in cur.fetchall():
            row = dict(zip(cols, r))
            # INFO strand-bias / low-complexity flags are a P2 refinement; empty for now so a
            # caller never mistakes absence-of-flags for a validated mosaic call.
            row["flags"] = []
            out.append(row)
        return out

    def _mosaicism_stats(self) -> dict:
        cand = self.mosaic_candidates()
        buckets = {"0.02-0.10": 0, "0.10-0.20": 0, "0.20-0.35": 0}
        dps = []
        for c in cand:
            v = c["vaf"]
            dps.append(c["dp"])
            if v < 0.10:
                buckets["0.02-0.10"] += 1
            elif v < 0.20:
                buckets["0.10-0.20"] += 1
            else:
                buckets["0.20-0.35"] += 1
        median_dp = None
        if dps:
            dps.sort()
            m = len(dps) // 2
            median_dp = dps[m] if len(dps) % 2 else (dps[m - 1] + dps[m]) / 2
        has_vaf = self.con.execute(
            "SELECT count(*) FROM variants WHERE vaf IS NOT NULL").fetchone()[0] > 0
        return {"n_candidates": len(cand), "vaf_buckets": buckets, "median_dp": median_dp,
                "reason": None if has_vaf else "no AD/DP in source VCF"}

    def stats(self) -> dict:
        return {"n_variants": self.count(), "n_pass": self.n_pass(),
                "pass_rate": self.pass_rate(), "ts_tv": self.ts_tv(),
                "mosaicism": self._mosaicism_stats()}

    # ---- E1 F3: QC trust-gate (refuses to advertise a bad call set) ---------- #
    def _het_hom(self) -> float | None:
        """Het / hom-alt ratio over PASS calls with a parsed GT (needs F1). Hom-ref (0/0) excluded."""
        rows = self.con.execute(
            "SELECT gt FROM variants WHERE filter='PASS' AND gt IS NOT NULL").fetchall()
        het = hom = 0
        for (gt,) in rows:
            alleles = [a for a in re.split(r"[/|]", gt) if a not in (".", "")]
            if len(alleles) != 2:
                continue
            a, b = alleles
            if a != b:
                het += 1
            elif a != "0":            # hom-alt (exclude hom-ref 0/0)
                hom += 1
        return round(het / hom, 3) if hom else None

    def _snv_indel(self) -> float | None:
        snv = self.con.execute(
            "SELECT count(*) FROM variants WHERE filter='PASS' AND length(ref)=1 AND length(alt)=1"
        ).fetchone()[0]
        indel = self.con.execute(
            "SELECT count(*) FROM variants WHERE filter='PASS' AND (length(ref)>1 OR length(alt)>1)"
        ).fetchone()[0]
        return round(snv / indel, 3) if indel else None

    def qc_report(self) -> dict:
        """A QC *verdict*, not just numbers: Ts/Tv (PASS biallelic), het/hom, SNV/indel, and a
        ``verdict in {pass, warn, fail}`` with the expected ranges inline. A false-positive-heavy
        call set (Ts/Tv far from ~2.0) is flagged BEFORE any interpretation. ``contamination`` is a
        hook — null until a real estimator is wired (honestly labeled)."""
        tstv = self.ts_tv()
        hethom = self._het_hom()
        flags = []
        if not (1.8 <= tstv <= 2.3):
            flags.append({"metric": "ts_tv", "value": tstv, "note": f"outside {_QC_EXPECT['ts_tv']}"})
        if hethom is not None and not (1.3 <= hethom <= 2.5):
            flags.append({"metric": "het_hom", "value": hethom, "note": f"outside {_QC_EXPECT['het_hom']}"})
        hard_fail = tstv < 1.6 or tstv > 2.6            # FP-heavy / broken call set
        verdict = "fail" if hard_fail else ("warn" if flags else "pass")
        return {"ts_tv_pass_biallelic": tstv, "het_hom": hethom, "snv_indel": self._snv_indel(),
                "n_pass": self.n_pass(), "contamination": None, "expected": dict(_QC_EXPECT),
                "verdict": verdict, "flags": flags}

    def qc_gate(self, min_ok: str = "warn") -> bool:
        """Trust-gate: True iff the call set is trustworthy enough to interpret. ``min_ok='warn'``
        allows pass|warn; ``'pass'`` requires pass. On a ``fail`` set, a caller must WITHHOLD
        interpretation and surface the red flag instead of the variants (the executable honesty
        guardrail: an engine that says 'I don't trust this data')."""
        return _VERDICT_RANK[self.qc_report()["verdict"]] >= _VERDICT_RANK[min_ok]
