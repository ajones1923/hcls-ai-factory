"""
Queryable variant store (E2) — load VCF into DuckDB for cohort/agent queries without a
warehouse. Columnar, SQL-queryable by region/filter/sample, with a real Ts/Tv QC metric.
"""
import gzip
from pathlib import Path

_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


class VariantStore:
    def __init__(self, db_path: str = ":memory:"):
        import duckdb
        self.con = duckdb.connect(db_path)
        self.con.execute(
            "CREATE TABLE IF NOT EXISTS variants("
            "chrom VARCHAR, pos BIGINT, id VARCHAR, ref VARCHAR, alt VARCHAR,"
            "qual DOUBLE, filter VARCHAR, sample VARCHAR)")

    def load_vcf(self, path, sample: str = "HG002", limit: int | None = None) -> int:
        opener = gzip.open if str(path).endswith(".gz") else open
        batch, n = [], 0
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
                batch.append((p[0], int(p[1]), p[2], p[3], p[4], qual, p[6], sample))
                n += 1
                if len(batch) >= 20000:
                    self.con.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?)", batch); batch = []
                if limit and n >= limit:
                    break
        if batch:
            self.con.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?,?,?)", batch)
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

    def stats(self) -> dict:
        return {"n_variants": self.count(), "n_pass": self.n_pass(),
                "pass_rate": self.pass_rate(), "ts_tv": self.ts_tv()}
