"""Tests for the variant store (E2). Real DuckDB, tiny synthetic VCF."""
import gzip, os, tempfile
from variant_store import VariantStore

VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t50\tPASS\t.
chr1\t200\t.\tC\tT\t50\tPASS\t.
chr1\t300\t.\tA\tC\t50\tPASS\t.
chr1\t400\t.\tG\tT\t40\tRefCall\t.
chr1\t500\t.\tT\tA\t30\tPASS\t.
"""

def _vcf_file():
    fd, path = tempfile.mkstemp(suffix=".vcf"); os.write(fd, VCF.encode()); os.close(fd); return path

class TestStore:
    def test_load_and_count(self):
        vs = VariantStore(); n = vs.load_vcf(_vcf_file(), sample="X")
        assert n == 5 and vs.count() == 5

    def test_n_pass_excludes_refcall(self):
        vs = VariantStore(); vs.load_vcf(_vcf_file())
        assert vs.n_pass() == 4                       # the RefCall is excluded

    def test_ts_tv_over_pass_only(self):
        vs = VariantStore(); vs.load_vcf(_vcf_file())
        # PASS SNVs: A>G(ts), C>T(ts), A>C(tv), T>A(tv) -> ts=2 tv=2 -> 1.0
        assert vs.ts_tv() == 1.0

    def test_region_query(self):
        vs = VariantStore(); vs.load_vcf(_vcf_file())
        hits = vs.query_region("chr1", 150, 450)
        assert [h["pos"] for h in hits] == [200, 300, 400]

    def test_stats_shape(self):
        vs = VariantStore(); vs.load_vcf(_vcf_file())
        s = vs.stats()
        assert s["n_variants"] == 5 and s["n_pass"] == 4 and "ts_tv" in s


# ── E1 F1: VAF / somatic-mosaicism flag ──────────────────────────────────────
import os as _os  # noqa: E402

_MINI = _os.path.join(_os.path.dirname(__file__), "fixtures", "mosaic_mini.vcf")


def test_vaf_parsed():
    v = VariantStore()
    v.load_vcf(_MINI)
    rows = {r[0]: r[1] for r in v.con.execute("SELECT pos, vaf FROM variants ORDER BY pos").fetchall()}
    assert rows[2098000] == round(22 / 420, 6)      # low-VAF PASS call
    assert rows[2099000] == round(205 / 415, 6)     # het PASS call
    assert rows[2100000] is None                    # GT-only row -> no VAF


def test_mosaic_candidates_bounds():
    v = VariantStore()
    v.load_vcf(_MINI)
    cand = v.mosaic_candidates()                    # defaults 0.02-0.35, dp>=30
    # only the ~0.05 PASS row qualifies
    assert [c["pos"] for c in cand] == [2098000]
    c = cand[0]
    assert c["filter"] == "PASS" and c["dp"] == 420 and c["ad_alt"] == 22
    assert 0.02 <= c["vaf"] <= 0.35 and c["flags"] == []
    # the ~0.49 het row is excluded (VAF too high); the RefCall row is excluded (not PASS)
    positions = {c["pos"] for c in cand}
    assert 2099000 not in positions and 2101000 not in positions


def test_vaf_null_when_no_ad():
    # a FORMAT-less (8-column) VCF still loads; VAF is null and /mosaic is empty with a reason
    import tempfile
    body = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\nchr1\t200\t.\tC\tT\t50\tPASS\t.\n")
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as f:
        f.write(body); path = f.name
    v = VariantStore()
    n = v.load_vcf(path)
    assert n == 2                                   # still loads
    assert v.con.execute("SELECT count(*) FROM variants WHERE vaf IS NOT NULL").fetchone()[0] == 0
    assert v.mosaic_candidates() == []
    assert v.stats()["mosaicism"]["reason"] == "no AD/DP in source VCF"


# ── E1 F3: QC trust-gate ─────────────────────────────────────────────────────
_GOOD_QC = _os.path.join(_os.path.dirname(__file__), "fixtures", "good_qc.vcf")
_BAD_QC = _os.path.join(_os.path.dirname(__file__), "fixtures", "bad_qc.vcf")


def test_qc_report_ranges():
    # a healthy set (Ts/Tv ~2.0, het/hom in range) -> pass, and the gate lets it through
    v = VariantStore()
    v.load_vcf(_GOOD_QC)
    r = v.qc_report()
    assert r["ts_tv_pass_biallelic"] == 2.0
    assert 1.3 <= r["het_hom"] <= 2.5
    assert r["contamination"] is None                 # hook, honestly null (no estimator yet)
    assert r["expected"]["ts_tv"] == "~2.0-2.1"
    assert r["verdict"] == "pass" and r["flags"] == []
    assert v.qc_gate() is True and v.qc_gate("pass") is True


def test_qc_gate_fails_bad_set():
    # a transversion-heavy (FP-looking) set -> fail; the gate WITHHOLDS interpretation
    v = VariantStore()
    v.load_vcf(_BAD_QC)
    r = v.qc_report()
    assert r["ts_tv_pass_biallelic"] < 1.6
    assert r["verdict"] == "fail"
    assert any(f["metric"] == "ts_tv" for f in r["flags"])
    assert v.qc_gate() is False and v.qc_gate("pass") is False
