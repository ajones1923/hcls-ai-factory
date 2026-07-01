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
