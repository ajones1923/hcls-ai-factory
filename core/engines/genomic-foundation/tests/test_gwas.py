"""E5: lightweight single-box GWAS (logistic-regression association)."""
import numpy as np
from gwas import single_variant_test, run_gwas


def _synthetic(n=300, seed=0):
    rng = np.random.default_rng(seed)
    causal = rng.integers(0, 3, n).astype(float)
    prob = 1.0 / (1.0 + np.exp(-(-1.0 + 1.4 * causal)))     # phenotype driven by 'causal'
    pheno = (rng.random(n) < prob).astype(int)
    null = rng.integers(0, 3, n).astype(float)              # unrelated variant
    return causal, null, pheno


class TestGWAS:
    def test_detects_causal_variant(self):
        causal, null, pheno = _synthetic()
        r = single_variant_test(causal, pheno)
        assert r["tested"] and r["p_value"] < 0.001 and r["odds_ratio"] > 1.0

    def test_null_variant_not_significant(self):
        causal, null, pheno = _synthetic()
        r = single_variant_test(null, pheno)
        assert r["tested"] and r["p_value"] > 0.05

    def test_run_gwas_ranks_causal_first(self):
        causal, null, pheno = _synthetic()
        out = run_gwas([causal, null], pheno, ["rs_causal", "rs_null"])
        assert out["n_tested"] == 2
        assert out["top"][0]["variant_id"] == "rs_causal"   # lowest p-value first
        assert len(out["manhattan"]) == 2

    def test_separation_handled_gracefully(self):
        # perfectly separating genotype -> untestable, not a crash
        pheno = [0] * 10 + [1] * 10
        geno = [0.0] * 10 + [2.0] * 10
        r = single_variant_test(geno, pheno)
        assert "tested" in r   # returns a verdict either way, never raises
