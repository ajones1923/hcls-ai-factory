"""
Lightweight single-box GWAS — association testing (E5).

Per-variant logistic-regression association of genotypes against a binary case/control
phenotype, with optional covariates — pure Python (statsmodels + scipy), no Spark/JVM, ARM-clean.
This is the single-box-appropriate substitute for a cluster-scale statistical-genomics stack:
the DuckDB variant store (E3) already provides columnar variant storage at single-box scale, and
this adds the association test on top. Separation/convergence failures are handled gracefully
(reported as untestable, never a crash).
"""
from __future__ import annotations

import math
from typing import Any, Sequence


def single_variant_test(genotypes: Sequence[float], phenotype: Sequence[int],
                        covariates: Sequence[Sequence[float]] | None = None) -> dict[str, Any]:
    """Logistic regression of one variant's genotypes (0/1/2) vs a binary phenotype."""
    import numpy as np
    import statsmodels.api as sm

    g = np.asarray(genotypes, dtype=float)
    y = np.asarray(phenotype, dtype=float)
    cols = [np.ones(len(g)), g]
    if covariates is not None:
        cov = np.asarray(covariates, dtype=float)
        if cov.ndim == 1:
            cov = cov.reshape(-1, 1)
        cols.extend(cov.T)
    X = np.column_stack(cols)
    try:
        res = sm.Logit(y, X).fit(disp=0, maxiter=100)
        beta = float(res.params[1])
        p = float(res.pvalues[1])
        if math.isnan(p):
            return {"tested": False, "reason": "non-finite p-value"}
        return {"tested": True, "beta": round(beta, 5), "p_value": p,
                "odds_ratio": round(math.exp(beta), 5)}
    except Exception as e:  # PerfectSeparation, singular matrix, non-convergence
        return {"tested": False, "reason": type(e).__name__}


def run_gwas(genotype_matrix: Sequence[Sequence[float]], phenotype: Sequence[int],
             variant_ids: Sequence[str], covariates: Sequence[Sequence[float]] | None = None,
             min_maf: float = 0.0) -> dict[str, Any]:
    """Association-test each variant (rows of genotype_matrix). Returns results sorted by p-value."""
    import numpy as np
    results = []
    skipped = 0
    for vid, geno in zip(variant_ids, genotype_matrix):
        g = np.asarray(geno, dtype=float)
        maf = min(g.mean() / 2.0, 1.0 - g.mean() / 2.0)     # allele freq from 0/1/2 dosage
        if maf < min_maf:
            skipped += 1
            continue
        r = single_variant_test(geno, phenotype, covariates)
        results.append({"variant_id": vid, "maf": round(maf, 4), **r})
    tested = [r for r in results if r.get("tested")]
    tested.sort(key=lambda r: r["p_value"])
    return {
        "n_variants": len(variant_ids),
        "n_tested": len(tested),
        "n_skipped_maf": skipped,
        "n_untestable": len(results) - len(tested),
        "results": results,
        "top": tested[:25],
        "manhattan": [{"variant_id": r["variant_id"], "neg_log10_p": round(-math.log10(r["p_value"]), 3)}
                      for r in tested if r["p_value"] > 0],
    }
