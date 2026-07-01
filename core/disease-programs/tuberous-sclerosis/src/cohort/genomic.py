"""
Genomic substrate (PRD §3 FR-CG-2/3; master paper §15).

Emits a VCFv4.3-conformant per-patient VCF carrying the TSC1/TSC2 variant at the specified
allele fraction (germline ~0.5, mosaic 0.04-0.12) with **read-level evidence** — strand-
resolved allelic depths (ADF/ADR), a strand-balance (SB), and depth — so the Variant
Curator's discrimination is real: a true variant is strand-balanced, and every sample also
carries a low-VAF, strongly strand-biased candidate (FILTER=strand_bias) that a naive
caller might surface but the curator MUST reject. That makes "0 false-positive Pathogenic"
a real result, not a freebie.

BUILD NOTE: the faithful path is BAMSurgeon spike-in into a control BAM -> Parabricks /
BWA-MEM + GATK HaplotypeCaller/Mutect2 (mosaic-aware) -> *blinded* VCF the agent must
discover variants in. That requires samtools/bcftools/BAMSurgeon + GPU (RunPod burst) and
is the W2 upgrade; this writer synthesizes realistic read-level fields deterministically so
the rest of the engine — and its read-level artifact rejection — runs on the Spark today.
"""
from __future__ import annotations

import random
from pathlib import Path

from src.cohort.spec import PatientSpec

VCF_HEADER = """##fileformat=VCFv4.3
##source=TSC-Intelligence-Engine-synthetic-cohort (SYNTHETIC; not real patient data)
##reference=GRCh38
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=strand_bias,Description="Alt reads strongly strand-imbalanced; likely artifact">
##INFO=<ID=GENE,Number=1,Type=String,Description="Affected gene">
##INFO=<ID=HGVSc,Number=1,Type=String,Description="cDNA change">
##INFO=<ID=HGVSp,Number=1,Type=String,Description="Protein change">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence">
##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description="Synthetic demonstration variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">
##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Forward-strand allelic depths (ref,alt)">
##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reverse-strand allelic depths (ref,alt)">
##FORMAT=<ID=SB,Number=1,Type=Float,Description="Alt-allele strand balance (0.5=balanced)">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

_FMT = "GT:DP:AD:ADF:ADR:SB:AF"


def _sample_field(gt: str, dp: int, alt: int, sb: float) -> str:
    """Build the FORMAT sample with strand-resolved allelic depths from a strand balance."""
    ref = dp - alt
    alt_f = round(alt * sb)
    alt_r = alt - alt_f
    ref_f = ref // 2
    ref_r = ref - ref_f
    af = alt / dp if dp else 0.0
    return f"{gt}:{dp}:{ref},{alt}:{ref_f},{alt_f}:{ref_r},{alt_r}:{sb:.3f}:{af:.3f}"


def write_vcf(spec: PatientSpec, out_path: str | Path, seed: int = 0) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(f"{seed}:{spec.patient_id}:genomic")
    lines = [VCF_HEADER.format(sample=spec.patient_id)]

    if spec.gene is not None and spec.variant is not None:
        v = spec.variant
        # high coverage so low-VAF mosaic is detectable (tissue ~250x; germline 120x)
        dp = 250 if spec.zygosity == "mosaic" else 120
        alt = max(1, round(spec.vaf * dp))
        qual = 200 if spec.zygosity == "germline" else 60  # lower QUAL for low-VAF
        sb = round(0.5 + rng.uniform(-0.06, 0.06), 3)       # true variant: strand-balanced
        info = (f"GENE={spec.gene};HGVSc={v['cdna']};HGVSp={v['protein']};"
                f"CONSEQUENCE={v['kind']};SYNTHETIC")
        lines.append(
            f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}\t{qual}\tPASS\t{info}\t"
            f"{_FMT}\t{_sample_field('0/1', dp, alt, sb)}\n"
        )

    # Strand-biased low-VAF artifact every sample carries — must be rejected, never reported.
    art_dp = 240
    art_alt = max(2, round(0.04 * art_dp))
    art_sb = 0.94 if rng.random() < 0.5 else 0.05
    lines.append(
        f"chr16\t2090117\t.\tG\tA\t38\tstrand_bias\tSYNTHETIC\t"
        f"{_FMT}\t{_sample_field('0/1', art_dp, art_alt, art_sb)}\n"
    )
    # NMI patients thus get a header + only the filtered artifact (no PASS variant) — the gap.

    out_path.write_text("".join(lines))
    return out_path


def _strand_balance(adf: tuple[int, int], adr: tuple[int, int]) -> float:
    alt_f, alt_r = adf[1], adr[1]
    tot = alt_f + alt_r
    return round(alt_f / tot, 3) if tot else 0.5


def parse_variants(vcf_path: str | Path, include_filtered: bool = False) -> list[dict]:
    """Minimal VCF reader (no pysam dependency) returning variant records with read-level
    fields. By default returns only PASS records; set include_filtered to inspect artifacts."""
    records = []
    for line in Path(vcf_path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        cols = line.split("\t")
        filt = cols[6]
        if not include_filtered and filt not in ("PASS", "."):
            continue
        info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv)
        fmt, sample = cols[8].split(":"), cols[9].split(":")
        s = dict(zip(fmt, sample))

        def _pair(key):
            if key in s and "," in s[key]:
                a, b = s[key].split(",")[:2]
                return (int(a), int(b))
            return (0, 0)

        ad, adf, adr = _pair("AD"), _pair("ADF"), _pair("ADR")
        sb = float(s["SB"]) if "SB" in s else _strand_balance(adf, adr)
        records.append({
            "chrom": cols[0], "pos": int(cols[1]), "ref": cols[3], "alt": cols[4],
            "qual": float(cols[5]) if cols[5] != "." else None, "filter": filt,
            "gene": info.get("GENE"), "hgvsc": info.get("HGVSc"), "hgvsp": info.get("HGVSp"),
            "consequence": info.get("CONSEQUENCE"),
            "af": float(s.get("AF", 0)), "dp": int(s.get("DP", 0)),
            "alt_reads": ad[1], "ref_reads": ad[0], "strand_balance": sb,
        })
    return records


def count_filtered_artifacts(vcf_path: str | Path) -> int:
    """How many non-PASS candidate records the caller emitted (artifacts the curator rejects)."""
    all_recs = parse_variants(vcf_path, include_filtered=True)
    pass_recs = parse_variants(vcf_path, include_filtered=False)
    return len(all_recs) - len(pass_recs)
