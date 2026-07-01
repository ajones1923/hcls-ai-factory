"""
Faithful genomic substrate (PRD §3 FR-CG-3; master paper §15) — the W2 upgrade.

BAMSurgeon spike-in into a control BAM (NA12878 / Genome in a Bottle) -> mosaic-aware
*blinded* variant calling (BWA-MEM + GATK HaplotypeCaller/Mutect2, or NVIDIA Parabricks
on GPU) -> a realistic VCF the TSC-Variant Curator must DISCOVER variants in (rather than
reading them back, as the direct writer in src/cohort/genomic.py does for the demo).

Requires samtools, bcftools, bwa, BAMSurgeon (addsnv.py/addindel.py), and GATK/Parabricks.
Run on a host with these installed — RunPod GPU burst for Parabricks-scale calling
(docker-compose.runpod.yml). `ensure_tools()` fails fast with a clear message otherwise.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from src.cohort.spec import PatientSpec

REQUIRED = ["samtools", "bcftools", "bwa"]   # + bamsurgeon scripts + gatk/pbrun on PATH


def available() -> dict[str, bool]:
    return {t: shutil.which(t) is not None for t in REQUIRED}


def ensure_tools() -> None:
    missing = [t for t, ok in available().items() if not ok]
    if missing:
        raise RuntimeError(
            f"Faithful blinded calling requires {missing} (plus BAMSurgeon and GATK/Parabricks). "
            "Run on a host with these installed — RunPod GPU for Parabricks. "
            "For the demo, use the direct VCF writer (src/cohort/genomic.py)."
        )


def _run(cmd: list[str]) -> None:  # pragma: no cover - requires external tools
    subprocess.run(cmd, check=True)


def spike_and_call(spec: PatientSpec, control_bam: str | Path, out_dir: str | Path,
                   reference: str = "GRCh38.fa") -> Path:  # pragma: no cover - external tools
    """Spike the patient's variant into the control BAM at its VAF, then call -> blinded VCF."""
    ensure_tools()
    if spec.variant is None:
        raise ValueError("NMI patient: no variant to spike (header-only VCF in the demo).")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    v = spec.variant
    is_indel = v["kind"] == "frameshift" or len(v["ref"]) != len(v["alt"])

    # 1) BAMSurgeon variant file: chrom start end VAF [bases]
    varfile = out_dir / "spike.txt"
    end = v["pos"] + max(len(v["ref"]), len(v["alt"]))
    if is_indel:
        op = "DEL" if len(v["ref"]) > len(v["alt"]) else f"INS;{v['alt']}"
        varfile.write_text(f"{v['chrom']} {v['pos']} {end} {spec.vaf} {op}\n")
        spiked = out_dir / "spiked.bam"
        _run(["addindel.py", "-v", str(varfile), "-f", str(control_bam),
              "-r", reference, "-o", str(spiked), "--mindepth", "200"])
    else:
        varfile.write_text(f"{v['chrom']} {v['pos']} {v['pos']} {spec.vaf} {v['alt']}\n")
        spiked = out_dir / "spiked.bam"
        _run(["addsnv.py", "-v", str(varfile), "-f", str(control_bam),
              "-r", reference, "-o", str(spiked), "--mindepth", "200"])

    _run(["samtools", "sort", "-o", str(out_dir / "sorted.bam"), str(spiked)])
    _run(["samtools", "index", str(out_dir / "sorted.bam")])

    # 2) mosaic-aware calling (Mutect2; substitute `pbrun` for Parabricks on GPU)
    vcf = out_dir / "variants.vcf"
    _run(["gatk", "Mutect2", "-R", reference, "-I", str(out_dir / "sorted.bam"),
          "-O", str(vcf)])
    return vcf
