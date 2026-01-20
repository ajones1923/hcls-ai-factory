# Quick Start Guide

## One-Line Commands

```bash
# Check prerequisites
./run.sh check

# Login to NGC
./run.sh login

# Download data (~200GB, several hours)
./run.sh download

# Setup reference genome
./run.sh reference

# Run fast test (chr20 only, ~5-20 min)
./run.sh test

# Run full genome (~120-240 min)
./run.sh full
```

## Complete Workflow (Copy-Paste)

```bash
# Navigate to project
cd /home/adam/transfer/genomics-pipeline

# Run all setup steps
./run.sh check && \
./run.sh login && \
./run.sh download && \
./run.sh reference && \
./run.sh test

# If test succeeds, run full genome
./run.sh full
```

## Minimal Test (Skip Large Download)

If you want to test the pipeline without downloading 200GB:

```bash
# 1. Check prerequisites
./run.sh check

# 2. Login to NGC
./run.sh login

# 3. Setup reference only
./run.sh reference

# 4. Use your own FASTQ files
# Copy them to: data/input/HG002_R1.fastq.gz
#               data/input/HG002_R2.fastq.gz

# 5. Run pipeline
./run.sh test   # or ./run.sh full
```

## Monitor Progress

```bash
# Watch GPU utilization
watch -n 1 nvidia-smi

# Follow logs in real-time
tail -f data/output/logs/genome_fq2bam.log

# Check output files
ls -lh data/output/
```

## Expected Output

After successful completion:

```
data/output/
├── HG002.chr20.bam            # Chr20 test BAM
├── HG002.chr20.vcf.gz         # Chr20 test VCF
├── HG002.genome.bam           # Full genome BAM
├── HG002.genome.vcf.gz        # Full genome VCF ← Main output!
├── HG002.genome.vcf.gz.tbi    # VCF index
└── logs/                      # All execution logs
```

## Troubleshooting

**GPU out of memory?**
```bash
# Edit config/pipeline.env
# Set: LOW_MEMORY=1
./run.sh full
```

**Download interrupted?**
```bash
# Re-run download script - aria2 will resume
./run.sh download
```

**Start over?**
```bash
# Clean outputs only (keeps downloaded data)
./run.sh clean

# Clean everything
./run.sh clean-all
```

## Time Estimates

| Operation | Duration |
|-----------|----------|
| Prerequisites check | < 1 min |
| NGC login | < 1 min |
| Data download | 2-6 hours |
| Reference setup | 5-15 min |
| Chr20 test | 5-20 min |
| Full genome | 120-240 min |

## Next Steps

Use the output VCF for:
- Pharmacogenomics analysis (Procedure 3)
- Variant annotation
- Clinical interpretation
- Research workflows

---

**Need help?** See [README.md](README.md) for detailed documentation.
