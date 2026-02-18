#!/usr/bin/env python3
"""
Pipeline Performance Benchmark.

Runs the drug discovery pipeline with mocked NIM services and reports
per-stage timing in a formatted table.

Usage:
    python scripts/benchmark.py
    python scripts/benchmark.py --output results.json
"""
import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock

# Add drug-discovery-pipeline to path
sys.path.insert(0, str(Path(__file__).parent.parent / "drug-discovery-pipeline"))

from src.models import PipelineConfig, TargetHypothesis
from src.pipeline import DrugDiscoveryPipeline


def create_mock_nim_manager():
    """Create a mock NIM service manager for benchmarking."""
    manager = MagicMock()
    manager.get_available_services.return_value = ["molmim", "diffdock"]

    # MolMIM returns mock molecules
    manager.molmim.generate.return_value = [
        {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "score": 0.85, "method": "MolMIM"},
        {"smiles": "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", "score": 0.72, "method": "MolMIM"},
        {"smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "score": 0.91, "method": "MolMIM"},
        {"smiles": "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "score": 0.68, "method": "MolMIM"},
        {"smiles": "c1ccc(cc1)C(=O)Nc2ccccc2", "score": 0.77, "method": "MolMIM"},
    ]

    # DiffDock returns mock poses
    manager.diffdock.dock.return_value = [
        {"docking_score": -8.5, "confidence": 0.82, "hydrogen_bonds": 3, "contacts": ["GLU305", "ASN348"]},
        {"docking_score": -7.2, "confidence": 0.65, "hydrogen_bonds": 2, "contacts": ["GLU305"]},
    ]

    return manager


def run_benchmark(output_file=None):
    """Run the benchmark and print results."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        config = PipelineConfig(
            target_gene="VCP",
            num_molecules=5,
            output_dir=tmpdir,
        )

        target = TargetHypothesis(
            gene="VCP",
            protein="Valosin-containing protein (p97)",
            rationale="Benchmark target",
            pdb_ids=["5FTK"],
        )

        mock_nim = create_mock_nim_manager()
        pipeline = DrugDiscoveryPipeline(config, nim_manager=mock_nim)
        run = pipeline.run_pipeline(target)

        # Print results
        timings = pipeline.stage_timings
        total = sum(timings.values())

        print("\n" + "=" * 60)
        print("  HCLS AI Factory - Pipeline Benchmark Results")
        print("=" * 60)
        print(f"  Run ID:     {pipeline.run_id}")
        print(f"  Target:     {config.target_gene}")
        print(f"  Molecules:  {config.num_molecules}")
        print(f"  Status:     {run.status}")
        print("-" * 60)
        print(f"  {'Stage':<5} {'Name':<25} {'Duration':>10}")
        print("-" * 60)

        for stage_num in sorted(timings.keys()):
            name = pipeline.STAGES[stage_num] if stage_num < len(pipeline.STAGES) else f"Stage {stage_num}"
            duration = timings[stage_num]
            pct = (duration / total * 100) if total > 0 else 0
            print(f"  {stage_num:<5} {name:<25} {duration:>8.3f}s ({pct:>4.1f}%)")

        print("-" * 60)
        print(f"  {'TOTAL':<31} {total:>8.3f}s")
        print("=" * 60 + "\n")

        # Build JSON output
        results = {
            "benchmark_date": datetime.now().isoformat(),
            "run_id": pipeline.run_id,
            "config": {
                "target_gene": config.target_gene,
                "num_molecules": config.num_molecules,
            },
            "stage_timings": {str(k): v for k, v in timings.items()},
            "total_seconds": round(total, 3),
            "status": run.status,
        }

        if output_file:
            with open(output_file, "w") as f:
                json.dump(results, f, indent=2)
            print(f"  Results saved to: {output_file}\n")

        return results


def main():
    parser = argparse.ArgumentParser(description="HCLS AI Factory Pipeline Benchmark")
    parser.add_argument("--output", "-o", help="Output JSON file for results")
    args = parser.parse_args()
    run_benchmark(args.output)


if __name__ == "__main__":
    main()
