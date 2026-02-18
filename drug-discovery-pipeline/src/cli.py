"""
Drug Discovery Pipeline CLI.

Typer-based command line interface for the ddpipe tool.
Based on phase-5-6.pdf specification.

Usage:
    ddpipe run --target VCP --config config.yaml
    ddpipe validate --input target.json
    ddpipe services --check
    ddpipe stage 4 --run-id abc123
    ddpipe report --run-id abc123
    ddpipe export --run-id abc123 --format sdf
"""
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import typer
from loguru import logger
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import PipelineConfig, TargetHypothesis
from src.nim_clients import create_nim_clients
from src.pipeline import DrugDiscoveryPipeline, run_vcp_demo_pipeline

app = typer.Typer(
    name="ddpipe",
    help="Drug Discovery Pipeline CLI - From target to molecule candidates",
    add_completion=False,
)
console = Console()


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
    )


@app.command()
def run(
    target: str = typer.Option(..., "--target", "-t", help="Target gene symbol (e.g., VCP)"),  # noqa: B008
    config_file: Path | None = typer.Option(None, "--config", "-c", help="Config YAML file"),  # noqa: B008
    output_dir: str = typer.Option("outputs", "--output", "-o", help="Output directory"),  # noqa: B008
    num_molecules: int = typer.Option(20, "--num-molecules", "-n", help="Number of molecules to generate"),  # noqa: B008
    seed_smiles: str | None = typer.Option(None, "--seed", "-s", help="Seed SMILES for generation"),  # noqa: B008
    use_mock: bool = typer.Option(False, "--mock", help="Use mock NIM services"),  # noqa: B008
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Verbose output"),  # noqa: B008
):
    """
    Run the drug discovery pipeline for a target.

    Example:
        ddpipe run --target VCP --num-molecules 50
    """
    setup_logging(verbose)

    console.print(Panel.fit(
        f"[bold blue]Drug Discovery Pipeline[/bold blue]\n"
        f"Target: [green]{target}[/green]\n"
        f"Molecules: [cyan]{num_molecules}[/cyan]",
        title="Starting Pipeline"
    ))

    # Build config
    config = PipelineConfig(
        target_gene=target.upper(),
        reference_compound_smiles=seed_smiles,
        num_molecules=num_molecules,
        output_dir=output_dir,
    )

    # Create pipeline
    nim = create_nim_clients(use_mock=use_mock)
    pipeline = DrugDiscoveryPipeline(
        config=config,
        nim_manager=nim,
        progress_callback=lambda stage, msg: console.print(f"  [{stage}] {msg}"),
    )

    # Run with progress
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Running pipeline...", total=10)

        try:
            result = pipeline.run_pipeline()

            for _i in range(10):
                progress.update(task, advance=1)

        except Exception as e:
            console.print(f"[red]Pipeline failed: {e}[/red]")
            raise typer.Exit(1) from None

    # Show results
    console.print("\n[bold green]Pipeline Complete![/bold green]\n")

    # Display top candidates
    if pipeline.ranked_candidates:
        table = Table(title=f"Top {len(pipeline.ranked_candidates)} Candidates")
        table.add_column("Rank", style="cyan")
        table.add_column("ID", style="dim")
        table.add_column("Score", style="green")
        table.add_column("MW")
        table.add_column("Docking")

        for c in pipeline.ranked_candidates:
            mw = f"{c.properties.molecular_weight:.0f}" if c.properties else "-"
            table.add_row(
                str(c.rank),
                c.molecule_id,
                f"{c.composite_score:.3f}",
                mw,
                f"{c.docking_score:.1f}",
            )

        console.print(table)

    console.print(f"\nResults saved to: [cyan]{output_dir}[/cyan]")


@app.command()
def validate(
    input_file: Path = typer.Argument(..., help="Target hypothesis JSON file"),  # noqa: B008
    verbose: bool = typer.Option(False, "--verbose", "-v"),  # noqa: B008
):
    """
    Validate a target hypothesis file.

    Example:
        ddpipe validate target.json
    """
    setup_logging(verbose)

    if not input_file.exists():
        console.print(f"[red]File not found: {input_file}[/red]")
        raise typer.Exit(1)

    try:
        with open(input_file) as f:
            data = json.load(f)

        target = TargetHypothesis(**data)
        console.print("[green]✓ Target hypothesis is valid[/green]")

        # Show summary
        console.print(f"\n  Gene: [cyan]{target.gene}[/cyan]")
        console.print(f"  Protein: {target.protein or 'Not specified'}")
        console.print(f"  Priority: {target.priority}")
        console.print(f"  Structures: {len(target.pdb_ids)}")
        console.print(f"  Status: {target.status}")

    except Exception as e:
        console.print(f"[red]✗ Validation failed: {e}[/red]")
        raise typer.Exit(1) from None


@app.command()
def services(
    check: bool = typer.Option(True, "--check", help="Check service availability"),  # noqa: B008
):
    """
    Check NIM service status.

    Example:
        ddpipe services --check
    """
    console.print("[bold]NIM Service Status[/bold]\n")

    nim = create_nim_clients(use_mock=False)

    table = Table()
    table.add_column("Service", style="cyan")
    table.add_column("Status")
    table.add_column("URL")

    # Check MolMIM
    molmim_ok = nim.molmim.check_health()
    molmim_status = "[green]✓ Available[/green]" if molmim_ok else "[red]✗ Unavailable[/red]"
    table.add_row("MolMIM", molmim_status, nim.molmim_config.base_url)

    # Check DiffDock
    diffdock_ok = nim.diffdock.check_health()
    diffdock_status = "[green]✓ Available[/green]" if diffdock_ok else "[red]✗ Unavailable[/red]"
    table.add_row("DiffDock", diffdock_status, nim.diffdock_config.base_url)

    console.print(table)

    if not molmim_ok and not diffdock_ok:
        console.print("\n[yellow]Note: Mock services will be used when running pipeline[/yellow]")


@app.command()
def stage(
    stage_num: int = typer.Argument(..., help="Stage number (0-9)"),  # noqa: B008
    run_id: str = typer.Option(..., "--run-id", "-r", help="Pipeline run ID"),  # noqa: B008
    output_dir: str = typer.Option("outputs", "--output", "-o"),  # noqa: B008
):
    """
    Run a specific pipeline stage.

    Example:
        ddpipe stage 4 --run-id abc123
    """
    stage_names = [
        "Initialize",
        "Normalize Target",
        "Structure Discovery",
        "Structure Prep",
        "Molecule Generation",
        "Chemistry QC",
        "Conformers",
        "Docking",
        "Ranking",
        "Reporting",
    ]

    if stage_num < 0 or stage_num > 9:
        console.print("[red]Stage must be 0-9[/red]")
        raise typer.Exit(1)

    console.print(f"Running stage {stage_num}: [cyan]{stage_names[stage_num]}[/cyan]")
    console.print(f"Run ID: {run_id}")

    # Would load state and run specific stage
    console.print("[yellow]Note: Single stage execution requires existing run state[/yellow]")


@app.command()
def report(
    run_id: str = typer.Option(..., "--run-id", "-r", help="Pipeline run ID"),  # noqa: B008
    format: str = typer.Option("json", "--format", "-f", help="Output format (json, html, csv)"),  # noqa: B008
    output_dir: str = typer.Option("outputs", "--output", "-o"),  # noqa: B008
):
    """
    Generate or view a pipeline report.

    Example:
        ddpipe report --run-id abc123 --format html
    """
    report_file = Path(output_dir) / f"report_{run_id}.json"

    if not report_file.exists():
        console.print(f"[red]Report not found: {report_file}[/red]")
        raise typer.Exit(1)

    with open(report_file) as f:
        report = json.load(f)

    console.print(Panel.fit(
        f"[bold]Pipeline Report[/bold]\n\n"
        f"Run ID: [cyan]{run_id}[/cyan]\n"
        f"Target: [green]{report.get('target', {}).get('gene', 'Unknown')}[/green]\n"
        f"Molecules Generated: {report.get('summary', {}).get('molecules_generated', 0)}\n"
        f"Molecules Docked: {report.get('summary', {}).get('molecules_docked', 0)}\n"
        f"Completed: {report.get('completed_at', 'Unknown')}",
    ))

    if format == "json":
        console.print_json(json.dumps(report, indent=2, default=str))


@app.command("export")
def export_results(
    run_id: str = typer.Option(..., "--run-id", "-r", help="Pipeline run ID"),  # noqa: B008
    format: str = typer.Option("sdf", "--format", "-f", help="Export format (sdf, csv, smiles)"),  # noqa: B008
    output_dir: str = typer.Option("outputs", "--output", "-o"),  # noqa: B008
    top_n: int = typer.Option(10, "--top", "-n", help="Number of candidates to export"),  # noqa: B008
):
    """
    Export pipeline results.

    Example:
        ddpipe export --run-id abc123 --format sdf --top 10
    """
    report_file = Path(output_dir) / f"report_{run_id}.json"

    if not report_file.exists():
        console.print(f"[red]Report not found: {report_file}[/red]")
        raise typer.Exit(1)

    with open(report_file) as f:
        report = json.load(f)

    candidates = report.get("top_candidates", [])[:top_n]

    if format == "smiles":
        output_file = Path(output_dir) / f"candidates_{run_id}.smi"
        with open(output_file, 'w') as f:
            for c in candidates:
                f.write(f"{c['smiles']}\t{c['molecule_id']}\n")
        console.print(f"[green]Exported to: {output_file}[/green]")

    elif format == "csv":
        output_file = Path(output_dir) / f"candidates_{run_id}.csv"
        with open(output_file, 'w') as f:
            f.write("rank,id,smiles,score,docking_score,mw,logp\n")
            for c in candidates:
                props = c.get('properties', {})
                f.write(f"{c['rank']},{c['molecule_id']},{c['smiles']},"
                        f"{c['composite_score']},{c['docking_score']},"
                        f"{props.get('molecular_weight', '')},{props.get('logP', '')}\n")
        console.print(f"[green]Exported to: {output_file}[/green]")

    elif format == "sdf":
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            output_file = Path(output_dir) / f"candidates_{run_id}.sdf"
            writer = Chem.SDWriter(str(output_file))

            for c in candidates:
                mol = Chem.MolFromSmiles(c['smiles'])
                if mol:
                    mol.SetProp("_Name", c['molecule_id'])
                    mol.SetProp("Rank", str(c['rank']))
                    mol.SetProp("CompositeScore", str(c['composite_score']))
                    mol.SetProp("DockingScore", str(c['docking_score']))

                    AllChem.EmbedMolecule(mol, randomSeed=42)
                    writer.write(mol)

            writer.close()
            console.print(f"[green]Exported to: {output_file}[/green]")

        except ImportError:
            console.print("[red]RDKit required for SDF export[/red]")
            raise typer.Exit(1) from None


@app.command()
def demo():
    """
    Run the VCP/FTD demo pipeline.

    This demonstrates the complete pipeline with the VCP target
    for frontotemporal dementia drug discovery.
    """
    console.print(Panel.fit(
        "[bold blue]VCP/FTD Demo Pipeline[/bold blue]\n\n"
        "Target: VCP (p97/Valosin-containing protein)\n"
        "Disease: Frontotemporal Dementia\n"
        "Reference: CB-5083 (VCP inhibitor)",
        title="Demo Mode"
    ))

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task("Running VCP demo pipeline...", total=None)

        try:
            result = run_vcp_demo_pipeline(output_dir="outputs/vcp_demo")
            console.print("\n[bold green]Demo Complete![/bold green]")
            console.print(f"Run ID: {result.run_id}")
            console.print(f"Status: {result.status}")
            console.print(f"Molecules: {result.molecules_generated}")
            console.print(f"Report: outputs/vcp_demo/report_{result.run_id}.json")

        except Exception as e:
            console.print(f"[red]Demo failed: {e}[/red]")
            raise typer.Exit(1) from None


@app.callback()
def main():
    """
    Drug Discovery Pipeline CLI (ddpipe)

    A command-line tool for AI-driven drug discovery,
    from genomic targets to molecule candidates.
    """
    pass


if __name__ == "__main__":
    app()
