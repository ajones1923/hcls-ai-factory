#!/usr/bin/env python3
"""
HLS Pipeline Runner
Python-based orchestrator for Genomics → RAG Chat → Drug Discovery

Alternative to Nextflow for environments with cgroup issues.
"""

import argparse
import json
import sys
import os
from pathlib import Path
from datetime import datetime
import subprocess
import shutil

# Add paths. Both the engine root AND its src/ go on the path: molecule_generator imports
# cleanly as a top-level module, but nim_clients uses package-relative imports
# (`from .models import ...`) and must be reached as `src.nim_clients`.
_ENGINE_DIR = Path(__file__).parent.parent / "core/engines/therapeutic-discovery"
sys.path.insert(0, str(_ENGINE_DIR / "src"))
sys.path.insert(0, str(_ENGINE_DIR))

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"


def print_banner():
    """Print pipeline banner."""
    print("""
╔═══════════════════════════════════════════════════════════════╗
║  Healthcare & Life Sciences Pipeline                          ║
║  ─────────────────────────────────────────────────────────── ║
║  Genomics → Target Discovery → Drug Design                    ║
╚═══════════════════════════════════════════════════════════════╝
    """)


def print_stage(stage_num: int, name: str, message: str = ""):
    """Print stage progress."""
    print(f"\n{'─' * 60}")
    print(f"  Stage {stage_num}: {name}")
    if message:
        print(f"  {message}")
    print(f"{'─' * 60}")


def _mol_to_dict(m):
    """Coerce a molecule to the dict shape the downstream stages expect. The real
    generate_vcp_molecules() returns GeneratedMolecule dataclass objects; the mock fallback
    returns dicts. (Pre-existing bug: stages used .get() on objects → AttributeError.)

    The generator deliberately labels its output -- `properties["provenance"] =
    "hand_written_analogue_library"`, `properties["score_is_placeholder"] = True`,
    `generation_method="analogue_library"` -- and this function used to keep only
    id/smiles/target_gene/score/properties, after which the candidate rows kept no provenance
    field at all. So report.json's `top_candidates` and report.html's "Top Drug Candidates"
    presented a hand-written CB-5083 analogue library, ranked by a fixed `0.9 - 0.05*i`
    placeholder, with nothing marking either fact. The labels are now carried through.
    """
    if isinstance(m, dict):
        return m
    props = getattr(m, 'properties', {}) or {}
    return {
        'id': getattr(m, 'name', None) or f"mol-{abs(hash(getattr(m, 'smiles', ''))) & 0xffffff:06x}",
        'smiles': getattr(m, 'smiles', ''),
        'target_gene': getattr(m, 'target_gene', ''),
        'score': getattr(m, 'score', 0.0),
        'generation_method': getattr(m, 'generation_method', None),
        'provenance': props.get('provenance', 'not recorded'),
        'score_is_placeholder': bool(props.get('score_is_placeholder', False)),
        'properties': props,
    }


def _qed(mol):
    """QED for a molecule, computed from its SMILES, or None.

    `_calculate_properties` in the generator computes MW, logP, hbd, hba, tpsa, rotatable bonds
    and Lipinski violations -- never `qed`. The ranking read `properties.get('qed', 0.4)`, so
    every candidate's QED was the literal default 0.4, half of every composite score with no
    docking was that constant, and report.html printed "0.400" for every row, to three decimals.
    """
    cached = (mol.get('properties') or {}).get('qed')
    if isinstance(cached, (int, float)):
        return float(cached)
    smiles = mol.get('smiles') or ''
    if not smiles:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import QED
    except ImportError:
        return None
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        return None
    return round(QED.qed(rdmol), 3)


def _find_structure_file(pdb_id: str):
    """Locate a cached PDB file for `pdb_id`, or None. Docking needs the coordinates; without
    them there is nothing to dock against and the stage must report that, not invent scores."""
    for root in (DATA_DIR / "structures", SCRIPT_DIR.parent / "data" / "structures"):
        for name in (f"{pdb_id}.pdb", f"{pdb_id.lower()}.pdb", f"pdb{pdb_id.lower()}.ent"):
            p = root / name
            if p.is_file():
                return p
    return None


def _run_docking(molecules, target):
    """Stage 7: dock with the real DiffDock client, or report that docking did not happen.

    Returns (docking_results, status). `status['performed']` is False whenever no docking
    actually ran; in that case `docking_results` is EMPTY -- never filled with placeholders.
    """
    status = {'performed': False, 'method': None, 'reason': None, 'message': ''}

    pdb_ids = target.get('pdb_ids') or []
    structure_file = next((f for f in (_find_structure_file(p) for p in pdb_ids) if f), None)
    if not structure_file:
        status['reason'] = (
            f"no cached structure file for {', '.join(pdb_ids) or 'the target'} under "
            f"{DATA_DIR / 'structures'}"
        )
        status['message'] = f"SKIPPED - {status['reason']}"
        return [], status

    try:
        from src.nim_clients import create_nim_clients, contains_mock
        nim = create_nim_clients()          # raises unless a real NIM/cloud backend is set up
    except Exception as e:
        status['reason'] = f"DiffDock unavailable ({type(e).__name__}: {e})"
        status['message'] = f"SKIPPED - {status['reason']}"
        return [], status

    pdb_content = structure_file.read_text()
    results = []
    for mol in molecules:
        mol_id = mol.get('id', mol.get('smiles', '')[:8])
        try:
            poses = nim.diffdock.dock(protein_pdb=pdb_content, ligand_smiles=mol.get('smiles', ''))
        except Exception as e:
            print(f"  docking failed for {mol_id}: {e}")
            continue
        if not poses:
            continue
        # NIM_ALLOW_MOCK_FALLBACK=true substitutes MockDiffDockClient, whose poses are random
        # numbers stamped `_provenance="mock"`. _run_docking used to read that flag into its
        # result dict and the ranking stage then dropped it, so report.html showed a green
        # "Docking: MockDiffDockClient" banner over simulated scores. A skipped stage is honest;
        # a simulated one presented as performed is not.
        if contains_mock(poses):
            status['method'] = type(nim.diffdock).__name__
            status['reason'] = (
                f"the docking backend is mock ({status['method']}): its poses are simulated "
                "numbers, not a docking. Set NIM_MODE=cloud with NVIDIA_API_KEY, or start the "
                "local NIM containers, for a docking term to be included."
            )
            status['message'] = f"SKIPPED - {status['reason']}"
            status['performed'] = False
            return [], status
        best = poses[0]
        results.append({
            'molecule_id': mol_id,
            'structure': structure_file.stem.upper(),
            'score_type': best.get('score_type', 'binding_energy'),
            'docking_score': best.get('docking_score'),
            'confidence_score': best.get('confidence_score'),
            'confidence': best.get('confidence'),
            'provenance': best.get('_provenance', 'real'),
        })

    status['performed'] = bool(results)
    status['method'] = type(nim.diffdock).__name__
    status['message'] = (f"Docked with {status['method']} against {structure_file.stem.upper()}"
                         if results else "SKIPPED - every docking call failed")
    return results, status


def _init_event_coordination():
    """BUS-10 (2026-06-29): wire the bidirectional cross-stage trigger engine into the
    orchestrator. `create_trigger_manager` registers the 4 built-in rules and SUBSCRIBES the
    manager to the shared EventBus; as the orchestrator publishes stage events below, the
    manager evaluates its rules and fires the cross-stage handoffs (RAG->Drug, Drug->Report).
    Returns (manager, emit_fn, EventType, PipelineStage). Fail-safe: no-ops if the lib is absent."""
    try:
        lib = str(Path(__file__).resolve().parent.parent / "lib")
        if lib not in sys.path:
            sys.path.insert(0, lib)
        from hcls_common.bidirectional_triggers import create_trigger_manager
        from hcls_common.event_bus import get_event_bus, PipelineEvent, EventType, PipelineStage
        bus = get_event_bus()
        manager = create_trigger_manager(event_bus=bus)   # 4 builtin rules + subscribe + start

        def emit(event_type, source_stage, payload=None):
            try:
                bus.publish(PipelineEvent(event_type=event_type, source_stage=source_stage,
                                          patient_id="HG002", payload=payload or {}))
            except Exception:
                pass

        print(f"  [triggers] cross-stage engine wired: "
              f"{manager.get_stats()['total_rules']} rules listening on the EventBus")
        return manager, emit, EventType, PipelineStage
    except Exception as e:  # never block the pipeline on the coordination layer
        print(f"  [triggers] event-coordination unavailable ({type(e).__name__}); "
              "running without cross-stage triggers")
        return None, (lambda *a, **k: None), None, None


def run_demo_pipeline(num_molecules: int = 20, output_dir: Path = None):
    """Run the VCP/FTD demo pipeline."""

    output_dir = output_dir or RESULTS_DIR / f"demo_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  Mode:          Demo (VCP/FTD)")
    print(f"  Output:        {output_dir}")
    print(f"  Molecules:     {num_molecules}")

    # Load demo target
    target_file = DATA_DIR / "demo" / "vcp_target.json"
    with open(target_file) as f:
        target = json.load(f)

    print_stage(0, "Initialize", f"Target: {target['gene']} ({target['protein']})")

    # BUS-10: stand up the cross-stage trigger engine; it reacts to the events published below.
    triggers, emit_event, ET, PS = _init_event_coordination()

    # Stage 1: Normalize target
    print_stage(1, "Normalize Target", f"UniProt: {target['uniprot_id']}")

    if ET:  # target identified -> fires the rag_to_drug_discovery rule (needs a druggable target)
        emit_event(ET.TARGETS_IDENTIFIED, PS.RAG_CHAT,
                   {"targets": [{"gene": target['gene'], "uniprot_id": target['uniprot_id'],
                                 "druggable": True, "pdb_ids": target.get('pdb_ids', [])}]})

    # Stage 2: Structure discovery
    print_stage(2, "Structure Discovery", f"PDB IDs: {', '.join(target['pdb_ids'])}")

    structures = []
    for pdb_id in target['pdb_ids']:
        local = _find_structure_file(pdb_id)
        structures.append({
            'pdb_id': pdb_id,
            # All four VCP entries used here are cryo-EM per their RCSB records (checked
            # 2026-11-03); 5FTK was singled out as X-ray and is not.
            'method': 'Cryo-EM',
            # `prepared: True` used to be hard-coded for every entry although this stage does
            # no preparation at all. Report what is actually on disk instead.
            'local_file': str(local) if local else None,
            'prepared': False,
            'preparation': 'not performed by this orchestrator',
        })

    # Stage 3: Structure prep
    n_local = len([s for s in structures if s['local_file']])
    print_stage(3, "Structure Preparation",
                f"{n_local}/{len(structures)} structures cached locally; "
                f"no preparation step is implemented")

    # Stage 4: Molecule generation
    print_stage(4, "Molecule Generation", f"Generating {num_molecules} molecules...")

    # The former fallback invented molecules by string-replacing the seed SMILES and drew
    # molecular_weight / logP / qed from random.uniform() -- numbers that belonged to no
    # molecule, printed as "(mock)" but written into report.json indistinguishably from real
    # output. If the generator cannot be imported, that is a setup fault; say so and stop.
    try:
        from molecule_generator import generate_vcp_molecules
    except ImportError as e:
        raise RuntimeError(
            "Cannot import the molecule generator from "
            "core/engines/therapeutic-discovery/src -- the demo pipeline cannot run. "
            f"({e}). Install the engine's requirements (RDKit in particular) and retry; "
            "this stage will not fabricate molecules."
        ) from e

    molecules = generate_vcp_molecules(num_molecules)
    if not molecules:
        raise RuntimeError(
            f"The molecule generator returned 0 molecules for {target['gene']}. "
            "Refusing to continue with invented structures."
        )
    print(f"  Generated {len(molecules)} molecules")

    # Coerce generator output to dicts (real generator returns GeneratedMolecule objects)
    molecules = [_mol_to_dict(m) for m in molecules]

    # Stage 5: Chemistry QC
    #
    # This was `m.get('properties', {}).get('molecular_weight', 500) < 550`: a molecule with no
    # computed MW got the literal default 500, silently passed, and the stage printed "N
    # molecules passed QC" -- a filter reporting a pass for a molecule it never examined. The
    # header also said "Filtering by drug-likeness" while applying one of the four Lipinski
    # terms. A molecule with no MW is now not examined and not passed.
    print_stage(5, "Chemistry QC", "Filtering on molecular weight (MW < 550 Da) only")
    passed, unexamined = [], []
    for m in molecules:
        mw = (m.get('properties') or {}).get('molecular_weight')
        if not isinstance(mw, (int, float)):
            unexamined.append(m.get('id', m.get('smiles', '')[:8]))
            continue
        if mw < 550:
            passed.append(m)
    print(f"  {len(passed)}/{len(molecules)} molecules passed the MW filter")
    if unexamined:
        print(f"  {len(unexamined)} molecules had no computed molecular weight and were NOT "
              f"examined (and so did not pass): {', '.join(unexamined[:5])}"
              f"{' ...' if len(unexamined) > 5 else ''}")
    if not passed:
        raise RuntimeError(
            f"No molecule could be examined by chemistry QC ({len(unexamined)} had no computed "
            "molecular weight). Install RDKit so properties are computed; this stage will not "
            "pass molecules it did not examine."
        )

    # Stage 6: Conformers
    print_stage(6, "Conformer Generation", f"Generating 3D conformers for {len(passed)} molecules")

    # Stage 7: Docking
    #
    # This stage used to print "Running DiffDock..." and then fill docking_score from
    # random.uniform(-12, -6) and confidence from random.uniform(0.6, 0.95) -- fabricated
    # binding numbers for a named tool that was never called, written into report.json and
    # report.html with no marker. Docking is now either real or absent.
    docking_results, docking_status = _run_docking(passed, target)
    print_stage(7, "Molecular Docking", docking_status['message'])
    print(f"  Docked {len(docking_results)} molecules")

    # Stage 8: Ranking
    print_stage(8, "Ranking Candidates", "Computing composite scores...")

    candidates = []
    dock_lookup = {d['molecule_id']: d for d in docking_results}

    # QED is computed, not defaulted. When RDKit cannot give one, the term is dropped from the
    # composite rather than filled in -- a composite built on a constant is not a ranking.
    qed_by_id = {}
    for mol in passed:
        qed_by_id[mol.get('id', mol.get('smiles', '')[:8])] = _qed(mol)
    have_qed = all(v is not None for v in qed_by_id.values()) and bool(qed_by_id)
    if not have_qed:
        print("  QED unavailable (RDKit missing or unparseable SMILES); "
              "ranking on generation score" + (" and docking" if docking_results else "") +
              " only — no QED term is included")

    # Without docking the composite cannot carry a docking term. Rather than substitute a
    # default score (the old code used -5 kcal/mol for every undocked molecule, which simply
    # added a constant), renormalize over the terms we actually have.
    if docking_results and have_qed:
        w_gen, w_dock, w_qed = 0.3, 0.4, 0.3
    elif docking_results:
        w_gen, w_dock, w_qed = 0.5, 0.5, 0.0
    elif have_qed:
        w_gen, w_dock, w_qed = 0.5, 0.0, 0.5
    else:
        w_gen, w_dock, w_qed = 1.0, 0.0, 0.0

    for mol in passed:
        mol_id = mol.get('id', mol.get('smiles', '')[:8])
        dock = dock_lookup.get(mol_id, {})

        gen_score = mol.get('score', 0.5)
        qed = qed_by_id.get(mol_id)

        dock_norm = 0.0
        dock_score = dock.get('docking_score')
        if dock.get('score_type') == 'diffdock_confidence':
            dock_norm = max(0.0, min(1.0, dock.get('confidence', 0.0)))
        elif dock_score is not None:
            # kcal/mol: -12 (excellent) -> 1.0, 0 (no binding) -> 0.0
            dock_norm = max(0.0, min(1.0, -dock_score / 12.0))

        composite = w_gen * gen_score + w_dock * dock_norm + w_qed * (qed or 0.0)

        candidates.append({
            'molecule_id': mol_id,
            'smiles': mol.get('smiles', ''),
            'score_type': dock.get('score_type'),
            'docking_score': dock_score,
            'pose_confidence': dock.get('confidence'),
            'docking_performed': bool(dock),
            'docking_provenance': dock.get('provenance', 'none'),
            'qed_score': qed,                       # None when it could not be computed
            # Carried from the generator so no surface can present a hand-written library with
            # a placeholder rank as generative-model output.
            'generation_method': mol.get('generation_method'),
            'provenance': mol.get('provenance', 'not recorded'),
            'score_is_placeholder': mol.get('score_is_placeholder', False),
            'composite_score': round(composite, 4),
            'composite_weights': {'generation': w_gen, 'docking': w_dock, 'qed': w_qed},
        })

    candidates.sort(key=lambda x: x['composite_score'], reverse=True)
    for i, c in enumerate(candidates):
        c['rank'] = i + 1

    print(f"  Ranked {len(candidates)} candidates")

    if ET:  # candidates ready -> fires the drug_to_report rule (needs a candidates list)
        emit_event(ET.DRUG_CANDIDATES_READY, PS.DRUG_DISCOVERY,
                   {"candidates": candidates[:10], "n_candidates": len(candidates)})

    # Stage 9: Reporting
    print_stage(9, "Report Generation", f"Writing to {output_dir}")

    report = {
        'pipeline': 'HLS Pipeline',
        'version': '1.0.0',
        'mode': 'demo',
        'target': target,
        'summary': {
            'molecules_generated': len(molecules),
            'passed_qc': len(passed),
            'top_candidates': min(10, len(candidates))
        },
        'docking': docking_status,          # performed / method / reason -- never implied
        'top_candidates': candidates[:10],
        'structures': structures,
        'completed_at': datetime.now().isoformat()
    }

    report_file = output_dir / 'report.json'
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)

    # Generate HTML report
    html_file = output_dir / 'report.html'
    generate_html_report(report, html_file)

    if ET:  # report done -> emit REPORT_GENERATED (terminal event)
        emit_event(ET.REPORT_GENERATED, PS.REPORTING, {"report": str(report_file)})

    # "SUCCESS" used to print unconditionally. A run that skipped docking completed, but it
    # did not do what the nine stage headers imply, so say which.
    skipped = [] if docking_status['performed'] else ['docking']
    print(f"\n{'═' * 60}")
    print(f"  PIPELINE COMPLETE")
    print(f"{'═' * 60}")
    if skipped:
        print(f"  Status:        COMPLETED WITH STAGES SKIPPED ({', '.join(skipped)})")
        print(f"  Skipped:       {docking_status['reason']}")
        print(f"  Ranking:       generation + QED only — NO binding prediction included")
    else:
        print(f"  Status:        SUCCESS")
    print(f"  Target:        {target['gene']}")
    print(f"  Molecules:     {len(candidates)}")
    print(f"  Top Score:     {candidates[0]['composite_score']:.4f}")
    print(f"  Report:        {report_file}")
    print(f"  HTML:          {html_file}")
    print(f"{'═' * 60}\n")

    if triggers:  # BUS-10: report what the cross-stage engine did, then unsubscribe cleanly
        st = triggers.get_stats()
        print(f"  Cross-stage triggers fired: {st['total_fired']} "
              f"({st['successful_fires']} ok / {st['failed_fires']} failed)")
        triggers.stop()

    return report


def generate_html_report(report: dict, output_file: Path):
    """Generate HTML report."""

    target = report.get('target', {})
    candidates = report.get('top_candidates', [])

    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>HLS Pipeline Report</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1000px; margin: 0 auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #1a1a2e; border-bottom: 3px solid #76b900; padding-bottom: 10px; }}
        h2 {{ color: #16213e; margin-top: 30px; }}
        .summary {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }}
        .card {{ background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%); color: white; padding: 20px; border-radius: 10px; text-align: center; }}
        .card h3 {{ margin: 0 0 10px 0; font-size: 14px; opacity: 0.9; }}
        .card .value {{ font-size: 32px; font-weight: bold; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #1a1a2e; color: white; }}
        tr:hover {{ background: #f5f5f5; }}
        .badge {{ display: inline-block; padding: 4px 12px; border-radius: 20px; font-size: 12px; background: #76b900; color: white; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Healthcare & Life Sciences Pipeline Report</h1>
        <p><strong>Generated:</strong> {report.get('completed_at', 'N/A')}</p>

        <div class="summary">
            <div class="card">
                <h3>Target Gene</h3>
                <div class="value">{target.get('gene', 'N/A')}</div>
            </div>
            <div class="card">
                <h3>Molecules Generated</h3>
                <div class="value">{report.get('summary', {}).get('molecules_generated', 0)}</div>
            </div>
            <div class="card">
                <h3>Top Candidates</h3>
                <div class="value">{report.get('summary', {}).get('top_candidates', 0)}</div>
            </div>
        </div>

        <h2>Target Information</h2>
        <table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Gene</td><td>{target.get('gene', 'N/A')}</td></tr>
            <tr><td>Protein</td><td>{target.get('protein', 'N/A')}</td></tr>
            <tr><td>UniProt ID</td><td>{target.get('uniprot_id', 'N/A')}</td></tr>
            <tr><td>Mechanism</td><td>{target.get('mechanism', 'N/A')}</td></tr>
            <tr><td>PDB Structures</td><td>{', '.join(target.get('pdb_ids', []))}</td></tr>
        </table>

        <h2>Top Drug Candidates</h2>'''

    # State plainly whether docking ran. The report used to show a "Docking Score" column
    # filled from random numbers with nothing distinguishing it from a real run.
    dock = report.get('docking', {})
    if dock.get('performed'):
        html += f'''
        <p style="padding:10px;border-left:4px solid #76b900;background:#f3f8e8;">
            Docking: <strong>{dock.get('method', 'unknown method')}</strong>.
        </p>'''
    else:
        html += f'''
        <p style="padding:10px;border-left:4px solid #b00020;background:#fdecee;">
            <strong>Docking was not performed</strong> ({dock.get('reason', 'reason not recorded')}).
            Candidates below are ranked on generation score and QED only; no binding
            prediction of any kind is included.
        </p>'''

    # The generator's provenance labels used to stop at this table: every row read as a ranked
    # drug candidate, and the QED column printed a hard-coded 0.400 to three decimals.
    if any(c.get('score_is_placeholder') for c in candidates):
        html += '''
        <p style="padding:10px;border-left:4px solid #e08700;background:#fff8ec;">
            <strong>Scores are placeholder ranks, not predictions.</strong> These candidates come
            from a hand-written analogue library; the generation score is a fixed rank order, so
            the ordering below carries no model confidence.
        </p>'''

    html += '''
        <table>
            <tr>
                <th>Rank</th>
                <th>Molecule ID</th>
                <th>Docking</th>
                <th>QED</th>
                <th>Composite Score</th>
                <th>Provenance</th>
            </tr>'''

    for c in candidates[:10]:
        if c.get('score_type') == 'diffdock_confidence':
            dock_cell = f"pose conf. {c.get('pose_confidence'):.3f}" if c.get('pose_confidence') is not None else 'n/a'
        elif c.get('docking_score') is not None:
            dock_cell = f"{c['docking_score']} kcal/mol"
        else:
            dock_cell = 'not docked'
        qed = c.get('qed_score')
        qed_cell = f"{qed:.3f}" if isinstance(qed, (int, float)) else 'not computed'
        composite = c.get('composite_score')
        composite_cell = f"{composite:.4f}" if isinstance(composite, (int, float)) else 'n/a'
        prov = c.get('provenance', 'not recorded')
        if c.get('score_is_placeholder'):
            prov += ' (placeholder rank)'
        html += f'''
            <tr>
                <td><span class="badge">{c.get('rank', 'N/A')}</span></td>
                <td>{c.get('molecule_id', 'N/A')}</td>
                <td>{dock_cell}</td>
                <td>{qed_cell}</td>
                <td><strong>{composite_cell}</strong></td>
                <td>{prov}</td>
            </tr>'''

    html += '''
        </table>

        <p style="margin-top: 40px; color: #666; font-size: 12px;">
            Generated by HLS Pipeline v1.0.0. Research use only - not for clinical decision-making.
        </p>
    </div>
</body>
</html>'''

    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(description='HLS Pipeline Runner')
    parser.add_argument('--mode', choices=['demo', 'full', 'target', 'drug'], default='demo',
                       help='Pipeline mode')
    parser.add_argument('--molecules', type=int, default=20,
                       help='Number of molecules to generate')
    parser.add_argument('--output', type=Path, default=None,
                       help='Output directory')
    parser.add_argument('--input', type=Path, default=None,
                       help='Input file (FASTQ, VCF, or target JSON)')

    args = parser.parse_args()

    print_banner()

    if args.mode == 'demo':
        run_demo_pipeline(num_molecules=args.molecules, output_dir=args.output)
    else:
        print(f"Mode '{args.mode}' requires additional implementation.")
        print("Use --mode demo for the VCP/FTD demonstration.")


if __name__ == '__main__':
    main()
