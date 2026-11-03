/*
 * Drug Discovery Pipeline Module
 * Target Hypothesis → Drug Candidates
 *
 * STATUS (2026-11-03): NOT RUNNABLE. Two processes below --
 * DRUG_DISCOVERY_GENERATE_MOLECULES and DRUG_DISCOVERY_DOCKING -- were mock
 * implementations that produced random numbers under the MolMIM and DiffDock names and
 * wrote versions.yml claiming those tools. The fabricated bodies have been replaced with
 * an explicit failure, so this workflow now fails fast instead of emitting invented
 * chemistry. Nextflow does not run on the DGX Spark in any case (see main.nf); the active
 * orchestrator is run_pipeline.py, which calls the real NIM clients in
 * core/engines/therapeutic-discovery/src/nim_clients.py.
 */

process DRUG_DISCOVERY_PREP_TARGET {
    tag "$sample_id"
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_id), path(targets_json)

    output:
    tuple val(sample_id), path("*.prepared.json"), emit: prepared
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json

    with open('${targets_json}') as f:
        data = json.load(f)

    # Prepare target for drug discovery
    targets = data.get('targets', [data]) if 'targets' in data else [data]
    primary_target = targets[0]

    prepared = {
        'sample': '${sample_id}',
        'gene': primary_target.get('gene', 'UNKNOWN'),
        'protein': primary_target.get('protein', 'Unknown'),
        'uniprot_id': primary_target.get('uniprot_id'),
        'pdb_ids': primary_target.get('pdb_ids', []),
        'reference_smiles': primary_target.get('reference_smiles', ''),
        'reference_drug': primary_target.get('reference_drug', ''),
        'config': {
            'num_molecules': ${params.num_molecules},
            'diversity': ${params.diversity},
            'max_mw': ${params.max_mw}
        }
    }

    with open('${sample_id}.prepared.json', 'w') as f:
        json.dump(prepared, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  drug_discovery: 1.0\\n')
    """
}

process DRUG_DISCOVERY_GENERATE_MOLECULES {
    tag "$sample_id"
    label 'gpu'
    // Was `ignore` on final attempt: a failed stage emitted nothing and the workflow carried
    // on as if it had succeeded. A missing model is a hard stop, not a warning.
    errorStrategy 'terminate'
    publishDir "${params.outdir}/drug_discovery/molecules", mode: 'copy'

    input:
    tuple val(sample_id), path(prepared)

    output:
    tuple val(sample_id), path("*.molecules.json"), emit: molecules
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    #
    # HONESTY FIX (2026-11-03). This process was labeled MolMIM and wrote
    # `versions.yml` claiming `molmim: 1.0`, but no MolMIM service was ever contacted. It
    # built "molecules" by string-replacing substrings of the seed SMILES -- producing
    # strings that are frequently not valid chemistry -- and filled molecular_weight, logP,
    # hbd, hba, tpsa, rotatable_bonds and qed from random.uniform()/random.randint(). The
    # downstream rank process then ordered candidates by those random numbers.
    #
    # There is a real MolMIM client in the repo:
    #   core/engines/therapeutic-discovery/src/nim_clients.py (local NIM + NVIDIA cloud)
    # It is driven by run_pipeline.py, which is the ACTIVE orchestrator (see main.nf). Wiring
    # it in here needs the engine's Python environment inside the task container, which this
    # workflow does not define. Until it does, this process FAILS rather than fabricating.
    import sys

    sys.stderr.write(
        "DRUG_DISCOVERY_GENERATE_MOLECULES is not implemented.\\n"
        "This process previously emitted randomly-generated 'molecules' labeled MolMIM.\\n"
        "That has been removed; there is no simulated path.\\n"
        "\\n"
        "To generate molecules for real, use the active orchestrator:\\n"
        "  cd hcls-orchestrator && python run_pipeline.py --mode demo\\n"
        "with a MolMIM backend configured (NIM_MODE=cloud + NVIDIA_API_KEY, or a local\\n"
        "NIM container on :8001). See core/engines/therapeutic-discovery/src/nim_clients.py.\\n"
    )
    sys.exit(1)
    """
}

process DRUG_DISCOVERY_DOCKING {
    tag "$sample_id"
    label 'gpu'
    // Was `ignore` on final attempt: a failed stage emitted nothing and the workflow carried
    // on as if it had succeeded. A missing model is a hard stop, not a warning.
    errorStrategy 'terminate'
    publishDir "${params.outdir}/drug_discovery/docking", mode: 'copy'

    input:
    tuple val(sample_id), path(molecules)

    output:
    tuple val(sample_id), path("*.docking.json"), emit: docking
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    #
    # HONESTY FIX (2026-11-03). This process was labeled DiffDock and wrote
    # `versions.yml` claiming `diffdock: 1.0`, but no DiffDock service was ever contacted.
    # Binding scores came from random.uniform(-12, -6) and were presented as kcal/mol; the
    # per-pose "contacts" were residue names with random.randint(100, 500) sequence numbers,
    # i.e. invented binding-site evidence for a protein whose structure was never read.
    #
    # The real client is core/engines/therapeutic-discovery/src/nim_clients.py
    # (DiffDockClient / CloudDiffDockClient), driven by run_pipeline.py. Note that DiffDock
    # predicts a POSE and a confidence, not a kcal/mol affinity -- the random numbers were
    # fabricating a quantity the named tool does not even produce.
    import sys

    sys.stderr.write(
        "DRUG_DISCOVERY_DOCKING is not implemented.\\n"
        "This process previously emitted random binding scores and invented contact\\n"
        "residues labeled DiffDock. That has been removed; there is no simulated path.\\n"
        "\\n"
        "To dock for real, use the active orchestrator:\\n"
        "  cd hcls-orchestrator && python run_pipeline.py --mode demo\\n"
        "with a DiffDock backend (NIM_MODE=cloud + NVIDIA_API_KEY, or a local NIM container\\n"
        "on :8002) and a cached PDB file for the target structure.\\n"
    )
    sys.exit(1)
    """
}

process DRUG_DISCOVERY_RANK {
    tag "$sample_id"
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.outdir}/drug_discovery/ranked", mode: 'copy'

    input:
    tuple val(sample_id), path(molecules), path(docking)

    output:
    tuple val(sample_id), path("*.ranked.json"), emit: ranked
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json

    with open('${molecules}') as f:
        mol_data = json.load(f)
    with open('${docking}') as f:
        dock_data = json.load(f)

    molecules = {m['id']: m for m in mol_data.get('molecules', [])}
    docking = {d['molecule_id']: d for d in dock_data.get('docking_results', [])}

    # Calculate composite scores and rank
    candidates = []
    for mol_id, mol in molecules.items():
        dock = docking.get(mol_id, {})

        gen_score = mol.get('generation_score', 0.5)
        dock_score = dock.get('docking_score', -5)
        qed = mol.get('properties', {}).get('qed', 0.4)

        # Normalize docking score (lower/more negative is better binding)
        # Range: -12 kcal/mol (excellent) -> 1.0, 0 kcal/mol (no binding) -> 0.0
        dock_norm = max(0, min(1, -dock_score / 12.0))

        # Composite score
        composite = 0.3 * gen_score + 0.4 * dock_norm + 0.3 * qed

        candidates.append({
            'molecule_id': mol_id,
            'smiles': mol['smiles'],
            'target_gene': mol.get('target_gene'),
            'generation_score': gen_score,
            'docking_score': dock_score,
            'qed_score': qed,
            'composite_score': round(composite, 4),
            'properties': mol.get('properties', {}),
            'contacts': dock.get('contacts', [])
        })

    # Sort by composite score
    candidates.sort(key=lambda x: x['composite_score'], reverse=True)

    # Assign ranks
    for i, c in enumerate(candidates):
        c['rank'] = i + 1

    output = {
        'sample': '${sample_id}',
        'target_gene': mol_data.get('target_gene'),
        'total_candidates': len(candidates),
        'ranked_candidates': candidates[:10]  # Top 10
    }

    with open('${sample_id}.ranked.json', 'w') as f:
        json.dump(output, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  ranking: 1.0\\n')
    """
}

// Main Drug Discovery workflow
workflow DRUG_DISCOVERY {
    take:
    targets  // tuple(sample_id, targets_json) or path

    main:
    ch_versions = Channel.empty()

    // Handle both tuple and path inputs
    ch_targets = targets.map { it ->
        if (it instanceof List) {
            return it
        } else {
            return [it.baseName.replaceAll(/\\..*/, ''), it]
        }
    }

    // Prepare target for drug discovery
    DRUG_DISCOVERY_PREP_TARGET(ch_targets)
    ch_versions = ch_versions.mix(DRUG_DISCOVERY_PREP_TARGET.out.versions)

    // Generate molecules (NOT IMPLEMENTED -- fails; see the process body)
    DRUG_DISCOVERY_GENERATE_MOLECULES(DRUG_DISCOVERY_PREP_TARGET.out.prepared)
    ch_versions = ch_versions.mix(DRUG_DISCOVERY_GENERATE_MOLECULES.out.versions)

    // Dock molecules (NOT IMPLEMENTED -- fails; see the process body)
    DRUG_DISCOVERY_DOCKING(DRUG_DISCOVERY_GENERATE_MOLECULES.out.molecules)
    ch_versions = ch_versions.mix(DRUG_DISCOVERY_DOCKING.out.versions)

    // Rank candidates
    DRUG_DISCOVERY_RANK(
        DRUG_DISCOVERY_GENERATE_MOLECULES.out.molecules
            .join(DRUG_DISCOVERY_DOCKING.out.docking)
    )
    ch_versions = ch_versions.mix(DRUG_DISCOVERY_RANK.out.versions)

    emit:
    molecules = DRUG_DISCOVERY_RANK.out.ranked
    versions = ch_versions
}
