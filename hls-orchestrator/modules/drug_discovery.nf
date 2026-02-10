/*
 * Drug Discovery Pipeline Module
 * Target Hypothesis → Drug Candidates
 */

process DRUG_DISCOVERY_PREP_TARGET {
    tag "$sample_id"
    label 'process_low'

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
    publishDir "${params.outdir}/drug_discovery/molecules", mode: 'copy'

    input:
    tuple val(sample_id), path(prepared)

    output:
    tuple val(sample_id), path("*.molecules.json"), emit: molecules
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import random
    import hashlib
    from datetime import datetime

    with open('${prepared}') as f:
        data = json.load(f)

    seed_smiles = data.get('reference_smiles', '')
    num_molecules = data.get('config', {}).get('num_molecules', 20)
    gene = data.get('gene', 'UNKNOWN')

    # In production, call MolMIM NIM service
    # For now, generate mock molecules from the seed compound
    # If no seed provided, use a generic drug-like scaffold
    if not seed_smiles:
        seed_smiles = 'c1ccc2c(c1)cc(cn2)NC(=O)c1ccccc1'  # Generic quinoline amide scaffold

    # Simple atom-swap diversification for mock data
    # Each swap produces a structurally distinct analogue
    atom_swaps = [
        ('F', 'Cl'), ('Cl', 'F'), ('C(=O)', 'C(=S)'), ('O', 'S'),
        ('N', 'C'), ('c1ccc', 'c1cnc'), ('CC', 'CCC'), ('NH', 'N(C)'),
        ('F', 'Br'), ('OC', 'OCC'), ('(C)', '(CC)'), ('NC', 'NCC'),
    ]

    molecules = []
    seen_smiles = set()

    # First molecule is always the seed
    seen_smiles.add(seed_smiles)
    mol_id = hashlib.md5(f'{sample_id}-0-{seed_smiles}'.encode()).hexdigest()[:8]
    molecules.append({
        'id': f'{sample_id}-mol-{mol_id}',
        'smiles': seed_smiles,
        'target_gene': gene,
        'generation_score': 1.0,
        'method': 'MolMIM-seed',
        'properties': {
            'molecular_weight': round(random.uniform(420, 520), 1),
            'logP': round(random.uniform(3.0, 5.0), 2),
            'hbd': random.randint(1, 3),
            'hba': random.randint(5, 8),
            'tpsa': round(random.uniform(70, 110), 1),
            'rotatable_bonds': random.randint(5, 9),
            'qed': round(random.uniform(0.35, 0.55), 3)
        },
        'generated_at': datetime.now().isoformat()
    })

    # Generate analogues via atom swaps
    for i in range(1, num_molecules):
        swap = atom_swaps[i % len(atom_swaps)]
        smiles = seed_smiles.replace(swap[0], swap[1], 1)
        # If swap produced no change, try a different approach
        if smiles == seed_smiles or smiles in seen_smiles:
            # Insert a methyl group at a random bond position
            insert_pos = min(i * 3, len(seed_smiles) - 1)
            smiles = seed_smiles[:insert_pos] + 'C' + seed_smiles[insert_pos:]
        seen_smiles.add(smiles)
        mol_id = hashlib.md5(f'{sample_id}-{i}-{smiles}'.encode()).hexdigest()[:8]

        molecules.append({
            'id': f'{sample_id}-mol-{mol_id}',
            'smiles': smiles,
            'target_gene': gene,
            'generation_score': round(random.uniform(0.7, 0.95), 3),
            'method': 'MolMIM',
            'properties': {
                'molecular_weight': round(random.uniform(420, 520), 1),
                'logP': round(random.uniform(3.0, 5.0), 2),
                'hbd': random.randint(1, 3),
                'hba': random.randint(5, 8),
                'tpsa': round(random.uniform(70, 110), 1),
                'rotatable_bonds': random.randint(5, 9),
                'qed': round(random.uniform(0.35, 0.55), 3)
            },
            'generated_at': datetime.now().isoformat()
        })

    output = {
        'sample': '${sample_id}',
        'target_gene': gene,
        'molecules': molecules,
        'generation_config': {**data.get('config', {}), 'pdb_ids': data.get('pdb_ids', [])}
    }

    with open('${sample_id}.molecules.json', 'w') as f:
        json.dump(output, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  molmim: 1.0\\n')
    """
}

process DRUG_DISCOVERY_DOCKING {
    tag "$sample_id"
    label 'gpu'
    publishDir "${params.outdir}/drug_discovery/docking", mode: 'copy'

    input:
    tuple val(sample_id), path(molecules)

    output:
    tuple val(sample_id), path("*.docking.json"), emit: docking
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import random

    with open('${molecules}') as f:
        data = json.load(f)

    molecules = data.get('molecules', [])
    target_gene = data.get('target_gene', 'UNKNOWN')

    # Determine structure ID from generation config or first available PDB
    gen_config = data.get('generation_config', {})
    pdb_ids = gen_config.get('pdb_ids', [])
    structure_id = pdb_ids[0] if pdb_ids else 'NONE'

    # In production, call DiffDock NIM service
    # For now, generate mock docking scores

    docking_results = []
    for mol in molecules:
        docking_score = round(random.uniform(-12, -6), 2)
        confidence = round(random.uniform(0.6, 0.95), 2)

        docking_results.append({
            'molecule_id': mol['id'],
            'smiles': mol['smiles'],
            'structure_id': structure_id,
            'docking_score': docking_score,
            'confidence': confidence,
            'poses': ${params.docking_poses},
            'contacts': [
                f'ALA{random.randint(100, 500)}',
                f'GLY{random.randint(100, 500)}',
                f'ASP{random.randint(100, 500)}'
            ]
        })

    output = {
        'sample': '${sample_id}',
        'target_gene': target_gene,
        'structure': structure_id,
        'docking_results': docking_results
    }

    with open('${sample_id}.docking.json', 'w') as f:
        json.dump(output, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  diffdock: 1.0\\n')
    """
}

process DRUG_DISCOVERY_RANK {
    tag "$sample_id"
    label 'process_low'
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

    // Generate molecules with MolMIM
    DRUG_DISCOVERY_GENERATE_MOLECULES(DRUG_DISCOVERY_PREP_TARGET.out.prepared)
    ch_versions = ch_versions.mix(DRUG_DISCOVERY_GENERATE_MOLECULES.out.versions)

    // Dock molecules with DiffDock
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
