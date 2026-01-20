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
        'gene': primary_target.get('gene', 'VCP'),
        'protein': primary_target.get('protein', 'Unknown'),
        'uniprot_id': primary_target.get('uniprot_id'),
        'pdb_ids': primary_target.get('pdb_ids', ['5FTK']),
        'reference_smiles': 'CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5',  # CB-5083
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

    seed_smiles = data.get('reference_smiles', 'CCO')
    num_molecules = data.get('config', {}).get('num_molecules', 20)
    gene = data.get('gene', 'VCP')

    # In production, call MolMIM NIM service
    # For now, generate mock molecules based on CB-5083 scaffold

    scaffold_variants = [
        seed_smiles,  # Original CB-5083
        'CC(C)c1ccc(Nc2ncc3c(ccn3C)n2)cc1C(=O)Nc1ccc(CN2CCOCC2)cc1',
        'CC(N)c1ccc(Nc2ncc3c(ccn3C)n2)cc1C(=O)Nc1ccc(CN2CCOCC2)cc1',
        'Cc1ccc(NC2=NC3=C(C=N2)N(C=C3)C)c(C(=O)Nc4ccc(CN5CCOCC5)cc4)c1',
        'CCc1ccc(Nc2ncc3c(ccn3C)n2)cc1C(=O)Nc1ccc(CN2CCOCC2)cc1',
        'CC(C)c1ccc(Nc2ncc3c(ccn3C)n2)c(F)c1C(=O)Nc1ccc(CN2CCOCC2)cc1',
    ]

    molecules = []
    for i in range(min(num_molecules, len(scaffold_variants) * 3)):
        smiles = scaffold_variants[i % len(scaffold_variants)]
        mol_id = hashlib.md5(f'{sample_id}-{i}-{smiles}'.encode()).hexdigest()[:8]

        molecules.append({
            'id': f'{sample_id}-mol-{mol_id}',
            'smiles': smiles,
            'target_gene': gene,
            'generation_score': round(random.uniform(0.7, 1.0), 3),
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
        'generation_config': data.get('config', {})
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
    target_gene = data.get('target_gene', 'VCP')

    # In production, call DiffDock NIM service
    # For now, generate mock docking scores

    docking_results = []
    for mol in molecules:
        docking_score = round(random.uniform(-12, -6), 2)
        confidence = round(random.uniform(0.6, 0.95), 2)

        docking_results.append({
            'molecule_id': mol['id'],
            'smiles': mol['smiles'],
            'structure_id': '5FTK',
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
        'structure': '5FTK',
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

        # Normalize docking score
        dock_norm = max(0, min(1, (10 + dock_score) / 20))

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
