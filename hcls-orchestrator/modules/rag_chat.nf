/*
 * Target Lookup Module
 * VCF → Target Hypothesis identification
 *
 * NOT RAG. The process below was named RAG_IDENTIFY_TARGETS, commented "Identify drug targets
 * using RAG", wired into main.nf as "Stage 2 (RAG Chat)" and wrote `gene_table_lookup: 1.0` into
 * versions.yml -- while contacting no Milvus, no embedding model and no LLM. It is an 11-gene
 * literal dictionary. Three further faults made its output fabricated rather than merely
 * limited, and all are fixed here (see the process body):
 *   1. `if gene in info` was a substring test against the raw INFO column, so VCP matched
 *      VCPIP1/VCPKMT, ALK matched ALKBH1/ALKBH5, GRN matched PGRN and KRAS matched KRASP1.
 *   2. The rsID clause sat inside the per-gene loop, so one variant carrying rs188935092
 *      emitted all eleven genes as targets.
 *   3. Every emitted target was stamped confidence 'high', priority 5, status 'validated'
 *      regardless of evidence, and a VCF with no hit at all defaulted to a VCP/FTD target with
 *      the same stamps and a CB-5083 reference compound -- so a healthy control's report read
 *      "Target Gene: VCP".
 * `params.confidence_threshold` and `params.max_targets` were never read.
 */

process RAG_EXTRACT_VARIANTS {
    tag "$sample_id"
    label 'process_low'
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("*.variants.json"), emit: variants
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import gzip

    variants = []
    vcf_file = "${vcf}"

    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\\t')
            if len(fields) >= 8:
                chrom, pos, rsid, ref, alt, qual, filt, info = fields[:8]
                variants.append({
                    'chrom': chrom,
                    'pos': int(pos),
                    'rsid': rsid if rsid != '.' else None,
                    'ref': ref,
                    'alt': alt,
                    'qual': float(qual) if qual != '.' else None,
                    'filter': filt,
                    'info': info
                })

    with open('${sample_id}.variants.json', 'w') as f:
        json.dump({'sample': '${sample_id}', 'variants': variants}, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  python: 3.10\\n')
    """
}

process LOOKUP_TARGETS_IN_GENE_TABLE {
    tag "$sample_id"
    label 'process_medium'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.outdir}/rag_chat/targets", mode: 'copy'

    input:
    tuple val(sample_id), path(variants)

    output:
    tuple val(sample_id), path("*.targets.json"), emit: targets
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    from pathlib import Path

    # Load variants
    with open('${variants}') as f:
        data = json.load(f)

    # Disease-gene associations (in production, use RAG with Claude + full knowledge.py)
    disease_genes = {
        'VCP': {
            'diseases': ['Frontotemporal Dementia', 'ALS', 'IBM'],
            'mechanism': 'AAA+ ATPase, protein quality control',
            'druggability': 'high',
            'pdb_ids': ['5FTK', '8OOI', '9DIL', '7K56'],
            'reference_smiles': 'CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5',
            'reference_drug': 'CB-5083',
        },
        'MAPT': {
            'diseases': ['Frontotemporal Dementia', 'Alzheimer'],
            'mechanism': 'Microtubule stabilization',
            'druggability': 'medium',
            'pdb_ids': ['6QJH', '5O3L'],
        },
        'GRN': {
            'diseases': ['Frontotemporal Dementia'],
            'mechanism': 'Growth factor signaling',
            'druggability': 'medium',
            'pdb_ids': ['2JYE'],
        },
        'EGFR': {
            'diseases': ['Non-small cell lung cancer', 'Glioblastoma'],
            'mechanism': 'Receptor tyrosine kinase signaling',
            'druggability': 'high',
            'pdb_ids': ['1M17', '4HJO', '5UG9'],
            'reference_smiles': 'C1=CC2=C(C=C1OCCOC)C(=NC=N2)NC3=CC=CC(=C3)C#C',
            'reference_drug': 'Erlotinib',
        },
        'KRAS': {
            'diseases': ['Pancreatic Cancer', 'Lung Cancer', 'Colorectal Cancer'],
            'mechanism': 'GTPase signaling',
            'druggability': 'high',
            'pdb_ids': ['4OBE', '6OIM'],
            'reference_smiles': 'C1CC(=O)N(C1)C2=NC(=NC(=C2F)N)N3CC(C3)O',
            'reference_drug': 'Sotorasib',
        },
        'BRAF': {
            'diseases': ['Melanoma', 'Colorectal Cancer', 'Thyroid Cancer'],
            'mechanism': 'Serine/threonine kinase',
            'druggability': 'high',
            'pdb_ids': ['1UWH', '4MNE', '5CSW'],
            'reference_smiles': 'CCCS(=O)(=O)NC1=CC(=C(C=C1F)C(=O)C2=CNC3=NC=C(C=C23)Cl)F',
            'reference_drug': 'Vemurafenib',
        },
        'ALK': {
            'diseases': ['Non-small cell lung cancer', 'Neuroblastoma'],
            'mechanism': 'Receptor tyrosine kinase',
            'druggability': 'high',
            'pdb_ids': ['2XP2', '4MKC'],
            'reference_smiles': 'CC(C1=CC(=CC=C1)C2=NN=C(O2)C3=CC(=CC=C3)Cl)N4CCN(CC4)CC5=CC=CC=N5',
            'reference_drug': 'Crizotinib',
        },
        'BRCA1': {
            'diseases': ['Breast Cancer', 'Ovarian Cancer'],
            'mechanism': 'Homologous recombination repair',
            'druggability': 'high',
            'pdb_ids': ['1JM7', '4Y18'],
        },
        'BTK': {
            'diseases': ['Chronic Lymphocytic Leukemia', 'Mantle Cell Lymphoma'],
            'mechanism': 'Non-receptor tyrosine kinase',
            'druggability': 'high',
            'pdb_ids': ['3GEN', '5P9J'],
            'reference_smiles': 'C=CC(=O)N1CCC(CC1)N2C3=C(C=C(C=C3)OC4=CC=CC=C4)C(=N2)C5=CC=C(C=C5)NC6=NC=CC(=N6)N7CCN(CC7)C',
            'reference_drug': 'Ibrutinib',
        },
        'JAK2': {
            'diseases': ['Myeloproliferative Neoplasms', 'Polycythemia Vera'],
            'mechanism': 'Non-receptor tyrosine kinase',
            'druggability': 'high',
            'pdb_ids': ['2B7A', '4YTH'],
        },
        'CFTR': {
            'diseases': ['Cystic Fibrosis'],
            'mechanism': 'Chloride/bicarbonate channel',
            'druggability': 'high',
            'pdb_ids': ['5UAK', '6MSM'],
            'reference_smiles': 'CC1=C(C=CC(=C1)OCC(C)(C)C2=CC(=C(C=C2)O)Br)C3=CC=CC=C3',
            'reference_drug': 'Ivacaftor',
        },
    }

    def genes_in(variant):
        '''Gene symbols annotated on a variant, parsed out of the structured field.

        `if gene in info` was a substring test against the whole raw INFO column, so VCP matched
        VCPIP1 and VCPKMT, ALK matched ALKBH1/ALKBH5, GRN matched PGRN and KRAS matched KRASP1 --
        a target for a gene the variant does not touch. snpEff ANN= and VEP CSQ= are
        pipe-delimited; the gene symbol is one field, so compare that field for equality.
        '''
        info = variant.get('info', '') or ''
        symbols = set()
        for key in ('ANN=', 'CSQ='):
            start = info.find(key)
            if start < 0:
                continue
            block = info[start + len(key):].split(';')[0]
            for annotation in block.split(','):
                fields = annotation.split('|')
                # snpEff ANN: allele|effect|impact|gene_name|... ; VEP CSQ: allele|effect|impact|symbol|...
                if len(fields) > 3 and fields[3].strip():
                    symbols.add(fields[3].strip().upper())
        if not symbols:
            # No structured annotation: a bare GENE=... key is the only other explicit form.
            for part in info.split(';'):
                if part.upper().startswith('GENE='):
                    symbols.add(part.split('=', 1)[1].strip().upper())
        return symbols

    # Identify target genes from variants. Every emitted target used to be stamped
    # confidence 'high', priority 5, status 'validated' regardless of evidence; those are now
    # derived from what the table and the annotation actually say, and params.max_targets /
    # params.confidence_threshold (which were never read) are applied.
    CONFIDENCE_SCORE = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
    threshold = float(${params.confidence_threshold})
    max_targets = int(${params.max_targets})

    by_gene = {}
    for variant in data.get('variants', []):
        for gene in genes_in(variant):
            if gene not in disease_genes:
                continue
            gene_info = disease_genes[gene]
            # Druggability is a property of the gene in the table, not of this variant, so it is
            # what the confidence derives from -- never an unconditional 'high'.
            druggability = gene_info['druggability']
            score = CONFIDENCE_SCORE.get(druggability, 0.0)
            if score < threshold:
                continue
            impact = str(variant.get('impact', '') or '').upper()
            entry = by_gene.setdefault(gene, {
                'gene': gene,
                'protein': f'{gene} protein',
                'uniprot_id': gene_info.get('uniprot_id'),
                'variants': [],
                'diseases': gene_info['diseases'],
                'mechanism': gene_info['mechanism'],
                'druggability': druggability,
                'pdb_ids': gene_info['pdb_ids'],
                'confidence': druggability,
                'confidence_score': score,
                # 'validated' is a claim about a target, not about a gene table lookup.
                'status': 'gene_table_match',
                'evidence': 'matched a gene symbol in the built-in 11-gene table; no retrieval, '
                            'no literature search and no model inference was performed',
            })
            entry['variants'].append(variant)
            if gene_info.get('reference_smiles'):
                entry['reference_smiles'] = gene_info['reference_smiles']
            if gene_info.get('reference_drug'):
                entry['reference_drug'] = gene_info['reference_drug']

    targets = sorted(by_gene.values(), key=lambda t: -t['confidence_score'])[:max_targets]

    # No match is no match. This used to append a VCP/FTD target stamped confidence 'high',
    # status 'validated' and seeded with CB-5083, so a VCF with no variant in any of the eleven
    # genes -- a healthy control -- produced a report reading "Target Gene: VCP" for that
    # patient. A built-in demo target now requires an explicit opt-in and says what it is.
    if not targets and '${params.demo_target_fallback}'.lower() == 'true':
        targets.append({
            'gene': 'VCP',
            'protein': 'Valosin-containing protein (p97)',
            'uniprot_id': 'P55072',
            'diseases': ['Frontotemporal Dementia'],
            'mechanism': 'AAA+ ATPase inhibition',
            'druggability': 'high',
            'pdb_ids': ['5FTK', '8OOI', '9DIL', '7K56'],
            'reference_smiles': 'CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5',
            'reference_drug': 'CB-5083',
            'confidence': 'not derived from this sample',
            'confidence_score': 0.0,
            'status': 'demo_default',
            'evidence': 'DEMO DEFAULT -- not identified from this VCF. No variant in this '
                        'sample matched any gene in the built-in table.',
            'rationale': 'VCP mutations cause FTD-ALS. CB-5083 is a known VCP inhibitor '
                         '(demo default, --demo_target_fallback true).'
        })

    output = {
        'sample': '${sample_id}',
        'targets': targets,
        'method': 'built-in 11-gene table lookup (not RAG: no vector search, no LLM)',
        'n_targets': len(targets),
        'confidence_threshold': threshold,
        'max_targets': max_targets,
        'analysis_version': '1.0'
    }

    with open('${sample_id}.targets.json', 'w') as f:
        json.dump(output, f, indent=2)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  gene_table_lookup: 1.0\\n')
    """
}

// Main target-lookup workflow (named RAG_CHAT_PIPELINE for the callers in main.nf)
workflow RAG_CHAT_PIPELINE {
    take:
    vcf  // tuple(sample_id, vcf_file)

    main:
    ch_versions = Channel.empty()

    // Handle both tuple and path inputs
    ch_vcf = vcf.map { it ->
        if (it instanceof List) {
            return it
        } else {
            return [it.baseName.replaceAll(/\\.vcf.*/, ''), it]
        }
    }

    // Extract variants from VCF
    RAG_EXTRACT_VARIANTS(ch_vcf)
    ch_versions = ch_versions.mix(RAG_EXTRACT_VARIANTS.out.versions)

    // Look the variants' genes up in the built-in gene table (no RAG involved)
    LOOKUP_TARGETS_IN_GENE_TABLE(RAG_EXTRACT_VARIANTS.out.variants)
    ch_versions = ch_versions.mix(LOOKUP_TARGETS_IN_GENE_TABLE.out.versions)

    emit:
    targets = LOOKUP_TARGETS_IN_GENE_TABLE.out.targets
    versions = ch_versions
}
