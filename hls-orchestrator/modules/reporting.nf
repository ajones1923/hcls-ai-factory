/*
 * Reporting Module
 * Generate final pipeline reports
 */

process GENERATE_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(ranked)
    tuple val(sample_id2), path(targets)

    output:
    path "*.report.json", emit: json
    path "*.report.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import json
    from datetime import datetime

    # Load data
    with open('${ranked}') as f:
        ranked_data = json.load(f)
    with open('${targets}') as f:
        targets_data = json.load(f)

    sample_id = '${sample_id}'

    # Build report
    report = {
        'pipeline': 'Healthcare & Life Sciences Pipeline',
        'version': '1.0.0',
        'sample_id': sample_id,
        'generated_at': datetime.now().isoformat(),

        'summary': {
            'target_gene': ranked_data.get('target_gene'),
            'total_molecules_generated': ranked_data.get('total_candidates', 0),
            'top_candidates': len(ranked_data.get('ranked_candidates', []))
        },

        'target': targets_data.get('targets', [{}])[0] if 'targets' in targets_data else targets_data,

        'top_candidates': ranked_data.get('ranked_candidates', [])[:5],

        'parameters': {
            'num_molecules': ${params.num_molecules},
            'diversity': ${params.diversity},
            'max_mw': ${params.max_mw},
            'docking_poses': ${params.docking_poses}
        }
    }

    # Save JSON report
    with open(f'{sample_id}.report.json', 'w') as f:
        json.dump(report, f, indent=2)

    # Generate HTML report
    html = f'''<!DOCTYPE html>
<html>
<head>
    <title>HLS Pipeline Report - {sample_id}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #1a1a2e; border-bottom: 3px solid #76b900; padding-bottom: 10px; }}
        h2 {{ color: #16213e; margin-top: 30px; }}
        .summary {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }}
        .card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 10px; }}
        .card.nvidia {{ background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%); }}
        .card h3 {{ margin: 0 0 10px 0; font-size: 14px; opacity: 0.9; }}
        .card .value {{ font-size: 32px; font-weight: bold; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #1a1a2e; color: white; }}
        tr:hover {{ background: #f5f5f5; }}
        .smiles {{ font-family: monospace; font-size: 11px; max-width: 300px; overflow: hidden; text-overflow: ellipsis; }}
        .score {{ font-weight: bold; }}
        .score.high {{ color: #76b900; }}
        .score.medium {{ color: #f0ad4e; }}
        .badge {{ display: inline-block; padding: 4px 12px; border-radius: 20px; font-size: 12px; background: #76b900; color: white; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Healthcare & Life Sciences Pipeline Report</h1>
        <p><strong>Sample:</strong> {sample_id} | <strong>Generated:</strong> {report["generated_at"]}</p>

        <div class="summary">
            <div class="card nvidia">
                <h3>Target Gene</h3>
                <div class="value">{report["summary"]["target_gene"]}</div>
            </div>
            <div class="card">
                <h3>Molecules Generated</h3>
                <div class="value">{report["summary"]["total_molecules_generated"]}</div>
            </div>
            <div class="card">
                <h3>Top Candidates</h3>
                <div class="value">{report["summary"]["top_candidates"]}</div>
            </div>
        </div>

        <h2>Target Information</h2>
        <table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Gene</td><td>{report["target"].get("gene", "N/A")}</td></tr>
            <tr><td>Protein</td><td>{report["target"].get("protein", "N/A")}</td></tr>
            <tr><td>UniProt ID</td><td>{report["target"].get("uniprot_id", "N/A")}</td></tr>
            <tr><td>Mechanism</td><td>{report["target"].get("mechanism", "N/A")}</td></tr>
            <tr><td>PDB Structures</td><td>{", ".join(report["target"].get("pdb_ids", []))}</td></tr>
        </table>

        <h2>Top Drug Candidates</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>Molecule ID</th>
                <th>SMILES</th>
                <th>Docking Score</th>
                <th>QED</th>
                <th>Composite</th>
            </tr>'''

    for candidate in report['top_candidates']:
        score_class = 'high' if candidate['composite_score'] > 0.4 else 'medium'
        html += f'''
            <tr>
                <td><span class="badge">{candidate["rank"]}</span></td>
                <td>{candidate["molecule_id"]}</td>
                <td class="smiles">{candidate["smiles"][:50]}...</td>
                <td>{candidate["docking_score"]:.2f}</td>
                <td>{candidate["qed_score"]:.3f}</td>
                <td class="score {score_class}">{candidate["composite_score"]:.4f}</td>
            </tr>'''

    html += '''
        </table>

        <h2>Pipeline Parameters</h2>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>'''

    for key, value in report['parameters'].items():
        html += f'<tr><td>{key}</td><td>{value}</td></tr>'

    html += '''
        </table>

        <p style="margin-top: 40px; color: #666; font-size: 12px;">
            Generated by HLS Pipeline v1.0.0 | Powered by NVIDIA BioNeMo
        </p>
    </div>
</body>
</html>'''

    with open(f'{sample_id}.report.html', 'w') as f:
        f.write(html)

    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n  reporting: 1.0\\n')
    """
}
