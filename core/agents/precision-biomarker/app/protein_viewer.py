"""3D Protein Structure Viewer using Mol* (Molstar).

Embeds an interactive WebGL protein structure viewer in Streamlit
for visualizing drug target structures with mutation sites and
ligand binding pockets.

Usage:
    from protein_viewer import render_protein_viewer
    render_protein_viewer("5FTK")  # VCP/p97 structure
"""

from __future__ import annotations

import streamlit as st
import streamlit.components.v1 as components


def get_molstar_html(pdb_id: str, highlight_residue: int = 155, highlight_chain: str = "A",
                     width: int = 800, height: int = 500) -> str:
    """Generate HTML for Mol* viewer with a specific PDB structure.

    Args:
        pdb_id: PDB accession code (e.g., "5FTK")
        highlight_residue: Residue number to highlight (e.g., 155 for R155H)
        highlight_chain: Chain ID for the highlight
        width: Viewer width in pixels
        height: Viewer height in pixels
    """
    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <link rel="stylesheet" type="text/css" href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdbe-molstar-3.1.0.css">
        <script type="text/javascript" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-3.1.0.js"></script>
        <style>
            body {{ margin: 0; padding: 0; background: #0a0a0f; overflow: hidden; }}
            #viewer {{ width: {width}px; height: {height}px; position: relative; }}
            .info-panel {{
                position: absolute; bottom: 10px; left: 10px; z-index: 100;
                background: rgba(27, 27, 47, 0.9); color: white; padding: 10px 15px;
                border-radius: 8px; font-family: 'Segoe UI', sans-serif; font-size: 12px;
                border: 1px solid rgba(118, 185, 0, 0.3);
            }}
            .info-panel .label {{ color: #76B900; font-weight: bold; }}
            .info-panel .mutation {{ color: #ef4444; font-weight: bold; }}
        </style>
    </head>
    <body>
        <div id="viewer">
            <div class="info-panel">
                <span class="label">PDB:</span> {pdb_id.upper()} &nbsp;|&nbsp;
                <span class="label">Target:</span> VCP/p97 &nbsp;|&nbsp;
                <span class="mutation">Mutation:</span> R{highlight_residue}H (Chain {highlight_chain})
            </div>
        </div>
        <script>
            var viewerInstance = new PDBeMolstarPlugin();
            var options = {{
                moleculeId: '{pdb_id.lower()}',
                hideControls: false,
                bgColor: {{r: 10, g: 10, b: 15}},
                lighting: 'metallic',
                pdbeUrl: 'https://www.ebi.ac.uk/pdbe/',
                encoding: 'bcif',
                selectInteraction: true,
                landscape: true,
                hideCanvasControls: ['animation', 'expand'],
            }};

            var viewerContainer = document.getElementById('viewer');
            viewerInstance.render(viewerContainer, options);

            // Highlight mutation site after structure loads
            viewerInstance.events.loadComplete.subscribe(function() {{
                // Select and highlight the mutation residue
                viewerInstance.visual.select({{
                    data: [{{
                        struct_asym_id: '{highlight_chain}',
                        start_residue_number: {highlight_residue},
                        end_residue_number: {highlight_residue},
                        color: {{r: 239, g: 68, b: 68}},  // Red for mutation
                        focus: true,
                    }}],
                    nonSelectedColor: {{r: 180, g: 180, b: 180}},
                }});
            }});
        </script>
    </body>
    </html>
    """


def render_protein_viewer(pdb_id: str = "5FTK", show_controls: bool = True):
    """Render the 3D protein structure viewer in Streamlit.

    Args:
        pdb_id: PDB accession code
        show_controls: Whether to show PDB selection controls
    """
    st.markdown("### 🧬 3D Protein Structure Viewer")
    st.caption("Interactive visualization of drug target structures with mutation sites")

    if show_controls:
        col1, col2, col3 = st.columns([2, 1, 1])
        with col1:
            pdb_options = {
                "5FTK": "5FTK — VCP/p97 Apo (2.3Å)",
                "8OOI": "8OOI — VCP/p97 + CB-5083 (3.1Å)",
                "9DIL": "9DIL — VCP/p97 Cryo-EM (2.8Å)",
                "7K56": "7K56 — VCP/p97 D2 Domain (2.1Å)",
            }
            selected = st.selectbox(
                "Select Structure",
                options=list(pdb_options.keys()),
                format_func=lambda x: pdb_options[x],
                key="pdb_select",
            )
            pdb_id = selected

        with col2:
            residue = st.number_input("Highlight Residue", value=155, min_value=1, max_value=999, key="highlight_res")

        with col3:
            chain = st.selectbox("Chain", ["A", "B", "C", "D", "E", "F"], key="highlight_chain")
    else:
        residue = 155
        chain = "A"

    # Render Mol* viewer
    html = get_molstar_html(pdb_id, highlight_residue=residue, highlight_chain=chain)
    components.html(html, height=520, scrolling=False)

    # Structure details
    structure_info = {
        "5FTK": {"resolution": "2.3 Å", "method": "X-ray", "chains": 6, "ligand": "None (Apo)", "year": 2016},
        "8OOI": {"resolution": "3.1 Å", "method": "Cryo-EM", "chains": 6, "ligand": "CB-5083", "year": 2023},
        "9DIL": {"resolution": "2.8 Å", "method": "Cryo-EM", "chains": 6, "ligand": "Allosteric inhibitor", "year": 2024},
        "7K56": {"resolution": "2.1 Å", "method": "X-ray", "chains": 2, "ligand": "ATP analog", "year": 2021},
    }

    info = structure_info.get(pdb_id, {})
    if info:
        cols = st.columns(5)
        with cols[0]:
            st.metric("Resolution", info["resolution"])
        with cols[1]:
            st.metric("Method", info["method"])
        with cols[2]:
            st.metric("Chains", info["chains"])
        with cols[3]:
            st.metric("Ligand", info["ligand"])
        with cols[4]:
            st.metric("Year", info["year"])

    st.markdown("""
    **VCP/p97 (Valosin-containing protein)** — AAA+ ATPase critical for protein quality control.
    Mutations at R155 cause multisystem proteinopathy including frontotemporal dementia (FTD),
    inclusion body myopathy (IBM), and Paget disease of bone.

    🔴 **R155H mutation** highlighted in red — located in the N-D1 domain interface,
    disrupts cofactor binding and impairs hexamer dynamics.
    """)


if __name__ == "__main__":
    st.set_page_config(page_title="VCP Structure Viewer", layout="wide")
    render_protein_viewer()
