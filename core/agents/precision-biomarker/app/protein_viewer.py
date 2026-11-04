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
                # Verified against data.rcsb.org 2026-11-03. None of these entries contains
                # CB-5083; 8OOI was labelled as a CB-5083 co-structure and is an ADP.Pi state.
                "5FTK": "5FTK — p97 + ADP, cryo-EM (2.4Å)",
                "8OOI": "8OOI — p97 ADP·Pi state, cryo-EM (2.61Å)",
                "9DIL": "9DIL — p97 + VCPIP1 complex, cryo-EM (3.3Å)",
                "7K56": "7K56 — p97 dodecamer, cryo-EM (3.9Å)",
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
    # Every field here was wrong, including the method for two entries and a CB-5083 ligand
    # that none of these structures carries. Values below are from the RCSB entry records
    # (data.rcsb.org/rest/v1/core/entry/<ID>), checked 2026-11-03: method from exptl,
    # resolution from resolution_combined, chains from deposited_polymer_entity_instance_count,
    # ligands from nonpolymer_bound_components, year from the initial release date.
    structure_info = {
        "5FTK": {"resolution": "2.4 Å", "method": "Cryo-EM", "chains": 6, "ligand": "ADP", "year": 2016},
        "8OOI": {"resolution": "2.61 Å", "method": "Cryo-EM", "chains": 6, "ligand": "ADP, ATP, Mg, K, PO4", "year": 2024},
        "9DIL": {"resolution": "3.3 Å", "method": "Cryo-EM", "chains": 3, "ligand": "None (VCPIP1 complex)", "year": 2025},
        "7K56": {"resolution": "3.9 Å", "method": "Cryo-EM", "chains": 12, "ligand": "None reported", "year": 2021},
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
