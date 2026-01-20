"""
Drug Discovery Pipeline UI - Streamlit interface for VCP/FTD demo.

End-to-end flow: Target Hypothesis → Cryo-EM Evidence → Molecule Generation
"""
import streamlit as st
import json
import sys
from pathlib import Path
from typing import List, Dict, Any

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from target_import import TargetImporter, ImportedTarget, get_vcp_target
from cryoem_evidence import CryoEMEvidenceManager, CryoEMStructure, get_vcp_structures
from molecule_generator import MoleculeGenerator, GeneratedMolecule, generate_vcp_molecules
from structure_viewer import StructureViewer, HAS_STMOL

# Page config
st.set_page_config(
    page_title="Drug Discovery Pipeline | VCP/FTD Demo",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 1.5rem 2rem;
        border-radius: 12px;
        margin-bottom: 1.5rem;
        border: 1px solid #30363d;
        border-left: 4px solid #76b900;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
    }
    .main-header h1 {
        color: #ffffff;
        margin: 0;
        font-size: 1.8rem;
    }
    .main-header p {
        color: #a0aec0;
        margin: 0.5rem 0 0 0;
    }
    .pipeline-badge {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
        color: #000;
        padding: 0.3rem 0.9rem;
        border-radius: 20px;
        font-size: 0.7rem;
        font-weight: 700;
        display: inline-block;
        margin-left: 1rem;
        text-transform: uppercase;
    }
    .evidence-card {
        background: #f8f9fa;
        border: 1px solid #e1e4e8;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .structure-card {
        background: linear-gradient(135deg, #e8f5e9 0%, #ffffff 100%);
        border: 1px solid #4caf50;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .molecule-card {
        background: linear-gradient(135deg, #e3f2fd 0%, #ffffff 100%);
        border: 1px solid #2196f3;
        border-radius: 10px;
        padding: 1rem;
        margin: 0.5rem 0;
    }
    .quote-box {
        background: #fffde7;
        border-left: 4px solid #ffc107;
        padding: 1rem;
        border-radius: 0 8px 8px 0;
        font-style: italic;
        margin: 1rem 0;
    }
    .step-indicator {
        display: flex;
        align-items: center;
        gap: 0.5rem;
        padding: 0.5rem 1rem;
        background: #f0f0f0;
        border-radius: 20px;
        margin-bottom: 1rem;
    }
    .step-active {
        background: #76b900;
        color: white;
    }
</style>
""", unsafe_allow_html=True)


def render_header():
    """Render the main header."""
    st.markdown("""
    <div class="main-header">
        <h1>💊 Drug Discovery Pipeline <span class="pipeline-badge">VCP/FTD Demo</span></h1>
        <p>Target Hypothesis → Cryo-EM Structural Evidence → BioNeMo Molecule Generation</p>
    </div>
    """, unsafe_allow_html=True)


def render_sidebar():
    """Render the sidebar."""
    with st.sidebar:
        st.markdown("### 🧬 Pipeline Steps")

        step = st.radio(
            "Current Step",
            ["1. Target Hypothesis", "2. Cryo-EM Evidence", "3. Molecule Generation", "4. Summary"],
            label_visibility="collapsed"
        )

        st.markdown("---")

        st.markdown("### 🎯 Demo: Frontotemporal Dementia")
        st.info("""
        **Target:** VCP/p97 (Valosin-containing protein)

        **Disease:** Frontotemporal Dementia (FTD)

        **Variant:** rs188935092 (R159H region)

        **PDB Structures:** 8OOI, 9DIL, 7K56, 5FTK
        """)

        st.markdown("---")

        st.markdown("### 📚 Key Framing")
        st.markdown("""
        <div class="quote-box">
        "Genomics tells us <strong>what changed</strong>,
        Cryo-EM shows us <strong>how it changed</strong>,
        and generative AI helps us design <strong>what could fix it</strong>."
        </div>
        """, unsafe_allow_html=True)

        return step


def render_target_hypothesis():
    """Render the target hypothesis step."""
    st.markdown("## Step 1: Target Hypothesis")
    st.markdown("*From RAG Chat genomic evidence to validated drug target*")

    # Get VCP target
    target = get_vcp_target()

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### VCP/p97 - Validated Drug Target")

        st.markdown(f"""
        <div class="evidence-card">
            <h4>🧬 Gene: {target.gene}</h4>
            <p><strong>Protein:</strong> {target.protein}</p>
            <p><strong>Therapeutic Area:</strong> {target.therapeutic_area}</p>
            <p><strong>Mechanism:</strong> {target.mechanism}</p>
            <p><strong>Confidence:</strong> {target.confidence.upper()} | <strong>Priority:</strong> {target.priority}/5</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("#### Rationale")
        st.write(target.rationale)

        st.markdown("#### Key Variant Evidence")
        for v in target.variants:
            st.markdown(f"""
            - **{v.get('rsid', 'Unknown')}** at {v.get('chrom')}:{v.get('pos')}
              - Consequence: {v.get('consequence')}
              - Impact: {v.get('impact')}
              - Disease: {v.get('disease_association')}
            """)

    with col2:
        st.markdown("### Target Metrics")

        st.metric("Variant Count", target.variant_count)
        st.metric("Known PDB Structures", len(target.pdb_ids))
        st.metric("Priority Score", f"{target.priority}/5")

        st.markdown("### Available Structures")
        for pdb_id in target.pdb_ids:
            st.markdown(f"- [PDB:{pdb_id}](https://www.rcsb.org/structure/{pdb_id})")

    # Store in session state
    st.session_state.current_target = target


def render_cryoem_evidence():
    """Render the Cryo-EM evidence step."""
    st.markdown("## Step 2: Cryo-EM Structural Evidence")
    st.markdown("*Structural ground truth for drug design*")

    st.markdown("""
    <div class="quote-box">
    "This is where data becomes physics. We are grounding AI decisions in real molecular structures."
    </div>
    """, unsafe_allow_html=True)

    # Get VCP structures
    structures = get_vcp_structures()

    if not structures:
        st.warning("No Cryo-EM structures loaded. Check data/structures/vcp_structures.json")
        return

    # View mode toggle
    view_col1, view_col2 = st.columns([3, 1])
    with view_col1:
        st.markdown(f"### {len(structures)} VCP/p97 Structures Available")
    with view_col2:
        if HAS_STMOL:
            view_mode = st.selectbox(
                "View Mode",
                ["Interactive 3D", "Static Images"],
                key="structure_view_mode",
                help="Interactive 3D requires stmol/py3Dmol"
            )
        else:
            view_mode = "Static Images"
            st.caption("Install stmol for 3D view")

    # Initialize structure viewer if using 3D mode
    structure_viewer = StructureViewer() if view_mode == "Interactive 3D" else None

    for i, structure in enumerate(structures):
        with st.expander(f"📊 {structure.get_display_name()}", expanded=(i == 0)):
            col1, col2, col3 = st.columns([1.5, 1.5, 1])

            with col1:
                # Extract PDB ID for viewer
                pdb_id = structure.structure_id.replace("PDB:", "").strip()

                if view_mode == "Interactive 3D" and structure_viewer:
                    # Interactive 3D visualization
                    st.markdown("**Interactive 3D Structure:**")

                    # Style selector for this structure
                    style = st.selectbox(
                        "Visualization Style",
                        ["cartoon", "stick", "sphere", "line"],
                        key=f"style_{pdb_id}",
                        index=0
                    )

                    # Check if this structure has a known ligand/inhibitor
                    if structure.inhibitor:
                        show_ligand = st.checkbox(
                            f"Highlight binding site",
                            value=True,
                            key=f"ligand_{pdb_id}"
                        )
                        if show_ligand:
                            # Use ligand view for structures with inhibitors
                            structure_viewer.show_structure_with_ligand(
                                pdb_id,
                                ligand_resname="LIG",  # Common residue name
                                width=400,
                                height=350
                            )
                        else:
                            structure_viewer.show_structure_3d(
                                pdb_id,
                                style=style,
                                color_scheme="spectrum",
                                width=400,
                                height=350
                            )
                    else:
                        structure_viewer.show_structure_3d(
                            pdb_id,
                            style=style,
                            color_scheme="spectrum",
                            width=400,
                            height=350
                        )
                else:
                    # Static image display (original behavior)
                    st.markdown("**3D Structure:**")
                    if structure.image_url:
                        try:
                            st.image(
                                structure.image_url,
                                caption=f"{structure.structure_id} - {structure.conformation}",
                                use_container_width=True
                            )
                        except Exception as e:
                            st.warning(f"Could not load structure image")
                            # Fallback to PDB ID extraction for alternative URL
                            fallback_url = f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}/{pdb_id.lower()}_assembly-1.jpeg"
                            st.image(fallback_url, caption=structure.structure_id, use_container_width=True)

                # Display EMDB density map image if available (for Cryo-EM)
                if structure.emdb_image_url and "Cryo-EM" in structure.method:
                    st.markdown("**Cryo-EM Density Map:**")
                    try:
                        st.image(
                            structure.emdb_image_url,
                            caption=f"{structure.emdb_id} density map",
                            use_container_width=True
                        )
                    except Exception:
                        st.caption(f"EMDB: {structure.emdb_id}")

            with col2:
                st.markdown(f"""
                <div class="structure-card">
                    <h4>{structure.structure_id}</h4>
                    <p><strong>{structure.full_name or structure.protein}</strong></p>
                    <p><strong>Method:</strong> {structure.method}</p>
                    <p><strong>Resolution:</strong> {structure.resolution}</p>
                    <p><strong>Conformation:</strong> {structure.conformation}</p>
                </div>
                """, unsafe_allow_html=True)

                st.markdown("**Structural Context:**")
                st.write(structure.summary_text)

            with col3:
                st.markdown("**Binding Sites:**")
                for site in structure.binding_sites:
                    st.markdown(f"- {site}")

                st.markdown("**Druggable Pockets:**")
                for pocket in structure.druggable_pockets:
                    st.markdown(f"- {pocket}")

                if structure.inhibitor:
                    st.success(f"**Known Inhibitor:** {structure.inhibitor}")

                if structure.pdb_url:
                    st.link_button("View in RCSB PDB", structure.pdb_url)

                if structure.emdb_id:
                    emdb_url = f"https://www.ebi.ac.uk/emdb/entry/{structure.emdb_id}"
                    st.link_button("View in EMDB", emdb_url)

    # Best structure for drug design
    manager = CryoEMEvidenceManager()
    best = manager.get_best_structure_for_drug_design("VCP")
    if best:
        st.markdown("---")
        st.markdown("### 🎯 Recommended Structure for Drug Design")

        rec_col1, rec_col2 = st.columns([1, 2])
        with rec_col1:
            best_pdb_id = best.structure_id.replace("PDB:", "").strip()
            if view_mode == "Interactive 3D" and structure_viewer:
                # Interactive 3D for recommended structure
                if best.inhibitor:
                    structure_viewer.show_structure_with_ligand(
                        best_pdb_id,
                        ligand_resname="LIG",
                        width=350,
                        height=300
                    )
                else:
                    structure_viewer.show_structure_3d(
                        best_pdb_id,
                        style="cartoon",
                        color_scheme="spectrum",
                        width=350,
                        height=300
                    )
            elif best.image_url:
                st.image(best.image_url, caption=f"{best.structure_id}", use_container_width=True)
        with rec_col2:
            st.success(f"**{best.structure_id}** - {best.resolution} resolution")
            st.write(f"**Conformation:** {best.conformation}")
            if best.inhibitor:
                st.write(f"**Inhibitor:** {best.inhibitor}")
            st.write("This structure is optimal for structure-based drug design due to its high resolution and inhibitor-bound state.")

        st.session_state.best_structure = best


def render_molecule_generation():
    """Render the molecule generation step."""
    st.markdown("## Step 3: BioNeMo Molecule Generation")
    st.markdown("*AI-generated drug candidates from structural knowledge*")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Generation Parameters")

        num_molecules = st.slider("Number of molecules", 5, 20, 10)
        diversity = st.slider("Diversity", 0.1, 0.5, 0.3)

        use_seed = st.checkbox("Use known inhibitor as seed (CB-5083)", value=True)

        if st.button("🚀 Generate Molecules", type="primary"):
            with st.spinner("Generating VCP inhibitor candidates..."):
                molecules = generate_vcp_molecules(num_molecules)
                st.session_state.generated_molecules = molecules

    with col2:
        st.markdown("### Seed Molecule")
        st.info("""
        **CB-5083** - Known VCP inhibitor

        Used as seed for generating structural analogues
        that may have improved properties.
        """)

    # Display generated molecules
    if 'generated_molecules' in st.session_state and st.session_state.generated_molecules:
        molecules = st.session_state.generated_molecules

        st.markdown(f"### Generated {len(molecules)} Candidates")

        # Try to render molecules
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            rdkit_available = True
        except ImportError:
            rdkit_available = False

        for i, mol in enumerate(molecules):
            with st.expander(f"💊 {mol.name or f'Molecule-{i+1}'} (Score: {mol.score:.2f})", expanded=(i < 3)):
                col1, col2 = st.columns([1, 1])

                with col1:
                    st.markdown(f"""
                    <div class="molecule-card">
                        <h4>{mol.name or f'Molecule-{i+1}'}</h4>
                        <p><strong>Method:</strong> {mol.generation_method}</p>
                        <p><strong>Target:</strong> {mol.target_gene}</p>
                        <p><strong>Score:</strong> {mol.score:.2f}</p>
                    </div>
                    """, unsafe_allow_html=True)

                    st.code(mol.smiles, language="text")

                with col2:
                    st.markdown("**Properties:**")
                    props = mol.properties
                    st.write(f"- MW: {props.get('molecular_weight', 'N/A')}")
                    st.write(f"- LogP: {props.get('logP', 'N/A')}")
                    st.write(f"- HBD: {props.get('hbd', 'N/A')}")
                    st.write(f"- HBA: {props.get('hba', 'N/A')}")
                    st.write(f"- Lipinski violations: {props.get('lipinski_violations', 'N/A')}")

                    if rdkit_available:
                        try:
                            rdkit_mol = Chem.MolFromSmiles(mol.smiles)
                            if rdkit_mol:
                                img = Draw.MolToImage(rdkit_mol, size=(250, 200))
                                st.image(img, caption=mol.name)
                        except Exception as e:
                            st.warning(f"Could not render molecule: {e}")


def render_summary():
    """Render the summary step."""
    st.markdown("## Step 4: Pipeline Summary")
    st.markdown("*End-to-end genomics to drug discovery*")

    st.markdown("""
    <div class="quote-box">
    "This completes the loop from genome → structure → molecule.
    What traditionally took months can now be explored in minutes."
    </div>
    """, unsafe_allow_html=True)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### 🧬 Genomics")
        target = st.session_state.get('current_target')
        if target:
            st.success(f"Target: **{target.gene}**")
            st.write(f"- {target.variant_count} variants analyzed")
            st.write(f"- Disease: Frontotemporal Dementia")
            st.write(f"- Confidence: {target.confidence}")

    with col2:
        st.markdown("### 📊 Cryo-EM")
        structure = st.session_state.get('best_structure')
        if structure:
            st.success(f"Structure: **{structure.structure_id}**")
            st.write(f"- Resolution: {structure.resolution}")
            st.write(f"- {len(structure.druggable_pockets)} druggable pockets")
            st.write(f"- Known inhibitor: {structure.inhibitor or 'CB-5083'}")

    with col3:
        st.markdown("### 💊 Molecules")
        molecules = st.session_state.get('generated_molecules', [])
        if molecules:
            st.success(f"Generated: **{len(molecules)} candidates**")
            st.write(f"- Best score: {max(m.score for m in molecules):.2f}")
            st.write(f"- Method: {molecules[0].generation_method}")
            st.write("- Ready for docking/screening")

    st.markdown("---")

    st.markdown("### 🎯 Demo Narrative")
    st.markdown("""
    1. **From the RAG Chat**, we identified VCP as a high-priority target based on the rs188935092 variant
       linked to frontotemporal dementia.

    2. **Cryo-EM structural evidence** grounds our understanding in physical reality - we can see the
       actual 3D shape of VCP, its binding pockets, and how inhibitors like CB-5083 interact with it.

    3. **BioNeMo molecule generation** produces novel candidate molecules informed by both the
       genomic context and structural knowledge - molecules designed to potentially modulate VCP
       and address the underlying disease mechanism.

    This is the promise of AI-accelerated drug discovery: **from variant to molecule in one session**.
    """)

    # Export options
    st.markdown("### 📤 Export Results")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("Export Molecules (JSON)"):
            molecules = st.session_state.get('generated_molecules', [])
            if molecules:
                generator = MoleculeGenerator()
                output_file = generator.save_molecules(molecules, "vcp_ftd_candidates.json")
                st.success(f"Saved to {output_file}")

    with col2:
        if st.button("Export Full Report (PDF)"):
            with st.spinner("Generating professional PDF report..."):
                try:
                    import sys
                    sys.path.insert(0, str(Path(__file__).parent.parent))
                    from generate_vcp_report_enhanced import VCPReportGeneratorEnhanced

                    generator = VCPReportGeneratorEnhanced(
                        output_path="outputs/VCP_Drug_Candidate_Report.pdf"
                    )
                    output_file = generator.generate()

                    # Read the PDF for download
                    with open(output_file, "rb") as f:
                        pdf_data = f.read()

                    st.success(f"Report generated: {output_file}")
                    st.download_button(
                        label="Download PDF Report",
                        data=pdf_data,
                        file_name="VCP_Drug_Candidate_Report.pdf",
                        mime="application/pdf"
                    )
                except Exception as e:
                    st.error(f"Report generation failed: {e}")


def main():
    """Main application."""
    render_header()
    step = render_sidebar()

    if step == "1. Target Hypothesis":
        render_target_hypothesis()
    elif step == "2. Cryo-EM Evidence":
        render_cryoem_evidence()
    elif step == "3. Molecule Generation":
        render_molecule_generation()
    elif step == "4. Summary":
        render_summary()


if __name__ == "__main__":
    main()
