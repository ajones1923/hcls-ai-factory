"""
Drug Discovery Pipeline UI - Streamlit interface for target-driven drug discovery.

Supports any target gene from Stage 2 RAG/Chat export, with VCP/FTD as default demo.
End-to-end flow: Target Hypothesis → Structural Evidence → Molecule Generation
"""
import streamlit as st
import json
import os
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from target_import import TargetImporter, ImportedTarget, get_vcp_target
from cryoem_evidence import CryoEMEvidenceManager, CryoEMStructure
from molecule_generator import MoleculeGenerator, GeneratedMolecule
from structure_viewer import StructureViewer, HAS_STMOL

# Default paths for Stage 2 export
RAG_EXPORT_PATH = Path(os.environ.get(
    "RAG_EXPORT_PATH",
    str(Path(__file__).parent.parent.parent / "rag-chat-pipeline" / "data" / "targets" / "targets_for_phase5.json")
))


def load_targets_from_export() -> List[Dict[str, Any]]:
    """Load targets from Stage 2 RAG/Chat export file."""
    if RAG_EXPORT_PATH.exists():
        with open(RAG_EXPORT_PATH) as f:
            data = json.load(f)
        return data.get('targets', [])
    return []


def get_active_target() -> ImportedTarget:
    """Get the currently active target from session state or load from export."""
    if 'current_target' in st.session_state and st.session_state.current_target:
        return st.session_state.current_target

    # Try loading from Stage 2 export
    exported = load_targets_from_export()
    if exported:
        t = exported[0]
        target = ImportedTarget(
            gene=t.get('gene', 'Unknown'),
            protein=t.get('protein'),
            rationale=t.get('rationale', ''),
            confidence=t.get('confidence', 'medium'),
            priority=t.get('priority', 3),
            therapeutic_area=t.get('therapeutic_area'),
            mechanism=t.get('mechanism'),
            pdb_ids=t.get('pdb_ids', []),
            variant_count=t.get('variant_count', 0),
            variants=t.get('variants', []),
            source_query=t.get('source_query'),
            status=t.get('status', 'imported'),
        )
        st.session_state.current_target = target
        ref_compounds = t.get('reference_compounds', [])
        if ref_compounds:
            st.session_state.reference_smiles = ref_compounds[0].get('smiles')
            st.session_state.reference_drug = ref_compounds[0].get('name')
        else:
            st.session_state.reference_smiles = t.get('reference_smiles')
            st.session_state.reference_drug = t.get('reference_drug')
        return target

    # Fall back to VCP demo
    target = get_vcp_target()
    st.session_state.current_target = target
    st.session_state.reference_smiles = "CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5"
    st.session_state.reference_drug = "CB-5083"
    return target


# Page config
st.set_page_config(
    page_title="Drug Discovery Pipeline",
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
    target = get_active_target()
    disease = target.therapeutic_area or "Drug Discovery"
    badge = f"{target.gene}"

    st.markdown(f"""
    <div class="main-header">
        <h1>💊 Drug Discovery Pipeline <span class="pipeline-badge">{badge}</span></h1>
        <p>Target Hypothesis → Structural Evidence → BioNeMo Molecule Generation</p>
    </div>
    """, unsafe_allow_html=True)


def render_sidebar():
    """Render the sidebar."""
    target = get_active_target()

    with st.sidebar:
        st.markdown("### 🧬 Pipeline Steps")

        step = st.radio(
            "Current Step",
            ["1. Target Hypothesis", "2. Structural Evidence", "3. Molecule Generation", "4. Summary"],
            label_visibility="collapsed"
        )

        st.markdown("---")

        # Target selection
        exported = load_targets_from_export()
        if len(exported) > 1:
            st.markdown("### 🎯 Select Target")
            gene_options = [t.get('gene', 'Unknown') for t in exported]
            selected_gene = st.selectbox("Target Gene", gene_options)
            if selected_gene != target.gene:
                for t in exported:
                    if t.get('gene') == selected_gene:
                        new_target = ImportedTarget(
                            gene=t.get('gene', 'Unknown'),
                            protein=t.get('protein'),
                            rationale=t.get('rationale', ''),
                            confidence=t.get('confidence', 'medium'),
                            priority=t.get('priority', 3),
                            therapeutic_area=t.get('therapeutic_area'),
                            mechanism=t.get('mechanism'),
                            pdb_ids=t.get('pdb_ids', []),
                            variant_count=t.get('variant_count', 0),
                            variants=t.get('variants', []),
                            source_query=t.get('source_query'),
                            status=t.get('status', 'imported'),
                        )
                        st.session_state.current_target = new_target
                        st.session_state.reference_smiles = t.get('reference_smiles')
                        st.session_state.reference_drug = t.get('reference_drug')
                        st.session_state.pop('generated_molecules', None)
                        st.session_state.pop('best_structure', None)
                        st.rerun()
            st.markdown("---")

        # Current target info
        st.markdown(f"### 🎯 Current Target")
        diseases = ', '.join(target.variants[0].get('disease_association', '').split('/')) if target.variants and target.variants[0].get('disease_association') else (target.therapeutic_area or 'N/A')
        st.info(f"""
        **Target:** {target.gene} ({target.protein or 'Unknown protein'})

        **Area:** {target.therapeutic_area or 'N/A'}

        **PDB Structures:** {', '.join(target.pdb_ids) if target.pdb_ids else 'None'}

        **Confidence:** {target.confidence}
        """)

        st.markdown("---")

        st.markdown("### 📚 Key Framing")
        st.markdown("""
        <div class="quote-box">
        "Genomics tells us <strong>what changed</strong>,
        structural biology shows us <strong>how it changed</strong>,
        and generative AI helps us design <strong>what could fix it</strong>."
        </div>
        """, unsafe_allow_html=True)

        return step


def render_target_hypothesis():
    """Render the target hypothesis step."""
    st.markdown("## Step 1: Target Hypothesis")
    st.markdown("*From RAG Chat genomic evidence to validated drug target*")

    target = get_active_target()

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown(f"### {target.gene} - Validated Drug Target")

        st.markdown(f"""
        <div class="evidence-card">
            <h4>🧬 Gene: {target.gene}</h4>
            <p><strong>Protein:</strong> {target.protein or 'Unknown'}</p>
            <p><strong>Therapeutic Area:</strong> {target.therapeutic_area or 'N/A'}</p>
            <p><strong>Mechanism:</strong> {target.mechanism or 'N/A'}</p>
            <p><strong>Confidence:</strong> {target.confidence.upper()} | <strong>Priority:</strong> {target.priority}/5</p>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("#### Rationale")
        st.write(target.rationale)

        if target.variants:
            st.markdown("#### Key Variant Evidence")
            for v in target.variants:
                if isinstance(v, dict):
                    st.markdown(f"""
                    - **{v.get('rsid', 'Unknown')}** at {v.get('chrom', '?')}:{v.get('pos', '?')}
                      - Consequence: {v.get('consequence', 'N/A')}
                      - Impact: {v.get('impact', 'N/A')}
                      - Disease: {v.get('disease_association', 'N/A')}
                    """)

    with col2:
        st.markdown("### Target Metrics")

        st.metric("Variant Count", target.variant_count)
        st.metric("Known PDB Structures", len(target.pdb_ids))
        st.metric("Priority Score", f"{target.priority}/5")

        if target.pdb_ids:
            st.markdown("### Available Structures")
            for pdb_id in target.pdb_ids:
                st.markdown(f"- [PDB:{pdb_id}](https://www.rcsb.org/structure/{pdb_id})")

        ref_drug = st.session_state.get('reference_drug')
        if ref_drug:
            st.markdown("### Seed Compound")
            st.success(f"**{ref_drug}**")


def render_structural_evidence():
    """Render the structural evidence step."""
    target = get_active_target()

    st.markdown("## Step 2: Structural Evidence")
    st.markdown("*Structural ground truth for drug design*")

    st.markdown("""
    <div class="quote-box">
    "This is where data becomes physics. We are grounding AI decisions in real molecular structures."
    </div>
    """, unsafe_allow_html=True)

    # Load structures: cached first, then auto-fetch from RCSB if PDB IDs available
    manager = CryoEMEvidenceManager()
    if target.pdb_ids:
        structures = manager.auto_load_structures(target.gene, target.pdb_ids)
    else:
        structures = manager.load_structures_for_gene(target.gene)

    # View mode toggle
    view_col1, view_col2 = st.columns([3, 1])
    with view_col1:
        if structures:
            st.markdown(f"### {len(structures)} {target.gene} Structures Available")
        elif target.pdb_ids:
            st.markdown(f"### {len(target.pdb_ids)} PDB Structures for {target.gene}")
        else:
            st.warning(f"No structures available for {target.gene}")
            return

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

    if structures:
        # Show cached structures with full metadata
        for i, structure in enumerate(structures):
            with st.expander(f"📊 {structure.get_display_name()}", expanded=(i == 0)):
                col1, col2, col3 = st.columns([1.5, 1.5, 1])

                with col1:
                    pdb_id = structure.structure_id.replace("PDB:", "").strip()
                    if view_mode == "Interactive 3D" and structure_viewer:
                        st.markdown("**Interactive 3D Structure:**")
                        style = st.selectbox("Style", ["cartoon", "stick", "sphere", "line"], key=f"style_{pdb_id}", index=0)
                        if structure.inhibitor:
                            show_ligand = st.checkbox(f"Highlight binding site", value=True, key=f"ligand_{pdb_id}")
                            if show_ligand:
                                structure_viewer.show_structure_with_ligand(pdb_id, ligand_resname="LIG", width=400, height=350)
                            else:
                                structure_viewer.show_structure_3d(pdb_id, style=style, color_scheme="spectrum", width=400, height=350)
                        else:
                            structure_viewer.show_structure_3d(pdb_id, style=style, color_scheme="spectrum", width=400, height=350)
                    else:
                        st.markdown("**3D Structure:**")
                        fallback_url = f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}/{pdb_id.lower()}_assembly-1.jpeg"
                        try:
                            if structure.image_url:
                                st.image(structure.image_url, caption=f"{structure.structure_id} - {structure.conformation}", use_container_width=True)
                            else:
                                st.image(fallback_url, caption=structure.structure_id, use_container_width=True)
                        except Exception:
                            st.caption(f"PDB: {pdb_id}")

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
                    if structure.binding_sites:
                        st.markdown("**Binding Sites:**")
                        for site in structure.binding_sites:
                            st.markdown(f"- {site}")
                    if structure.druggable_pockets:
                        st.markdown("**Druggable Pockets:**")
                        for pocket in structure.druggable_pockets:
                            st.markdown(f"- {pocket}")
                    if structure.inhibitor:
                        st.success(f"**Known Inhibitor:** {structure.inhibitor}")
                    if structure.pdb_url:
                        st.link_button("View in RCSB PDB", structure.pdb_url)

        # Best structure recommendation
        best = manager.get_best_structure_for_drug_design(target.gene)
        if best:
            st.session_state.best_structure = best

    elif target.pdb_ids:
        # No cached structures — show PDB IDs with 3D viewer
        st.info(f"No pre-cached structure data for {target.gene}. Showing PDB structures directly.")
        for i, pdb_id in enumerate(target.pdb_ids):
            with st.expander(f"📊 PDB: {pdb_id}", expanded=(i == 0)):
                col1, col2 = st.columns([2, 1])
                with col1:
                    if view_mode == "Interactive 3D" and structure_viewer:
                        structure_viewer.show_structure_3d(pdb_id, style="cartoon", color_scheme="spectrum", width=500, height=400)
                    else:
                        img_url = f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}/{pdb_id.lower()}_assembly-1.jpeg"
                        try:
                            st.image(img_url, caption=f"PDB: {pdb_id}", use_container_width=True)
                        except Exception:
                            st.caption(f"PDB: {pdb_id}")
                with col2:
                    st.link_button("View in RCSB PDB", f"https://www.rcsb.org/structure/{pdb_id}")


def render_molecule_generation():
    """Render the molecule generation step."""
    target = get_active_target()
    ref_smiles = st.session_state.get('reference_smiles', '')
    ref_drug = st.session_state.get('reference_drug', '')

    st.markdown("## Step 3: BioNeMo Molecule Generation")
    st.markdown("*AI-generated drug candidates from structural knowledge*")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Generation Parameters")

        num_molecules = st.slider("Number of molecules", 5, 20, 10)
        diversity = st.slider("Diversity", 0.1, 0.5, 0.3)

        seed_label = f"Use known compound as seed ({ref_drug})" if ref_drug else "Use seed compound"
        use_seed = st.checkbox(seed_label, value=bool(ref_smiles))

        if st.button("🚀 Generate Molecules", type="primary"):
            seed = ref_smiles if use_seed and ref_smiles else None
            with st.spinner(f"Generating {target.gene} inhibitor candidates..."):
                generator = MoleculeGenerator()
                if seed:
                    molecules = generator.generate_from_seed(
                        seed_smiles=seed,
                        target_gene=target.gene,
                        num_molecules=num_molecules,
                        diversity=diversity,
                    )
                else:
                    molecules = generator.generate_from_seed(
                        seed_smiles="c1ccccc1",  # Benzene fallback
                        target_gene=target.gene,
                        num_molecules=num_molecules,
                        diversity=diversity,
                    )
                st.session_state.generated_molecules = molecules

    with col2:
        st.markdown("### Seed Molecule")
        if ref_drug and ref_smiles:
            st.info(f"""
            **{ref_drug}** — Known {target.gene} compound

            Used as seed for generating structural analogues
            that may have improved properties.
            """)
        else:
            st.warning("No known seed compound available for this target. Generic scaffolds will be used.")

    # Display generated molecules
    if 'generated_molecules' in st.session_state and st.session_state.generated_molecules:
        molecules = st.session_state.generated_molecules

        st.markdown(f"### Generated {len(molecules)} Candidates")

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
    target = get_active_target()
    ref_drug = st.session_state.get('reference_drug', '')

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
        st.success(f"Target: **{target.gene}**")
        st.write(f"- {target.variant_count} variants analyzed")
        st.write(f"- Area: {target.therapeutic_area or 'N/A'}")
        st.write(f"- Confidence: {target.confidence}")

    with col2:
        st.markdown("### 📊 Structures")
        structure = st.session_state.get('best_structure')
        if structure:
            st.success(f"Structure: **{structure.structure_id}**")
            st.write(f"- Resolution: {structure.resolution}")
            st.write(f"- {len(structure.druggable_pockets)} druggable pockets")
            if structure.inhibitor:
                st.write(f"- Known inhibitor: {structure.inhibitor}")
        elif target.pdb_ids:
            st.success(f"Structures: **{len(target.pdb_ids)} PDB entries**")
            for pdb_id in target.pdb_ids[:3]:
                st.write(f"- {pdb_id}")
        else:
            st.warning("No structures available")

    with col3:
        st.markdown("### 💊 Molecules")
        molecules = st.session_state.get('generated_molecules', [])
        if molecules:
            st.success(f"Generated: **{len(molecules)} candidates**")
            st.write(f"- Best score: {max(m.score for m in molecules):.2f}")
            st.write(f"- Method: {molecules[0].generation_method}")
            st.write("- Ready for docking/screening")

    st.markdown("---")

    st.markdown("### 🎯 Pipeline Narrative")
    st.markdown(f"""
    1. **From the RAG Chat**, we identified **{target.gene}** as a high-priority target
       based on genomic variant analysis{f' linked to {target.therapeutic_area}' if target.therapeutic_area else ''}.

    2. **Structural evidence** grounds our understanding in physical reality — we can see the
       actual 3D shape of {target.gene}, its binding pockets, and how compounds interact with it.

    3. **BioNeMo molecule generation** produces novel candidate molecules informed by both the
       genomic context and structural knowledge — molecules designed to potentially modulate {target.gene}
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
                filename = f"{target.gene.lower()}_candidates.json"
                output_file = generator.save_molecules(molecules, filename)
                st.success(f"Saved to {output_file}")

    with col2:
        if st.button("Export Full Report (PDF)"):
            with st.spinner("Generating professional PDF report..."):
                try:
                    sys.path.insert(0, str(Path(__file__).parent.parent))
                    from generate_vcp_report_enhanced import VCPReportGeneratorEnhanced

                    output_path = f"outputs/{target.gene}_Drug_Candidate_Report.pdf"
                    generator = VCPReportGeneratorEnhanced(output_path=output_path)
                    output_file = generator.generate()

                    with open(output_file, "rb") as f:
                        pdf_data = f.read()

                    st.success(f"Report generated: {output_file}")
                    st.download_button(
                        label="Download PDF Report",
                        data=pdf_data,
                        file_name=f"{target.gene}_Drug_Candidate_Report.pdf",
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
    elif step == "2. Structural Evidence":
        render_structural_evidence()
    elif step == "3. Molecule Generation":
        render_molecule_generation()
    elif step == "4. Summary":
        render_summary()


if __name__ == "__main__":
    main()
