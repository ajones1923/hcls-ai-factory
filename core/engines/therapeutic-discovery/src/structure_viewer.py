"""
3D Structure Viewer for Drug Discovery Pipeline
Uses stmol/py3Dmol for interactive molecular visualization in Streamlit
"""

import streamlit as st
from pathlib import Path
import requests
from typing import Optional
import os

# Try to import 3D visualization libraries
try:
    import py3Dmol
    from stmol import showmol
    HAS_STMOL = True
except ImportError:
    HAS_STMOL = False
    print("Note: Install stmol for 3D visualization: pip install stmol py3Dmol")


class StructureViewer:
    """
    Handles fetching, caching, and displaying protein structures.
    """

    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the structure viewer.

        Args:
            cache_dir: Directory for caching PDB files. Defaults to data/structures/pdb_cache
        """
        if cache_dir is None:
            cache_dir = Path(__file__).parent.parent / "data" / "structures" / "pdb_cache"
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_pdb(self, pdb_id: str) -> Optional[str]:
        """
        Fetch PDB file content, using cache if available.

        Args:
            pdb_id: 4-character PDB ID (e.g., "8OOI")

        Returns:
            PDB file content as string, or None if fetch failed
        """
        pdb_id = pdb_id.upper().replace("PDB:", "")
        cache_file = self.cache_dir / f"{pdb_id.lower()}.pdb"

        # Check cache first
        if cache_file.exists():
            return cache_file.read_text()

        # Fetch from RCSB
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            pdb_content = response.text

            # Cache it
            cache_file.write_text(pdb_content)
            return pdb_content
        except Exception as e:
            print(f"Failed to fetch PDB {pdb_id}: {e}")
            return None

    def get_cif(self, pdb_id: str) -> Optional[str]:
        """
        Fetch mmCIF file content (alternative to PDB format).

        Args:
            pdb_id: 4-character PDB ID

        Returns:
            mmCIF file content as string, or None if fetch failed
        """
        pdb_id = pdb_id.upper().replace("PDB:", "")
        cache_file = self.cache_dir / f"{pdb_id.lower()}.cif"

        if cache_file.exists():
            return cache_file.read_text()

        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            cif_content = response.text
            cache_file.write_text(cif_content)
            return cif_content
        except Exception as e:
            print(f"Failed to fetch CIF {pdb_id}: {e}")
            return None

    def show_structure_3d(
        self,
        pdb_id: str,
        style: str = "cartoon",
        color_scheme: str = "spectrum",
        width: int = 700,
        height: int = 500,
        highlight_residues: Optional[list] = None,
        show_surface: bool = False,
        background_color: str = "white",
    ):
        """
        Display an interactive 3D structure in Streamlit.

        Args:
            pdb_id: PDB ID to display
            style: Visualization style - "cartoon", "stick", "sphere", "line"
            color_scheme: Color scheme - "spectrum", "chain", "residue", "element"
            width: Viewer width in pixels
            height: Viewer height in pixels
            highlight_residues: List of residue numbers to highlight
            show_surface: Whether to show molecular surface
            background_color: Background color
        """
        if not HAS_STMOL:
            st.warning("3D visualization requires stmol: `pip install stmol py3Dmol`")
            # Fall back to showing a link
            st.markdown(f"[View {pdb_id} on RCSB](https://www.rcsb.org/3d-view/{pdb_id})")
            return

        pdb_content = self.get_pdb(pdb_id)
        if pdb_content is None:
            st.error(f"Could not load structure {pdb_id}")
            return

        # Create the 3Dmol viewer
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(pdb_content, "pdb")

        # Apply style
        style_spec = {}
        if style == "cartoon":
            style_spec = {"cartoon": {"color": color_scheme}}
        elif style == "stick":
            style_spec = {"stick": {"colorscheme": color_scheme}}
        elif style == "sphere":
            style_spec = {"sphere": {"colorscheme": color_scheme}}
        elif style == "line":
            style_spec = {"line": {"colorscheme": color_scheme}}

        viewer.setStyle(style_spec)

        # Highlight specific residues if requested
        if highlight_residues:
            for resi in highlight_residues:
                viewer.addStyle(
                    {"resi": resi},
                    {"stick": {"color": "red"}, "sphere": {"color": "red", "radius": 0.5}}
                )

        # Add surface if requested
        if show_surface:
            viewer.addSurface(py3Dmol.VDW, {"opacity": 0.7, "color": "white"})

        # Set background and zoom
        viewer.setBackgroundColor(background_color)
        viewer.zoomTo()

        # Render in Streamlit
        showmol(viewer, height=height, width=width)

    def show_structure_with_ligand(
        self,
        pdb_id: str,
        ligand_resname: str = "LIG",
        width: int = 700,
        height: int = 500,
    ):
        """
        Display structure with highlighted ligand binding site.

        Args:
            pdb_id: PDB ID
            ligand_resname: Residue name of the ligand
            width: Viewer width
            height: Viewer height
        """
        if not HAS_STMOL:
            st.warning("3D visualization requires stmol: `pip install stmol py3Dmol`")
            return

        pdb_content = self.get_pdb(pdb_id)
        if pdb_content is None:
            st.error(f"Could not load structure {pdb_id}")
            return

        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(pdb_content, "pdb")

        # Show protein as cartoon
        viewer.setStyle({"cartoon": {"color": "spectrum"}})

        # Highlight ligand as sticks
        viewer.addStyle(
            {"resn": ligand_resname},
            {"stick": {"colorscheme": "greenCarbon", "radius": 0.3}}
        )

        # Show binding site residues (within 5Å of ligand)
        viewer.addStyle(
            {"byres": True, "within": {"distance": 5, "sel": {"resn": ligand_resname}}},
            {"stick": {"colorscheme": "whiteCarbon", "radius": 0.15}}
        )

        viewer.setBackgroundColor("white")
        viewer.zoomTo({"resn": ligand_resname})

        showmol(viewer, height=height, width=width)


def render_vcp_structure_gallery():
    """
    Render an interactive gallery of VCP structures for the Drug Discovery UI.
    """
    st.markdown("### 🔬 VCP Protein Structures")

    viewer = StructureViewer()

    structures = [
        {
            "pdb_id": "8OOI",
            "title": "VCP Hexamer (Wild-type)",
            "description": "Reference structure for FTD mutations",
        },
        {
            "pdb_id": "5FTK",
            "title": "VCP + CB-5083 Inhibitor",
            "description": "Drug binding site template",
            "ligand": "5N4",  # CB-5083 residue name in PDB
        },
        {
            "pdb_id": "7K56",
            "title": "VCP D2 Domain",
            "description": "ATP-bound active state",
        },
    ]

    # Create tabs for each structure
    tabs = st.tabs([s["title"] for s in structures])

    for tab, struct in zip(tabs, structures):
        with tab:
            col1, col2 = st.columns([2, 1])

            with col1:
                if "ligand" in struct:
                    viewer.show_structure_with_ligand(
                        struct["pdb_id"],
                        ligand_resname=struct["ligand"],
                        height=400
                    )
                else:
                    viewer.show_structure_3d(
                        struct["pdb_id"],
                        style="cartoon",
                        color_scheme="spectrum",
                        height=400
                    )

            with col2:
                st.markdown(f"**{struct['title']}**")
                st.markdown(struct["description"])
                st.markdown(f"[View on RCSB](https://www.rcsb.org/structure/{struct['pdb_id']})")

                # Style controls
                st.selectbox(
                    "Style",
                    ["Cartoon", "Stick", "Surface"],
                    key=f"style_{struct['pdb_id']}"
                )


# Example usage in Streamlit
if __name__ == "__main__":
    st.set_page_config(page_title="Structure Viewer", layout="wide")
    render_vcp_structure_gallery()
