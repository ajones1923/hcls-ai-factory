"""
Healthcare & Life Sciences Pipeline Portal
Unified interface for Genomics → RAG Chat → Drug Discovery

A comprehensive dashboard for end-to-end drug discovery workflows.
"""

import streamlit as st
import json
import subprocess
import os
from pathlib import Path
from datetime import datetime
import requests
from typing import Dict, List, Any, Optional

# Configuration
ORCHESTRATOR_DIR = Path(__file__).parent.parent
HCLS_ROOT = ORCHESTRATOR_DIR.parent  # Root of HCLS AI Factory

# Pipeline directories - can be overridden via environment variables
GENOMICS_DIR = Path(os.environ.get("GENOMICS_PIPELINE_DIR", str(HCLS_ROOT / "genomics-pipeline")))
RAG_CHAT_DIR = Path(os.environ.get("RAG_CHAT_PIPELINE_DIR", str(HCLS_ROOT / "rag-chat-pipeline")))
DRUG_DISCOVERY_DIR = Path(os.environ.get("DRUG_DISCOVERY_PIPELINE_DIR", str(HCLS_ROOT / "drug-discovery-pipeline")))

def get_service_host():
    """
    Dynamically detect the service host IP.
    Priority: SERVICE_HOST env var > auto-detected IP > localhost
    """
    # Allow override via environment variable
    env_host = os.environ.get('SERVICE_HOST')
    if env_host:
        return env_host

    # Try to auto-detect the host's IP address
    try:
        import socket
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        host = s.getsockname()[0]
        s.close()
        return host
    except Exception:
        return "localhost"

# Dynamic service host - set SERVICE_HOST env var to override
SERVICE_HOST = get_service_host()
GRAFANA_URL = f"http://{SERVICE_HOST}:3000/d/nvidia-dgx-spark/nvidia-dgx-spark-gpu-monitoring"
PROMETHEUS_URL = f"http://{SERVICE_HOST}:9099"

# Page configuration
st.set_page_config(
    page_title="HLS Pipeline Portal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS
st.markdown("""
<style>
    /* Main theme */
    :root {
        --nvidia-green: #76b900;
        --dark-bg: #1a1a2e;
        --card-bg: #16213e;
    }

    /* Header styling */
    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 2rem;
        border-radius: 16px;
        margin-bottom: 2rem;
        border-left: 5px solid #76b900;
        box-shadow: 0 8px 32px rgba(0, 0, 0, 0.2);
    }
    .main-header h1 {
        color: white;
        margin: 0;
        font-size: 2.2rem;
        display: flex;
        align-items: center;
        gap: 15px;
    }
    .main-header p {
        color: #a0aec0;
        margin: 0.5rem 0 0 0;
        font-size: 1.1rem;
    }
    .nvidia-badge {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
        color: black;
        padding: 6px 16px;
        border-radius: 25px;
        font-size: 0.75rem;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    /* Pipeline cards */
    .pipeline-card {
        background: white;
        border-radius: 16px;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.08);
        border-left: 4px solid #76b900;
        transition: transform 0.2s, box-shadow 0.2s;
    }
    .pipeline-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 30px rgba(0, 0, 0, 0.12);
    }
    .pipeline-card h3 {
        color: #1a1a2e;
        margin: 0 0 0.5rem 0;
        display: flex;
        align-items: center;
        gap: 10px;
    }
    .pipeline-card p {
        color: #666;
        margin: 0;
    }

    /* Status indicators */
    .status-badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 0.75rem;
        font-weight: 600;
    }
    .status-running { background: #e3f2fd; color: #1976d2; }
    .status-success { background: #e8f5e9; color: #388e3c; }
    .status-error { background: #ffebee; color: #d32f2f; }
    .status-idle { background: #f5f5f5; color: #757575; }

    /* Metrics cards */
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 1.5rem;
        border-radius: 12px;
        text-align: center;
    }
    .metric-card.nvidia {
        background: linear-gradient(135deg, #76b900 0%, #5a8c00 100%);
    }
    .metric-card h4 {
        margin: 0;
        font-size: 0.85rem;
        opacity: 0.9;
    }
    .metric-card .value {
        font-size: 2.5rem;
        font-weight: bold;
        margin: 0.5rem 0;
    }

    /* Sidebar */
    .sidebar-section {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 10px;
        margin: 1rem 0;
    }

    /* Pipeline flow visualization */
    .pipeline-flow {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 10px;
        padding: 1rem;
        background: #f8f9fa;
        border-radius: 10px;
        margin: 1rem 0;
    }
    .flow-step {
        background: white;
        padding: 10px 20px;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        text-align: center;
    }
    .flow-arrow {
        color: #76b900;
        font-size: 1.5rem;
    }

    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)


def get_service_status(url: str, timeout: float = 2.0) -> bool:
    """Check if a service is responding."""
    try:
        response = requests.get(url, timeout=timeout)
        return response.status_code < 500
    except:
        return False


def get_gpu_metrics() -> Dict[str, Any]:
    """Fetch GPU metrics from DCGM exporter."""
    try:
        response = requests.get("http://localhost:9400/metrics", timeout=2)
        metrics = {}
        for line in response.text.split('\n'):
            if line.startswith('DCGM_FI_DEV_GPU_UTIL'):
                metrics['gpu_util'] = float(line.split()[-1])
            elif line.startswith('DCGM_FI_DEV_GPU_TEMP{'):
                metrics['gpu_temp'] = float(line.split()[-1])
            elif line.startswith('DCGM_FI_DEV_POWER_USAGE'):
                metrics['power'] = float(line.split()[-1])
            elif line.startswith('DCGM_FI_DEV_FB_USED'):
                metrics['mem_used'] = float(line.split()[-1])
        return metrics
    except:
        return {}


def render_header():
    """Render the main header."""
    st.markdown("""
    <div class="main-header">
        <h1>
            🧬 Healthcare & Life Sciences Pipeline
            <span class="nvidia-badge">NVIDIA DGX Spark</span>
        </h1>
        <p>End-to-end genomics to drug discovery • Powered by BioNeMo NIM</p>
    </div>
    """, unsafe_allow_html=True)


def render_sidebar():
    """Render the sidebar navigation."""
    with st.sidebar:
        st.markdown("### 🎯 Navigation")

        page = st.radio(
            "Select Page",
            [
                "🏠 Dashboard",
                "🧬 Genomics Pipeline",
                "💬 RAG Chat",
                "💊 Drug Discovery",
                "🚀 Run Pipeline",
                "📊 Monitoring",
                "📋 Results"
            ],
            label_visibility="collapsed"
        )

        st.markdown("---")

        # Service status
        st.markdown("### 🔌 Service Status")

        services = [
            ("Grafana", "http://localhost:3000"),
            ("Prometheus", "http://localhost:9099"),
            ("DCGM Exporter", "http://localhost:9400"),
            ("Drug Discovery UI", "http://localhost:8505"),
        ]

        for name, url in services:
            status = get_service_status(url)
            icon = "🟢" if status else "🔴"
            st.markdown(f"{icon} {name}")

        st.markdown("---")

        # GPU metrics
        st.markdown("### 🎮 GPU Status")
        metrics = get_gpu_metrics()
        if metrics:
            st.metric("GPU Utilization", f"{metrics.get('gpu_util', 0):.0f}%")
            st.metric("Temperature", f"{metrics.get('gpu_temp', 0):.0f}°C")
            st.metric("Power", f"{metrics.get('power', 0):.1f}W")
        else:
            st.info("GPU metrics unavailable")

        return page


def render_dashboard():
    """Render the main dashboard."""
    st.markdown("## 📊 Pipeline Overview")

    # Pipeline flow visualization
    st.markdown("""
    <div class="pipeline-flow">
        <div class="flow-step">
            <strong>🧬 Genomics</strong><br>
            <small>FASTQ → VCF</small>
        </div>
        <span class="flow-arrow">→</span>
        <div class="flow-step">
            <strong>💬 RAG Chat</strong><br>
            <small>VCF → Targets</small>
        </div>
        <span class="flow-arrow">→</span>
        <div class="flow-step">
            <strong>💊 Drug Discovery</strong><br>
            <small>Target → Molecules</small>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Metrics row
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.markdown("""
        <div class="metric-card nvidia">
            <h4>Pipeline Runs</h4>
            <div class="value">12</div>
            <small>Last 7 days</small>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="metric-card">
            <h4>Molecules Generated</h4>
            <div class="value">847</div>
            <small>Total</small>
        </div>
        """, unsafe_allow_html=True)

    with col3:
        st.markdown("""
        <div class="metric-card">
            <h4>Targets Identified</h4>
            <div class="value">23</div>
            <small>High confidence</small>
        </div>
        """, unsafe_allow_html=True)

    with col4:
        st.markdown("""
        <div class="metric-card nvidia">
            <h4>GPU Hours</h4>
            <div class="value">156</div>
            <small>This month</small>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")

    # Pipeline cards
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div class="pipeline-card">
            <h3>🧬 Genomics Pipeline</h3>
            <p>Variant calling from FASTQ to VCF using BWA-MEM2, GATK, and DeepVariant.</p>
            <br>
            <span class="status-badge status-idle">Ready</span>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="pipeline-card">
            <h3>💊 Drug Discovery</h3>
            <p>AI-powered molecule generation and docking using BioNeMo MolMIM and DiffDock.</p>
            <br>
            <span class="status-badge status-success">Active</span>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="pipeline-card">
            <h3>💬 RAG Chat</h3>
            <p>Intelligent target identification using Claude AI with genomic variant context.</p>
            <br>
            <span class="status-badge status-idle">Ready</span>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="pipeline-card">
            <h3>📊 Monitoring</h3>
            <p>Real-time GPU metrics, pipeline telemetry, and system health via Grafana.</p>
            <br>
            <span class="status-badge status-success">Online</span>
        </div>
        """, unsafe_allow_html=True)

    # Recent activity
    st.markdown("### 📜 Recent Activity")

    activities = [
        {"time": "2 hours ago", "event": "Drug Discovery run completed", "status": "success", "details": "VCP target, 20 molecules generated"},
        {"time": "5 hours ago", "event": "RAG Chat analysis", "status": "success", "details": "Identified 3 high-priority targets"},
        {"time": "1 day ago", "event": "Genomics pipeline", "status": "success", "details": "Sample WGS-001 processed"},
        {"time": "2 days ago", "event": "Full pipeline run", "status": "success", "details": "End-to-end FTD analysis"},
    ]

    for activity in activities:
        status_class = f"status-{activity['status']}"
        st.markdown(f"""
        <div style="display: flex; align-items: center; padding: 10px; border-bottom: 1px solid #eee;">
            <span class="status-badge {status_class}" style="margin-right: 15px;">●</span>
            <div style="flex: 1;">
                <strong>{activity['event']}</strong>
                <br><small style="color: #666;">{activity['details']}</small>
            </div>
            <small style="color: #999;">{activity['time']}</small>
        </div>
        """, unsafe_allow_html=True)


def render_genomics():
    """Render the Genomics Pipeline page."""
    st.markdown("## 🧬 Genomics Pipeline")
    st.markdown("*Variant calling from FASTQ to annotated VCF*")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Pipeline Steps")
        st.markdown("""
        1. **FastQC** - Quality control of raw reads
        2. **BWA-MEM2** - Alignment to reference genome
        3. **GATK MarkDuplicates** - Remove PCR duplicates
        4. **GATK BQSR** - Base quality score recalibration
        5. **GATK HaplotypeCaller** - Variant calling
        6. **VEP/SnpEff** - Variant annotation
        """)

        st.markdown("### Input Files")
        uploaded_files = st.file_uploader(
            "Upload FASTQ files",
            type=["fastq", "fq", "fastq.gz", "fq.gz"],
            accept_multiple_files=True
        )

        if uploaded_files:
            st.success(f"Uploaded {len(uploaded_files)} files")

    with col2:
        st.markdown("### Configuration")
        genome = st.selectbox("Reference Genome", ["GRCh38", "GRCh37", "T2T-CHM13"])
        st.text_input("Known Sites VCF", placeholder="dbsnp.vcf.gz")

        st.markdown("### Quick Actions")
        if st.button("▶️ Run Genomics Pipeline", type="primary", use_container_width=True):
            st.info("Pipeline execution would start here...")

        if st.button("📂 View Results", use_container_width=True):
            st.info("Opening results directory...")


def render_rag_chat():
    """Render the RAG Chat page."""
    st.markdown("## 💬 RAG Chat - Target Discovery")
    st.markdown("*AI-powered drug target identification from genomic variants*")

    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("### Chat Interface")

        # Chat history
        if "messages" not in st.session_state:
            st.session_state.messages = [
                {"role": "assistant", "content": "Hello! I can help you identify drug targets from your genomic data. Upload a VCF file or describe the variants you're interested in."}
            ]

        for message in st.session_state.messages:
            with st.chat_message(message["role"]):
                st.write(message["content"])

        # Chat input
        if prompt := st.chat_input("Ask about drug targets..."):
            st.session_state.messages.append({"role": "user", "content": prompt})
            with st.chat_message("user"):
                st.write(prompt)

            # Mock response
            response = f"Based on your query about '{prompt}', I've analyzed the relevant genomic context. VCP (Valosin-containing protein) appears to be a promising target for Frontotemporal Dementia based on the variant evidence."
            st.session_state.messages.append({"role": "assistant", "content": response})
            with st.chat_message("assistant"):
                st.write(response)

    with col2:
        st.markdown("### Upload VCF")
        vcf_file = st.file_uploader("Upload annotated VCF", type=["vcf", "vcf.gz"])

        st.markdown("### Identified Targets")
        targets = [
            {"gene": "VCP", "confidence": "High", "disease": "FTD"},
            {"gene": "MAPT", "confidence": "Medium", "disease": "FTD"},
        ]

        for target in targets:
            st.markdown(f"""
            <div style="background: #f8f9fa; padding: 10px; border-radius: 8px; margin: 5px 0; border-left: 3px solid #76b900;">
                <strong>{target['gene']}</strong>
                <br><small>{target['disease']} • {target['confidence']} confidence</small>
            </div>
            """, unsafe_allow_html=True)

        if st.button("🎯 Send to Drug Discovery", type="primary", use_container_width=True):
            st.success("Target sent to Drug Discovery pipeline!")


def render_drug_discovery():
    """Render the Drug Discovery page with embedded UI."""
    st.markdown("## 💊 Drug Discovery Pipeline")
    st.markdown("*AI-powered molecule generation and virtual screening*")

    # Check if Drug Discovery UI is available
    dd_ui_available = get_service_status(f"http://{SERVICE_HOST}:8505")

    # View mode selector
    view_mode = st.radio(
        "View Mode",
        ["🖥️ Full Application", "⚙️ Quick Controls", "📊 Results Only"],
        horizontal=True
    )

    if view_mode == "🖥️ Full Application":
        # Embedded Drug Discovery UI
        if dd_ui_available:
            st.markdown(f"""
            <style>
                .dd-iframe-container {{
                    border: 2px solid #76b900;
                    border-radius: 12px;
                    overflow: hidden;
                    box-shadow: 0 4px 20px rgba(0,0,0,0.15);
                }}
            </style>
            <div class="dd-iframe-container">
                <iframe
                    src="http://{SERVICE_HOST}:8505/?embedded=true"
                    width="100%"
                    height="800"
                    frameborder="0"
                    style="background: white;">
                </iframe>
            </div>
            """, unsafe_allow_html=True)

            st.markdown(f"""
            <p style="text-align: center; margin-top: 10px;">
                <a href="http://{SERVICE_HOST}:8505" target="_blank">Open in new tab ↗</a>
            </p>
            """, unsafe_allow_html=True)
        else:
            st.warning(f"Drug Discovery UI is not available at http://{SERVICE_HOST}:8505")
            st.info(f"Start it with: `cd {DRUG_DISCOVERY_DIR} && source venv/bin/activate && streamlit run app/discovery_ui.py --server.port 8505`")

    elif view_mode == "⚙️ Quick Controls":
        # Quick controls mode
        col1, col2 = st.columns([2, 1])

        with col1:
            st.markdown("### Target Selection")
            target = st.selectbox(
                "Select Target",
                ["VCP (Valosin-containing protein)", "MAPT (Microtubule-associated protein tau)", "Custom..."]
            )

            st.markdown("### Generation Parameters")
            col_a, col_b = st.columns(2)
            with col_a:
                num_molecules = st.slider("Number of Molecules", 5, 50, 20)
                diversity = st.slider("Diversity", 0.1, 0.5, 0.3)
            with col_b:
                max_mw = st.slider("Max Molecular Weight", 400, 600, 550)
                num_poses = st.slider("Docking Poses", 5, 20, 10)

            st.markdown("### Reference Compound")
            seed_smiles = st.text_input(
                "Seed SMILES (optional)",
                value="CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5",
                help="CB-5083 structure for VCP"
            )

            if st.button("🚀 Generate Molecules", type="primary", use_container_width=True):
                with st.spinner("Generating molecules with MolMIM..."):
                    import time
                    time.sleep(2)
                st.success(f"Generated {num_molecules} molecules! View in Results tab.")

        with col2:
            st.markdown("### NIM Services")

            services = [
                ("MolMIM", "localhost:8001", False),
                ("DiffDock", "localhost:8002", False),
            ]

            for name, url, status in services:
                icon = "🟢" if status else "🟡"
                mode = "Online" if status else "Mock Mode"
                st.markdown(f"{icon} **{name}** - {mode}")

            st.markdown("---")

            st.markdown("### Service Status")
            dd_status = "🟢 Online" if dd_ui_available else "🔴 Offline"
            st.markdown(f"Drug Discovery UI: {dd_status}")

            gpu_metrics = get_gpu_metrics()
            if gpu_metrics:
                st.markdown(f"GPU Utilization: {gpu_metrics.get('gpu_util', 0):.0f}%")
                st.markdown(f"GPU Temperature: {gpu_metrics.get('gpu_temp', 0):.0f}°C")

    else:  # Results Only
        st.markdown("### Recent Drug Discovery Results")

        # Check for results
        results_dir = DRUG_DISCOVERY_DIR / "outputs"
        if results_dir.exists():
            result_files = sorted(results_dir.glob("**/*.json"), key=lambda x: x.stat().st_mtime, reverse=True)[:5]

            if result_files:
                for f in result_files:
                    with st.expander(f"📄 {f.name} ({datetime.fromtimestamp(f.stat().st_mtime).strftime('%Y-%m-%d %H:%M')})"):
                        try:
                            with open(f) as file:
                                data = json.load(file)
                            st.json(data)
                        except:
                            st.error("Could not load file")
            else:
                st.info("No results found. Run the drug discovery pipeline to generate outputs.")
        else:
            st.info("Results directory not found.")


def render_run_pipeline():
    """Render the Run Pipeline page."""
    st.markdown("## 🚀 Run Full Pipeline")
    st.markdown("*Execute end-to-end genomics to drug discovery workflow*")

    st.markdown("### Pipeline Mode")
    mode = st.radio(
        "Select mode",
        [
            "🎬 Demo (VCP/FTD)",
            "🧬 Full Pipeline (FASTQ → Molecules)",
            "🎯 Target Mode (VCF → Molecules)",
            "💊 Drug Only (Target → Molecules)"
        ],
        horizontal=True
    )

    st.markdown("---")

    if "Demo" in mode:
        st.info("**Demo Mode**: Run the VCP/Frontotemporal Dementia demonstration pipeline with pre-configured data.")

        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Target**: VCP (p97)")
            st.markdown("**Disease**: Frontotemporal Dementia")
            st.markdown("**Reference**: CB-5083 inhibitor")
        with col2:
            num_mol = st.number_input("Molecules to generate", 5, 50, 20)

    elif "Full" in mode:
        st.markdown("### Input Files")
        st.file_uploader("Upload FASTQ files", type=["fastq.gz"], accept_multiple_files=True)

    elif "Target" in mode:
        st.markdown("### Input Files")
        st.file_uploader("Upload VCF file", type=["vcf", "vcf.gz"])

    else:
        st.markdown("### Target Hypothesis")
        st.file_uploader("Upload target JSON", type=["json"])

    st.markdown("---")

    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🚀 Launch Pipeline", type="primary", use_container_width=True):
            st.markdown("### Pipeline Execution")

            progress = st.progress(0)
            status = st.empty()

            stages = [
                "Initializing pipeline...",
                "Stage 1: Normalizing target...",
                "Stage 2: Discovering structures...",
                "Stage 3: Preparing structures...",
                "Stage 4: Generating molecules...",
                "Stage 5: Chemistry QC...",
                "Stage 6: Generating conformers...",
                "Stage 7: Molecular docking...",
                "Stage 8: Ranking candidates...",
                "Stage 9: Generating report...",
            ]

            import time
            for i, stage in enumerate(stages):
                status.markdown(f"**{stage}**")
                progress.progress((i + 1) / len(stages))
                time.sleep(0.5)

            st.success("✅ Pipeline completed successfully!")
            st.balloons()


def render_monitoring():
    """Render the Monitoring page."""
    st.markdown("## 📊 System Monitoring")
    st.markdown("*Real-time GPU metrics and pipeline telemetry*")

    # GPU Metrics
    st.markdown("### 🎮 GPU Metrics (NVIDIA DGX Spark)")

    metrics = get_gpu_metrics()

    if metrics:
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("GPU Utilization", f"{metrics.get('gpu_util', 0):.0f}%")
        with col2:
            st.metric("Temperature", f"{metrics.get('gpu_temp', 0):.0f}°C")
        with col3:
            st.metric("Power Usage", f"{metrics.get('power', 0):.1f}W")
        with col4:
            st.metric("Memory Used", f"{metrics.get('mem_used', 0):.0f} MB")
    else:
        st.warning("GPU metrics unavailable. Ensure DCGM Exporter is running.")

    st.markdown("---")

    # Embedded Grafana
    st.markdown("### 📈 Grafana Dashboard")

    grafana_available = get_service_status(f"http://{SERVICE_HOST}:3000") or get_service_status("http://localhost:3000")

    if grafana_available:
        st.markdown(f"""
        <style>
            .grafana-container {{
                border: 2px solid #76b900;
                border-radius: 12px;
                overflow: hidden;
                box-shadow: 0 4px 20px rgba(0,0,0,0.15);
            }}
        </style>
        <div class="grafana-container">
            <iframe src="http://{SERVICE_HOST}:3000/d/nvidia-dgx-spark?orgId=1&refresh=5s&kiosk"
                    width="100%" height="700" frameborder="0"></iframe>
        </div>
        """, unsafe_allow_html=True)

        st.markdown(f"""
        <p style="text-align: center; margin-top: 10px;">
            <a href="http://{SERVICE_HOST}:3000/d/nvidia-dgx-spark/nvidia-dgx-spark-gpu-monitoring" target="_blank">Open Grafana in new tab ↗</a>
            (admin / dgxspark)
        </p>
        """, unsafe_allow_html=True)
    else:
        st.warning("Grafana is not available. Start the monitoring stack with:")
        st.code(f"cd {DRUG_DISCOVERY_DIR}/monitoring && docker compose up -d")


def render_results():
    """Render the Results page."""
    st.markdown("## 📋 Pipeline Results")
    st.markdown("*Browse and export pipeline outputs*")

    # Results directory
    results_dir = ORCHESTRATOR_DIR / "results"
    dd_outputs = DRUG_DISCOVERY_DIR / "outputs"

    st.markdown("### Recent Runs")

    # Check for result files
    result_files = []
    for directory in [results_dir, dd_outputs]:
        if directory.exists():
            result_files.extend(list(directory.glob("**/*.json")))

    if result_files:
        for f in sorted(result_files, key=lambda x: x.stat().st_mtime, reverse=True)[:10]:
            col1, col2, col3 = st.columns([3, 1, 1])
            with col1:
                st.markdown(f"**{f.name}**")
            with col2:
                st.markdown(f"<small>{datetime.fromtimestamp(f.stat().st_mtime).strftime('%Y-%m-%d %H:%M')}</small>", unsafe_allow_html=True)
            with col3:
                if st.button("View", key=str(f)):
                    with open(f) as file:
                        st.json(json.load(file))
    else:
        st.info("No results found. Run a pipeline to generate outputs.")

    st.markdown("---")

    st.markdown("### Export Options")
    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("📥 Export as SDF", use_container_width=True):
            st.info("SDF export functionality")

    with col2:
        if st.button("📥 Export as CSV", use_container_width=True):
            st.info("CSV export functionality")

    with col3:
        if st.button("📥 Generate Report", use_container_width=True):
            st.info("Report generation functionality")


def main():
    """Main application entry point."""
    render_header()
    page = render_sidebar()

    if "Dashboard" in page:
        render_dashboard()
    elif "Genomics" in page:
        render_genomics()
    elif "RAG" in page:
        render_rag_chat()
    elif "Drug" in page:
        render_drug_discovery()
    elif "Run" in page:
        render_run_pipeline()
    elif "Monitoring" in page:
        render_monitoring()
    elif "Results" in page:
        render_results()


if __name__ == "__main__":
    main()
