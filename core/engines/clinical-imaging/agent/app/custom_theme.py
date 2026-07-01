"""Custom CSS theme for the Clinical Imaging Engine.

Provides a polished, professional dark theme with NVIDIA green
accents and medical-grade visual hierarchy.
"""

CUSTOM_CSS = """
<style>
/* ═══════════════════════════════════════════════════════════════
   CLINICAL IMAGING ENGINE — CUSTOM THEME
   NVIDIA Green: #76B900  |  Dark BG: #0E1117  |  Card BG: #1A1D23
   ═══════════════════════════════════════════════════════════════ */

/* ── Global Font & Spacing ── */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

html, body, [class*="st-"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
}

/* ── Hero Header ── */
.hero-banner {
    background: linear-gradient(135deg, #0a1628 0%, #122040 40%, #1a3a2a 100%);
    border: 1px solid rgba(118, 185, 0, 0.3);
    border-radius: 12px;
    padding: 28px 36px;
    margin-bottom: 24px;
    position: relative;
    overflow: hidden;
}

.hero-banner::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    height: 3px;
    background: linear-gradient(90deg, #76B900, #4CAF50, #76B900);
}

.hero-banner h1 {
    color: #ffffff;
    font-size: 2.2rem;
    font-weight: 700;
    margin: 0 0 6px 0;
    letter-spacing: -0.5px;
}

.hero-banner .subtitle {
    color: #76B900;
    font-size: 0.95rem;
    font-weight: 500;
    margin: 0 0 4px 0;
    letter-spacing: 0.5px;
}

.hero-banner .stats {
    color: rgba(255,255,255,0.6);
    font-size: 0.82rem;
    margin: 0;
}

.hero-banner .engine-badge {
    display: inline-block;
    background: rgba(118, 185, 0, 0.15);
    border: 1px solid rgba(118, 185, 0, 0.4);
    color: #76B900;
    padding: 3px 12px;
    border-radius: 20px;
    font-size: 0.75rem;
    font-weight: 600;
    letter-spacing: 1px;
    text-transform: uppercase;
    margin-bottom: 10px;
}

/* ── Clinical Disclaimer ── */
.clinical-disclaimer {
    background: rgba(255, 152, 0, 0.08);
    border: 1px solid rgba(255, 152, 0, 0.25);
    border-radius: 8px;
    padding: 10px 16px;
    margin-bottom: 16px;
    font-size: 0.78rem;
    color: rgba(255, 255, 255, 0.7);
    line-height: 1.4;
}

.clinical-disclaimer strong {
    color: #FFB74D;
}

/* ── Tab Styling ── */
.stTabs [data-baseweb="tab-list"] {
    gap: 4px;
    background: rgba(26, 29, 35, 0.5);
    border-radius: 10px;
    padding: 4px;
}

.stTabs [data-baseweb="tab"] {
    border-radius: 8px;
    padding: 8px 16px;
    font-size: 0.85rem;
    font-weight: 500;
    color: rgba(255, 255, 255, 0.65);
    transition: all 0.2s ease;
}

.stTabs [data-baseweb="tab"]:hover {
    color: rgba(255, 255, 255, 0.9);
    background: rgba(118, 185, 0, 0.1);
}

.stTabs [aria-selected="true"] {
    background: rgba(118, 185, 0, 0.15) !important;
    color: #76B900 !important;
    font-weight: 600;
    border-bottom: 2px solid #76B900 !important;
}

/* ── Metric Cards ── */
[data-testid="stMetric"] {
    background: linear-gradient(135deg, #1A1D23 0%, #1E2229 100%);
    border: 1px solid rgba(255, 255, 255, 0.08);
    border-radius: 10px;
    padding: 16px 18px;
    transition: border-color 0.2s ease;
}

[data-testid="stMetric"]:hover {
    border-color: rgba(118, 185, 0, 0.3);
}

[data-testid="stMetricLabel"] {
    color: rgba(255, 255, 255, 0.55) !important;
    font-size: 0.8rem !important;
    font-weight: 500 !important;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}

[data-testid="stMetricValue"] {
    color: #ffffff !important;
    font-size: 1.6rem !important;
    font-weight: 700 !important;
}

[data-testid="stMetricDelta"] > div {
    font-weight: 600 !important;
}

/* ── Sidebar ── */
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #0a0d12 0%, #111419 100%);
    border-right: 1px solid rgba(255, 255, 255, 0.06);
}

[data-testid="stSidebar"] .block-container {
    padding-top: 24px;
}

/* ── Sidebar Section Headers ── */
.sidebar-section {
    color: #76B900;
    font-size: 0.7rem;
    font-weight: 700;
    letter-spacing: 1.5px;
    text-transform: uppercase;
    margin: 20px 0 8px 0;
    padding-bottom: 4px;
    border-bottom: 1px solid rgba(118, 185, 0, 0.2);
}

/* ── Buttons ── */
.stButton > button[kind="primary"] {
    background: linear-gradient(135deg, #76B900, #5a9400);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    padding: 8px 24px;
    transition: all 0.2s ease;
}

.stButton > button[kind="primary"]:hover {
    background: linear-gradient(135deg, #8ad400, #76B900);
    box-shadow: 0 4px 12px rgba(118, 185, 0, 0.3);
}

.stButton > button:not([kind="primary"]) {
    border-radius: 8px;
    border: 1px solid rgba(255, 255, 255, 0.15);
    transition: all 0.2s ease;
}

.stButton > button:not([kind="primary"]):hover {
    border-color: rgba(118, 185, 0, 0.4);
    background: rgba(118, 185, 0, 0.08);
}

/* ── Expanders ── */
.streamlit-expanderHeader {
    font-weight: 600;
    border-radius: 8px;
}

/* ── Input Fields ── */
.stTextInput > div > div > input,
.stTextArea > div > div > textarea,
.stSelectbox > div > div > div {
    border-radius: 8px !important;
    border-color: rgba(255, 255, 255, 0.12) !important;
    transition: border-color 0.2s ease;
}

.stTextInput > div > div > input:focus,
.stTextArea > div > div > textarea:focus {
    border-color: #76B900 !important;
    box-shadow: 0 0 0 1px rgba(118, 185, 0, 0.3) !important;
}

/* ── Containers / Cards ── */
[data-testid="stExpander"] {
    border: 1px solid rgba(255, 255, 255, 0.08);
    border-radius: 10px;
    overflow: hidden;
}

div[data-testid="stVerticalBlock"] > div[style*="border"] {
    border-radius: 10px !important;
    border-color: rgba(255, 255, 255, 0.08) !important;
}

/* ── Charts ── */
.stBarChart, .stLineChart {
    border-radius: 10px;
    overflow: hidden;
}

/* ── DataFrames ── */
.stDataFrame {
    border-radius: 10px;
    overflow: hidden;
}

/* ── Dividers ── */
hr {
    border-color: rgba(255, 255, 255, 0.06) !important;
    margin: 20px 0 !important;
}

/* ── Success/Warning/Error/Info Alerts ── */
.stAlert {
    border-radius: 8px;
}

/* ── NIM Status Indicators ── */
.nim-status-card {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    background: rgba(26, 29, 35, 0.8);
    border: 1px solid rgba(255, 255, 255, 0.08);
    border-radius: 6px;
    padding: 4px 10px;
    font-size: 0.8rem;
}

/* ── Workflow Result Cards ── */
.workflow-result {
    background: linear-gradient(135deg, #1A1D23, #1E2229);
    border: 1px solid rgba(118, 185, 0, 0.2);
    border-radius: 12px;
    padding: 20px;
    margin: 12px 0;
}

.workflow-result .classification {
    font-size: 1.4rem;
    font-weight: 700;
    color: #76B900;
}

/* ── Severity Colors ── */
.severity-critical { color: #EF5350; font-weight: 700; }
.severity-urgent { color: #FF9800; font-weight: 700; }
.severity-significant { color: #FFC107; font-weight: 600; }
.severity-routine { color: #4CAF50; font-weight: 500; }
.severity-normal { color: #66BB6A; font-weight: 400; }

/* ── ACR Rating Bar ── */
.acr-bar {
    display: flex;
    gap: 2px;
    margin: 8px 0;
}

.acr-bar .segment {
    flex: 1;
    height: 8px;
    border-radius: 4px;
    opacity: 0.3;
}

.acr-bar .segment.active {
    opacity: 1;
}

.acr-bar .segment.green { background: #4CAF50; }
.acr-bar .segment.yellow { background: #FFC107; }
.acr-bar .segment.red { background: #EF5350; }

/* ── Footer ── */
.footer {
    text-align: center;
    padding: 24px 0;
    margin-top: 40px;
    border-top: 1px solid rgba(255, 255, 255, 0.06);
    color: rgba(255, 255, 255, 0.35);
    font-size: 0.75rem;
}

.footer a {
    color: #76B900;
    text-decoration: none;
}

/* ── Hide Streamlit Branding ── */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}

/* ── Scrollbar ── */
::-webkit-scrollbar {
    width: 6px;
    height: 6px;
}

::-webkit-scrollbar-track {
    background: transparent;
}

::-webkit-scrollbar-thumb {
    background: rgba(255, 255, 255, 0.15);
    border-radius: 3px;
}

::-webkit-scrollbar-thumb:hover {
    background: rgba(118, 185, 0, 0.4);
}
</style>
"""


def inject_theme():
    """Inject the custom CSS theme into the Streamlit page."""
    import streamlit as st
    st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


def render_hero_header():
    """Render the branded hero header banner."""
    import streamlit as st
    st.markdown("""
    <div class="hero-banner">
        <span class="engine-badge">Engine 4</span>
        <h1>Clinical Imaging Engine</h1>
        <p class="subtitle">HCLS AI Factory &mdash; Precision Medicine Platform</p>
        <p class="stats">9 Workflows &bull; 7 Scoring Systems &bull; 13 Collections &bull; 20 NVIDIA Technologies &bull; Apache 2.0</p>
    </div>
    """, unsafe_allow_html=True)


def render_clinical_disclaimer():
    """Render a styled clinical disclaimer."""
    import streamlit as st
    st.markdown("""
    <div class="clinical-disclaimer">
        <strong>Clinical Decision Support</strong> &mdash; This system provides evidence-based guidance for research and clinical decision support.
        All recommendations must be verified by a qualified healthcare professional. Not FDA-cleared.
    </div>
    """, unsafe_allow_html=True)


def render_footer():
    """Render the page footer."""
    import streamlit as st
    st.markdown("""
    <div class="footer">
        Clinical Imaging Engine v2.0 &bull; Part of the <a href="https://github.com/ajones1923/hcls-ai-factory">HCLS AI Factory</a> &bull; Apache 2.0 License<br>
        Patient DNA to Drug Candidates in &lt;5 hours on a single NVIDIA DGX Spark ($4,699)
    </div>
    """, unsafe_allow_html=True)
