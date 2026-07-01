"""Dose Intelligence tab -- radiation dose tracking and optimization.

Calls the FastAPI backend at /dose/* for dose recording, cumulative
tracking, DRL comparison, and population analytics. Also shows static
dose reference comparison from local data.
"""

import json
from datetime import datetime
from pathlib import Path

import streamlit as st
import requests

API_BASE = "http://localhost:8524"
DOSE_DATA_PATH = Path(__file__).parent.parent.parent / "data" / "reference" / "dose_reference.json"


def _load_dose_data():
    """Load the static dose reference comparison data."""
    try:
        with open(DOSE_DATA_PATH) as f:
            return json.load(f)
    except Exception:
        return []


def _alert_color(level: str) -> str:
    """Return a CSS color for dose alert levels."""
    return {
        "normal": "#2e7d32",
        "elevated": "#f9a825",
        "high": "#e65100",
        "critical": "#c62828",
    }.get(level, "#666")


def _alert_icon(level: str) -> str:
    """Return a status indicator for dose alert levels."""
    return {
        "normal": "[OK]",
        "elevated": "[!]",
        "high": "[!!]",
        "critical": "[!!!]",
    }.get(level, "[-]")


def _drl_status_color(status: str) -> str:
    """Return a CSS color for DRL comparison status."""
    return {
        "below_achievable": "#2e7d32",
        "below_drl": "#4caf50",
        "above_drl": "#e65100",
        "significantly_above": "#c62828",
        "no_drl_available": "#666",
    }.get(status, "#666")


def _drl_status_label(status: str) -> str:
    """Return a human-readable label for DRL status."""
    return {
        "below_achievable": "Below Achievable Level (Excellent)",
        "below_drl": "Below DRL (Acceptable)",
        "above_drl": "Above DRL (Review Needed)",
        "significantly_above": "Significantly Above DRL (Urgent Review)",
        "no_drl_available": "No DRL Available",
    }.get(status, status)


def render():
    st.header("Dose Intelligence")
    st.caption(
        "Radiation dose tracking, DRL comparison, and AI-driven optimization "
        "following ALARA principles"
    )

    tab_record, tab_cumulative, tab_drl, tab_reference, tab_pediatric = st.tabs([
        "Record Dose",
        "Cumulative Dose",
        "DRL Comparison",
        "Reference Data",
        "Pediatric Considerations",
    ])

    # ══════════════════════════════════════════════════════════════════
    # TAB 1: Record Dose
    # ══════════════════════════════════════════════════════════════════
    with tab_record:
        st.subheader("Record a Dose")
        st.markdown("Log a radiation dose for a patient study into the tracking registry.")

        # Generate Demo Data button
        col_demo, col_space = st.columns([1, 2])
        with col_demo:
            if st.button("Generate Demo Data (200 records)", type="secondary"):
                try:
                    resp = requests.get(f"{API_BASE}/dose/population", timeout=10)
                    if resp.status_code == 200:
                        data = resp.json()
                        st.success(
                            f"Demo data loaded: {data.get('total_records', 0)} records, "
                            f"{data.get('unique_patients', 0)} patients"
                        )
                    else:
                        st.error(f"Error: {resp.text}")
                except requests.ConnectionError:
                    st.error(f"Cannot connect to API at {API_BASE}.")
                except Exception as e:
                    st.error(f"Error: {e}")

        st.divider()

        with st.form("dose_record_form"):
            col1, col2 = st.columns(2)

            with col1:
                patient_id = st.text_input(
                    "Patient ID", value="PAT-00001",
                    placeholder="e.g. PAT-00001",
                )
                study_date = st.date_input("Study Date", value=datetime.now())
                modality = st.selectbox(
                    "Modality",
                    ["CT", "XR", "NM", "MAMMO", "FLUORO"],
                )
                protocol = st.text_input(
                    "Protocol",
                    placeholder="e.g. CT Head without contrast",
                )

            with col2:
                body_region = st.selectbox(
                    "Body Region",
                    ["head", "chest", "abdomen", "pelvis", "spine", "cardiac", "breast", "neck", "extremity"],
                )
                effective_dose = st.number_input(
                    "Effective Dose (mSv)", min_value=0.0, max_value=200.0,
                    value=5.0, step=0.1,
                )
                dlp = st.number_input(
                    "DLP (mGy-cm)", min_value=0.0, value=0.0, step=1.0,
                    help="Dose-Length Product (CT only)",
                )
                ctdi = st.number_input(
                    "CTDIvol (mGy)", min_value=0.0, value=0.0, step=0.1,
                    help="CT Dose Index, Volume (CT only)",
                )

            pediatric = st.checkbox("Pediatric patient")

            submitted = st.form_submit_button("Record Dose", type="primary", use_container_width=True)

            if submitted:
                if not protocol:
                    st.warning("Please enter a protocol name.")
                else:
                    payload = {
                        "patient_id": patient_id,
                        "study_date": study_date.strftime("%Y-%m-%d"),
                        "modality": modality,
                        "protocol": protocol,
                        "body_region": body_region,
                        "effective_dose_msv": effective_dose,
                        "pediatric": pediatric,
                    }
                    if dlp > 0:
                        payload["dlp_mgy_cm"] = dlp
                    if ctdi > 0:
                        payload["ctdi_vol_mgy"] = ctdi

                    try:
                        resp = requests.post(
                            f"{API_BASE}/dose/record",
                            json=payload,
                            timeout=10,
                        )
                        if resp.status_code == 200:
                            result = resp.json()
                            st.success(
                                f"Dose recorded: {result.get('protocol')} -- "
                                f"{result.get('effective_dose_msv')} mSv | "
                                f"Registry size: {result.get('registry_size')}"
                            )
                        else:
                            st.error(f"Error ({resp.status_code}): {resp.text}")
                    except requests.ConnectionError:
                        st.error(f"Cannot connect to API at {API_BASE}.")
                    except Exception as e:
                        st.error(f"Error recording dose: {e}")

    # ══════════════════════════════════════════════════════════════════
    # TAB 2: Cumulative Dose
    # ══════════════════════════════════════════════════════════════════
    with tab_cumulative:
        st.subheader("Cumulative Dose Lookup")
        st.markdown("View cumulative radiation exposure for a patient over time.")

        col_pid, col_period = st.columns([2, 1])
        with col_pid:
            lookup_id = st.text_input(
                "Patient ID", value="PAT-00001", key="cum_pid",
                placeholder="Enter patient ID",
            )
        with col_period:
            period = st.number_input(
                "Lookback period (days)", min_value=30, max_value=3650,
                value=365, step=30, key="cum_period",
            )

        if st.button("Look Up Cumulative Dose", type="primary", use_container_width=True):
            try:
                resp = requests.get(
                    f"{API_BASE}/dose/cumulative/{lookup_id}",
                    params={"period_days": period},
                    timeout=10,
                )
                if resp.status_code != 200:
                    st.error(f"Error ({resp.status_code}): {resp.text}")
                else:
                    data = resp.json()
                    _display_cumulative_dose(data)
            except requests.ConnectionError:
                st.error(f"Cannot connect to API at {API_BASE}.")
            except Exception as e:
                st.error(f"Error: {e}")

    # ══════════════════════════════════════════════════════════════════
    # TAB 3: DRL Comparison
    # ══════════════════════════════════════════════════════════════════
    with tab_drl:
        st.subheader("Diagnostic Reference Level Comparison")
        st.markdown(
            "Compare a patient dose to national Diagnostic Reference Levels (DRLs) "
            "to identify optimization opportunities."
        )

        with st.form("drl_form"):
            col_a, col_b = st.columns(2)
            with col_a:
                drl_protocol = st.selectbox(
                    "Protocol",
                    [
                        "CT Head without contrast",
                        "CTA Head/Neck",
                        "CT Chest",
                        "Low-dose CT Chest",
                        "CT Abdomen/Pelvis with contrast",
                        "CT Abdomen/Pelvis without contrast",
                        "CT Abdomen multiphase",
                        "CTA Coronary",
                        "CT Chest PE protocol",
                        "CT Spine Lumbar",
                        "CT Spine Cervical",
                    ],
                )
                drl_dose = st.number_input(
                    "Patient Dose (mSv)", min_value=0.1, max_value=200.0,
                    value=8.0, step=0.5,
                )
            with col_b:
                drl_modality = st.selectbox("Modality", ["CT", "XR", "NM"], key="drl_mod")
                drl_region = st.selectbox(
                    "Body Region",
                    ["head", "chest", "abdomen", "spine", "cardiac"],
                    key="drl_region",
                )
                drl_pediatric = st.checkbox("Pediatric", key="drl_pedi")

            drl_submit = st.form_submit_button("Compare to DRL", type="primary", use_container_width=True)

        if drl_submit:
            payload = {
                "protocol": drl_protocol,
                "effective_dose_msv": drl_dose,
                "modality": drl_modality,
                "body_region": drl_region,
                "pediatric": drl_pediatric,
            }
            try:
                resp = requests.post(
                    f"{API_BASE}/dose/compare-drl",
                    json=payload,
                    timeout=10,
                )
                if resp.status_code != 200:
                    st.error(f"Error ({resp.status_code}): {resp.text}")
                else:
                    _display_drl_comparison(resp.json())
            except requests.ConnectionError:
                st.error(f"Cannot connect to API at {API_BASE}.")
            except Exception as e:
                st.error(f"Error: {e}")

    # ══════════════════════════════════════════════════════════════════
    # TAB 4: Reference Data (static dose comparison chart)
    # ══════════════════════════════════════════════════════════════════
    with tab_reference:
        _render_reference_data()

    # ══════════════════════════════════════════════════════════════════
    # TAB 5: Pediatric Considerations
    # ══════════════════════════════════════════════════════════════════
    with tab_pediatric:
        _render_pediatric_section()


# ══════════════════════════════════════════════════════════════════════
# DISPLAY HELPERS
# ══════════════════════════════════════════════════════════════════════


def _display_cumulative_dose(data: dict):
    """Display cumulative dose results with alert visualization."""
    total = data.get("total_effective_dose_msv", 0)
    count = data.get("study_count", 0)
    alert_level = data.get("alert_level", "normal")
    alert_msg = data.get("alert_message")
    date_range = data.get("date_range", {})

    # Alert level banner
    color = _alert_color(alert_level)
    icon = _alert_icon(alert_level)
    st.markdown(
        f'<div style="padding:12px;border-radius:8px;border:2px solid {color};'
        f'background:{color}15;margin-bottom:16px;">'
        f'<span style="font-size:1.2em;font-weight:bold;color:{color};">'
        f"{icon} Alert Level: {alert_level.upper()}</span></div>",
        unsafe_allow_html=True,
    )

    if alert_msg:
        if alert_level in ("critical", "high"):
            st.error(alert_msg)
        elif alert_level == "elevated":
            st.warning(alert_msg)
        else:
            st.info("Cumulative dose within normal range.")

    # Metrics
    m1, m2, m3 = st.columns(3)
    with m1:
        st.metric("Total Effective Dose", f"{total:.1f} mSv")
    with m2:
        st.metric("Studies", count)
    with m3:
        first = date_range.get("first", "N/A")
        last = date_range.get("last", "N/A")
        st.metric("Date Range", f"{first} to {last}" if first != "N/A" else "No data")

    # Breakdown by modality and body region
    by_modality = data.get("by_modality", {})
    by_region = data.get("by_body_region", {})

    if by_modality or by_region:
        col1, col2 = st.columns(2)
        with col1:
            if by_modality:
                st.markdown("**Dose by Modality**")
                import pandas as pd
                df_mod = pd.DataFrame(
                    [{"Modality": k, "Dose (mSv)": v} for k, v in sorted(by_modality.items(), key=lambda x: -x[1])]
                )
                st.dataframe(df_mod, use_container_width=True, hide_index=True)

        with col2:
            if by_region:
                st.markdown("**Dose by Body Region**")
                import pandas as pd
                df_reg = pd.DataFrame(
                    [{"Region": k, "Dose (mSv)": v} for k, v in sorted(by_region.items(), key=lambda x: -x[1])]
                )
                st.dataframe(df_reg, use_container_width=True, hide_index=True)

    # Bar chart
    if by_modality:
        import pandas as pd
        st.markdown("**Dose Distribution by Modality**")
        chart_df = pd.DataFrame({"Modality": list(by_modality.keys()), "Dose (mSv)": list(by_modality.values())})
        chart_df = chart_df.set_index("Modality")
        st.bar_chart(chart_df)


def _display_drl_comparison(data: dict):
    """Display DRL comparison results."""
    status = data.get("status", "unknown")
    ratio = data.get("ratio", 0)
    patient_dose = data.get("patient_dose_msv", 0)
    drl = data.get("drl_msv", 0)
    achievable = data.get("achievable_dose_msv", 0)
    suggestions = data.get("optimization_suggestions", [])

    # Status banner
    color = _drl_status_color(status)
    label = _drl_status_label(status)
    st.markdown(
        f'<div style="padding:12px;border-radius:8px;border:2px solid {color};'
        f'background:{color}15;margin-bottom:16px;">'
        f'<span style="font-size:1.1em;font-weight:bold;color:{color};">'
        f"{label}</span></div>",
        unsafe_allow_html=True,
    )

    # Metrics
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Patient Dose", f"{patient_dose:.1f} mSv")
    with m2:
        st.metric("National DRL", f"{drl:.1f} mSv")
    with m3:
        st.metric("Achievable", f"{achievable:.1f} mSv")
    with m4:
        delta_color = "inverse" if ratio <= 1.0 else "normal"
        st.metric("Dose/DRL Ratio", f"{ratio:.2f}", delta=f"{(ratio - 1.0) * 100:+.0f}%", delta_color=delta_color)

    # Visual comparison bar
    if drl > 0:
        st.markdown("**Dose Comparison**")
        max_val = max(patient_dose, drl, achievable) * 1.2
        _render_dose_comparison_bar(patient_dose, achievable, drl, max_val)

    # Optimization suggestions
    if suggestions:
        st.markdown("---")
        st.markdown("**Optimization Suggestions:**")
        for s in suggestions:
            if s.startswith("PEDIATRIC:"):
                st.info(s)
            elif "reduce" in s.lower() or "consider" in s.lower():
                st.warning(s)
            else:
                st.markdown(f"- {s}")


def _render_dose_comparison_bar(patient: float, achievable: float, drl: float, max_val: float):
    """Render a simple visual dose comparison bar."""
    def pct(val):
        return min(val / max_val * 100, 100) if max_val > 0 else 0

    html = f"""
    <div style="margin:8px 0;font-size:0.85em;">
        <div style="margin-bottom:4px;">
            <span style="display:inline-block;width:100px;">Patient:</span>
            <div style="display:inline-block;width:{pct(patient):.0f}%;
                background:{'#c62828' if patient > drl else '#4caf50' if patient <= achievable else '#f9a825'};
                height:20px;border-radius:4px;min-width:2px;"></div>
            <span style="margin-left:4px;">{patient:.1f} mSv</span>
        </div>
        <div style="margin-bottom:4px;">
            <span style="display:inline-block;width:100px;">Achievable:</span>
            <div style="display:inline-block;width:{pct(achievable):.0f}%;
                background:#81c784;height:20px;border-radius:4px;min-width:2px;"></div>
            <span style="margin-left:4px;">{achievable:.1f} mSv</span>
        </div>
        <div>
            <span style="display:inline-block;width:100px;">DRL:</span>
            <div style="display:inline-block;width:{pct(drl):.0f}%;
                background:#ffb74d;height:20px;border-radius:4px;min-width:2px;"></div>
            <span style="margin-left:4px;">{drl:.1f} mSv</span>
        </div>
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)


def _render_reference_data():
    """Render the static dose reference comparison chart (standard vs AI-optimized)."""
    st.subheader("Standard vs AI-Optimized Dose")
    st.caption("Reference comparison of standard protocols versus AI-enhanced imaging")

    data = _load_dose_data()
    if not data:
        st.info("Dose reference data not available.")
        return

    # Summary metrics
    reductions = [d["reduction_pct"] for d in data]
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Protocols Analyzed", len(data))
    with col2:
        st.metric("Avg Dose Reduction", f"{sum(reductions) / len(reductions):.0f}%")
    with col3:
        st.metric("Max Dose Reduction", f"{max(reductions)}%")
    with col4:
        ct_data = [d for d in data if d["modality"] == "ct"]
        if ct_data:
            avg_ct = sum(d["reduction_pct"] for d in ct_data) / len(ct_data)
            st.metric("Avg CT Reduction", f"{avg_ct:.0f}%")

    st.divider()

    # Filter
    modality_filter = st.selectbox(
        "Filter by Modality",
        ["All"] + sorted(set(d["modality"] for d in data)),
        key="ref_mod_filter",
    )
    filtered = data if modality_filter == "All" else [d for d in data if d["modality"] == modality_filter]

    # Comparison chart
    chart_data = {}
    for d in filtered:
        name = d["protocol"][:30]
        chart_data[name] = {
            "Standard (mGy)": d["standard_ctdi_mgy"],
            "AI-Optimized (mGy)": d["ai_optimized_ctdi_mgy"],
        }

    if chart_data:
        import pandas as pd
        df = pd.DataFrame(chart_data).T
        st.bar_chart(df)

    # Detailed table
    st.subheader("Protocol Details")
    for d in filtered:
        pct = d["reduction_pct"]
        indicator = "[+++]" if pct >= 40 else "[++]" if pct >= 25 else "[+]"
        with st.expander(f"{indicator} **{d['protocol']}** -- {pct}% dose reduction"):
            c1, c2 = st.columns(2)
            with c1:
                st.metric("Standard Dose", f"{d['standard_ctdi_mgy']} mGy")
                st.metric("DLP", f"{d['dlp_mgy_cm']} mGy-cm")
            with c2:
                st.metric(
                    "AI-Optimized Dose",
                    f"{d['ai_optimized_ctdi_mgy']} mGy",
                    delta=f"-{pct}%",
                    delta_color="inverse",
                )
                st.metric("Effective Dose", f"{d['effective_dose_msv']} mSv")
            st.markdown(f"**Technique:** {d['technique']}")
            st.markdown(f"**Image Quality:** {d['image_quality']}")


def _render_pediatric_section():
    """Render pediatric dose considerations."""
    st.subheader("Pediatric Dose Considerations")

    st.markdown("""
**Key Principles for Pediatric Imaging:**

1. **ALARA Principle** -- As Low As Reasonably Achievable. Children are significantly
   more radiosensitive than adults, with 3-5x higher lifetime cancer risk per unit dose.

2. **Image Gently Campaign** -- Follow Image Gently guidelines for all pediatric CT:
   - Child-size the technique (kV and mAs)
   - Scan only the indicated area
   - Scan only when indicated
   - Use one scan phase (multiphase rarely needed)

3. **Size-Specific Dose Estimates (SSDE)** -- Use SSDE rather than CTDIvol for
   accurate pediatric dose assessment. SSDE accounts for patient size.

4. **Age-Based kV Reduction:**
""")

    import pandas as pd
    kv_data = pd.DataFrame({
        "Age Group": ["< 1 year", "1-5 years", "5-12 years", "12-18 years", "Adult"],
        "Recommended kV": ["70-80", "80", "80-100", "100-120", "120"],
        "mAs Approach": ["Weight-based", "Weight-based", "AEC with pedi protocol", "AEC with pedi protocol", "AEC standard"],
        "Dose Multiplier": ["0.3x adult", "0.5x adult", "0.6x adult", "0.8x adult", "1.0x (reference)"],
    })
    st.dataframe(kv_data, use_container_width=True, hide_index=True)

    st.markdown("""
5. **Alternative Modalities** -- Always consider non-ionizing alternatives:
   - **Ultrasound** -- First-line for appendicitis, pyloric stenosis, intussusception
   - **MRI** -- Brain, spine, musculoskeletal (no radiation)
   - **Rapid MRI** -- Emerging alternative for head CT in many indications

6. **Dose Tracking** -- Pediatric patients have lower cumulative dose thresholds
   (50% of adult limits). Use the Dose Intelligence tracker with the pediatric
   flag enabled for appropriate alerting.
""")

    st.info(
        "**ALARA Principle:** As Low As Reasonably Achievable. AI-enhanced "
        "reconstruction maintains diagnostic image quality while significantly "
        "reducing patient radiation exposure, particularly impactful for "
        "pediatric patients and repeat screening examinations."
    )
