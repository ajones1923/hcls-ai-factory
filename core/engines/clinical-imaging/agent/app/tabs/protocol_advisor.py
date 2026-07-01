"""Protocol Advisor tab -- AI-powered imaging protocol optimization.

Calls the FastAPI backend at /protocol/recommend with patient-specific
factors and displays ACR Appropriateness Criteria-based recommendations.
"""

import streamlit as st
import requests

API_BASE = "http://localhost:8524"


def _acr_color(rating: int) -> str:
    """Return a CSS color for the ACR appropriateness rating."""
    if rating >= 7:
        return "#2e7d32"  # green
    elif rating >= 4:
        return "#f9a825"  # yellow/amber
    else:
        return "#c62828"  # red


def _acr_label(rating: int) -> str:
    """Return appropriateness label for an ACR rating."""
    if rating >= 7:
        return "Usually Appropriate"
    elif rating >= 4:
        return "May Be Appropriate"
    else:
        return "Usually Not Appropriate"


def _render_acr_scale(rating: int):
    """Render a visual 1-9 ACR appropriateness scale."""
    blocks = []
    for i in range(1, 10):
        if i <= 3:
            bg = "#c62828" if i == rating else "#ffcdd2"
        elif i <= 6:
            bg = "#f9a825" if i == rating else "#fff9c4"
        else:
            bg = "#2e7d32" if i == rating else "#c8e6c9"
        border = "2px solid #333" if i == rating else "1px solid #ccc"
        font_weight = "bold" if i == rating else "normal"
        color = "white" if i == rating else "#333"
        blocks.append(
            f'<div style="display:inline-block;width:32px;height:32px;'
            f"line-height:32px;text-align:center;border-radius:4px;"
            f"border:{border};background:{bg};color:{color};"
            f'font-weight:{font_weight};margin:0 2px;">{i}</div>'
        )
    label = _acr_label(rating)
    color = _acr_color(rating)
    html = (
        '<div style="margin:8px 0;">'
        + "".join(blocks)
        + f'<div style="margin-top:4px;font-weight:bold;color:{color};">'
        f"Rating {rating}/9 -- {label}</div></div>"
    )
    st.markdown(html, unsafe_allow_html=True)


def _fetch_indications() -> list:
    """Fetch supported indications from the API."""
    try:
        resp = requests.get(f"{API_BASE}/protocol/indications", timeout=5)
        if resp.status_code == 200:
            return resp.json()
    except Exception:
        pass
    return []


def render():
    st.header("Protocol Advisor")
    st.caption(
        "Evidence-based imaging protocol recommendations using "
        "ACR Appropriateness Criteria with patient-specific adjustments"
    )

    # ── Example indications ──────────────────────────────────────────
    st.markdown("**Try an example indication:**")
    examples = [
        "Acute chest pain, rule out PE",
        "New onset seizure, rule out mass",
        "Lung cancer screening",
        "Abdominal pain, rule out appendicitis",
    ]
    example_cols = st.columns(len(examples))
    for i, (col, example) in enumerate(zip(example_cols, examples)):
        with col:
            if st.button(example, key=f"proto_ex_{i}", use_container_width=True):
                st.session_state["protocol_indication"] = example

    # Show supported indications
    indications = _fetch_indications()
    if indications:
        with st.expander("Supported clinical indications"):
            st.write(", ".join(f"`{ind}`" for ind in indications))

    st.divider()

    # ── Input form ───────────────────────────────────────────────────
    col_left, col_right = st.columns([2, 1])

    with col_left:
        indication = st.text_input(
            "Clinical Indication",
            value=st.session_state.get("protocol_indication", ""),
            placeholder="e.g. headache, chest pain, lung screening, stroke acute",
        )

    with col_right:
        st.markdown("**Patient Factors**")
        p_age = st.number_input("Age", min_value=0, max_value=120, value=50, key="proto_age")
        p_weight = st.number_input("Weight (kg)", min_value=1.0, max_value=300.0, value=70.0, key="proto_wt")
        p_sex = st.selectbox("Sex", ["", "M", "F"], key="proto_sex")
        p_egfr = st.number_input(
            "eGFR (mL/min/1.73m2)", min_value=0.0, max_value=200.0,
            value=90.0, key="proto_egfr",
            help="Glomerular filtration rate. <30 = severe impairment, <60 = moderate.",
        )

        p_pregnant = st.checkbox("Pregnant", key="proto_preg")
        p_allergy = st.checkbox("Contrast allergy", key="proto_allergy")
        p_allergy_severity = None
        if p_allergy:
            p_allergy_severity = st.selectbox(
                "Allergy severity", ["mild", "moderate", "severe"], key="proto_allsev"
            )
        p_pediatric = st.checkbox("Pediatric protocol", key="proto_pedi")

    # ── Submit ───────────────────────────────────────────────────────
    if st.button("Recommend Protocol", type="primary", use_container_width=True):
        if not indication:
            st.warning("Please enter a clinical indication.")
            return

        patient_payload = {
            "age": p_age,
            "weight_kg": p_weight,
            "sex": p_sex or None,
            "pregnant": p_pregnant,
            "renal_function_egfr": p_egfr,
            "contrast_allergy": p_allergy,
            "contrast_allergy_severity": p_allergy_severity,
            "pediatric": p_pediatric,
        }

        with st.spinner("Analyzing clinical parameters..."):
            try:
                resp = requests.post(
                    f"{API_BASE}/protocol/recommend",
                    json={"indication": indication, "patient": patient_payload},
                    timeout=10,
                )
                if resp.status_code != 200:
                    st.error(f"API error ({resp.status_code}): {resp.text}")
                    return

                rec = resp.json()
            except requests.ConnectionError:
                st.error(
                    "Cannot connect to the Clinical Imaging Engine API at "
                    f"{API_BASE}. Ensure the backend is running."
                )
                return
            except Exception as e:
                st.error(f"Protocol optimization error: {e}")
                return

        # ── Display results ──────────────────────────────────────────
        st.subheader("Recommended Protocol")

        # Key metrics
        m1, m2, m3 = st.columns(3)
        with m1:
            st.metric("Protocol", rec.get("recommended_protocol", "N/A"))
        with m2:
            st.metric("Modality", rec.get("recommended_modality", "N/A"))
        with m3:
            dose = rec.get("dose_estimate_msv")
            st.metric("Est. Dose", f"{dose} mSv" if dose is not None else "N/A")

        # ACR rating scale
        rating = rec.get("acr_appropriateness_rating", 5)
        _render_acr_scale(rating)

        # Warnings
        warnings = rec.get("warnings", [])
        if warnings:
            st.markdown("---")
            for w in warnings:
                if "PREGNANCY" in w.upper():
                    st.warning(w)
                elif "RENAL" in w.upper():
                    st.warning(w)
                elif "ALLERGY" in w.upper() or "CONTRAST" in w.upper():
                    st.warning(w)
                elif "PEDIATRIC" in w.upper():
                    st.info(w)
                else:
                    st.info(w)

        # Protocol parameters
        params = rec.get("parameters", {})
        contrast = rec.get("contrast")
        if params or contrast:
            with st.expander("Protocol Parameters", expanded=True):
                if params:
                    param_cols = st.columns(min(len(params), 4))
                    for idx, (k, v) in enumerate(params.items()):
                        with param_cols[idx % len(param_cols)]:
                            st.markdown(f"**{k.replace('_', ' ').title()}:** {v}")
                if contrast:
                    st.markdown("---")
                    st.markdown("**Contrast:**")
                    for k, v in contrast.items():
                        st.markdown(f"- **{k.replace('_', ' ').title()}:** {v}")

        # Rationale
        rationale = rec.get("rationale", "")
        if rationale:
            with st.expander("Clinical Rationale"):
                st.markdown(rationale)

        # Alternatives table
        alternatives = rec.get("alternatives", [])
        if alternatives:
            st.markdown("---")
            st.subheader("Alternative Protocols")
            import pandas as pd

            alt_rows = []
            for alt in alternatives:
                alt_rating = alt.get("rating", 0)
                alt_rows.append({
                    "Modality": alt.get("modality", ""),
                    "Protocol": alt.get("protocol", ""),
                    "ACR Rating": f"{alt_rating}/9",
                    "Appropriateness": _acr_label(alt_rating),
                    "Note": alt.get("note", ""),
                })
            df = pd.DataFrame(alt_rows)
            st.dataframe(df, use_container_width=True, hide_index=True)
