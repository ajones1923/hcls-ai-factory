"""Analytics tab — population-level imaging analytics."""
import streamlit as st


def render():
    st.header("Population Analytics")
    st.caption("GPU-accelerated imaging analytics powered by NVIDIA RAPIDS")

    import requests
    api_url = "http://localhost:8524"

    # Generate demo data button
    col1, col2 = st.columns([1, 3])
    with col1:
        if st.button("Generate Demo Data", type="primary"):
            try:
                r = requests.post(f"{api_url}/api/analytics/generate-demo-data?n_studies=500", timeout=10)
                st.success(f"Generated {r.json().get('generated', 0)} studies")
                st.rerun()
            except Exception as e:
                st.error(f"Error: {e}")

    # Population summary
    try:
        r = requests.get(f"{api_url}/api/analytics/population", timeout=10)
        data = r.json()

        if data.get("total_studies", 0) == 0:
            st.info("No analytics data yet. Click 'Generate Demo Data' to create synthetic imaging studies for demonstration.")
            return

        # Key metrics
        st.subheader("Population Overview")
        m1, m2, m3, m4 = st.columns(4)
        with m1:
            st.metric("Total Studies", f"{data['total_studies']:,}")
        with m2:
            st.metric("Critical Finding Rate", f"{data.get('critical_finding_rate',0):.1%}")
        with m3:
            st.metric("Studies w/ Critical", data.get("studies_with_critical_findings", 0))
        with m4:
            st.metric("Avg Processing", f"{data.get('mean_processing_time_ms',0):.0f}ms")

        st.divider()

        # Modality distribution
        st.subheader("Modality Distribution")
        mod_dist = data.get("modality_distribution", {})
        if mod_dist:
            import pandas as pd
            df = pd.DataFrame(list(mod_dist.items()), columns=["Modality", "Count"])
            st.bar_chart(df.set_index("Modality"))

        # Severity distribution
        st.subheader("Severity Distribution")
        sev_dist = data.get("severity_distribution", {})
        if sev_dist:
            import pandas as pd
            df = pd.DataFrame(list(sev_dist.items()), columns=["Severity", "Count"])
            st.bar_chart(df.set_index("Severity"))

        st.divider()

        # Cohort query
        st.subheader("Cohort Query")
        cq1, cq2, cq3 = st.columns(3)
        with cq1:
            cohort_mod = st.selectbox("Modality", [""] + list(mod_dist.keys()), key="cohort_mod")
        with cq2:
            cohort_sev = st.selectbox("Severity", [""] + list(sev_dist.keys()), key="cohort_sev")
        with cq3:
            if st.button("Search Cohort"):
                criteria = {}
                if cohort_mod:
                    criteria["modality"] = cohort_mod
                if cohort_sev:
                    criteria["severity"] = cohort_sev
                cr = requests.post(f"{api_url}/api/analytics/cohort", json=criteria, timeout=10).json()
                st.metric("Matching Studies", f"{cr.get('matching_studies',0)} / {cr.get('total_studies',0)}")
                st.metric("Match Rate", f"{cr.get('match_rate',0):.1%}")

        st.divider()

        # Trends
        st.subheader("Temporal Trends")
        try:
            tr = requests.get(f"{api_url}/api/analytics/trends/study_count?granularity=monthly", timeout=10).json()
            points = tr.get("data_points", [])
            if points:
                import pandas as pd
                df = pd.DataFrame(points)
                st.line_chart(df.set_index("period")["count"])
                st.caption(f"Trend direction: {tr.get('trend_direction', 'N/A')}")
        except Exception:
            st.info("Trend data unavailable")

        # Severity by modality
        st.subheader("Severity by Modality")
        try:
            sm = requests.get(f"{api_url}/api/analytics/severity-by-modality", timeout=10).json()
            if sm:
                import pandas as pd
                df = pd.DataFrame(sm).fillna(0).T
                st.dataframe(df, use_container_width=True)
        except Exception:
            st.info("Cross-tabulation unavailable")

    except Exception as e:
        st.warning(f"Analytics API unavailable: {e}")
        st.info("Start the FastAPI server to enable analytics.")
