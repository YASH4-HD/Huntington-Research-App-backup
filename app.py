 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/app.py b/app.py
index f42547d53f8dbf70d3f7ad76a1275f515f13c521..979a0d8acacc7c1208227bd4b1bba7998bfa332c 100644
--- a/app.py
+++ b/app.py
@@ -1,33 +1,34 @@
 import streamlit as st
 import pandas as pd
 import requests
 import networkx as nx
 import matplotlib.pyplot as plt
 import numpy as np
 from scipy.stats import fisher_exact
 import io
+from itertools import combinations
 
 # --- PAGE CONFIGURATION ---
 st.set_page_config(page_title="NeuroMetabolic Framework", page_icon="🧬", layout="wide")
 
 # --- GLOBAL SETTINGS ---
 role_colors = {
     "⭐ Core Gene": "#FF4B4B", 
     "🔋 Mitochondrial Dysfunction": "#FFA500", 
     "💀 Apoptosis": "#7D3C98", 
     "🧠 Synaptic / Excitotoxicity": "#2E86C1", 
     "♻️ Autophagy": "#28B463", 
     "📦 Proteostasis / PSMC": "#D4AC0D", 
     "🧬 Pathway Component": "#D5D8DC"
 }
 
 # --- DATA ACQUISITION (KEGG API) ---
 @st.cache_data
 def get_kegg_genes(pathway_id):
     url = f"https://rest.kegg.jp/get/{pathway_id}"
     response = requests.get(url)
     genes = []
     if response.status_code == 200:
         lines = response.text.split('\n')
         is_gene_section = False
         for line in lines:
@@ -174,56 +175,87 @@ st.sidebar.markdown("""
         This tool is intended for research hypothesis generation. Network edges represent inferred functional associations derived from KEGG pathway co-occurrence and do not imply direct physical, regulatory, or causal protein-protein interactions.<br><br>
         <i>"Findings guide experimental prioritization rather than replace wet-lab validation."</i>
     </p>
 </div>
 """, unsafe_allow_html=True)
 
 # --- LOAD DATA ---
 df = get_kegg_genes(pathway_id)
 if not df.empty:
     df["Functional Role"] = df.apply(lambda row: assign_role(row["Symbol"], row["Description"], disease_choice), axis=1)
     
     def calculate_validation(symbol):
         high_lit = ["HTT", "BDNF", "APP", "MAPT", "SNCA", "PRKN", "SOD1", "INS", "BDNF", "CASP3", "TP53"]
         if symbol in high_lit: return 95
         np.random.seed(sum(ord(c) for c in symbol))
         return np.random.randint(20, 60)
 
     def calculate_priority(row):
         base = 100 if "Core" in row['Functional Role'] else 50
         lit = calculate_validation(row['Symbol'])
         return (base * 0.6) + (lit * 0.4)
 
     df['Lit_Score'] = df['Symbol'].apply(calculate_validation)
     df['Score'] = df.apply(calculate_priority, axis=1)
 
+st.sidebar.caption("Data: KEGG API | System: Streamlit")
+
+# --- ACCESS CONTROL: HIDE INNOVATION LAB UNLESS PIN IS UNLOCKED ---
+if "innovation_lab_unlocked" not in st.session_state:
+    st.session_state.innovation_lab_unlocked = False
+
+configured_pin = str(st.secrets.get("innovation_lab_pin", ""))
+
+with st.sidebar.expander("🔒 Innovation Lab Access", expanded=False):
+    entered_pin = st.text_input("Enter PIN", type="password", key="innovation_lab_pin_input")
+    unlock_clicked = st.button("Unlock Innovation Lab", key="innovation_lab_unlock_btn", use_container_width=True)
+
+    if unlock_clicked:
+        if configured_pin and entered_pin == configured_pin:
+            st.session_state.innovation_lab_unlocked = True
+            st.success("Innovation Lab unlocked for this session.")
+        elif not configured_pin:
+            st.warning("No PIN configured. Add `innovation_lab_pin` to Streamlit secrets.")
+        else:
+            st.error("Incorrect PIN")
+
+    if st.session_state.innovation_lab_unlocked:
+        if st.button("Lock Innovation Lab", key="innovation_lab_lock_btn", use_container_width=True):
+            st.session_state.innovation_lab_unlocked = False
+            st.info("Innovation Lab locked.")
+
+
 # --- MAIN CONTENT ---
 st.title(f"🧬 {disease_choice} Metabolic Framework")
 
 st.markdown(f"*This resource guide serves as a foundational reference for computational hypothesis generation, validation, and extension of the {disease_choice} metabolic mechanisms.*")
 
-tab1, tab2, tab3 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Manuscript"])
+if st.session_state.get("innovation_lab_unlocked", False):
+    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Manuscript", "🚀 Innovation Lab", "🧪 Digital Twin", "🔁 Counterfactual Engine"])
+else:
+    tab1, tab2, tab3, tab5, tab6 = st.tabs(["📊 Target Discovery", "🕸️ Interaction Network", "🔬 Enrichment & Manuscript", "🧪 Digital Twin", "🔁 Counterfactual Engine"])
+    tab4 = None
 
 with tab1:
     col_a, col_b = st.columns([2, 1])
     with col_a:
         st.subheader("Genetic Components")
         search_query = st.text_input("🔍 Search genes or mechanisms:", placeholder="Type to filter...")
     with col_b:
         st.subheader("Deep Dive")
         selected_gene = st.selectbox("External Research:", ["Select a Gene"] + list(df['Symbol'].unique()))
         if selected_gene != "Select a Gene":
             st.markdown(f"**[View {selected_gene} on GeneCards ↗️](https://www.genecards.org/cgi-bin/carddisp.pl?gene={selected_gene})**")
 
     mask = df['Symbol'].str.contains(search_query.upper(), na=False) | df['Functional Role'].str.contains(search_query, case=False, na=False)
     filtered_df = df[mask] if search_query else df
     st.dataframe(filtered_df[['Symbol', 'Functional Role', 'Lit_Score', 'Score', 'Description']].sort_values('Score', ascending=False), use_container_width=True, height=300)
 
     with st.expander("ℹ️ Understanding the Scoring System", expanded=False):
         col1, col2 = st.columns(2)
         with col1:
             st.markdown("""
             **📊 Lit_Score (Literature Prevalence)**
             - **What it is:** A numerical value representing the research density for a specific gene.
             - **How it's calculated:** High-confidence core genes are assigned a baseline of **95**.
             """)
         with col2:
@@ -375,26 +407,274 @@ indicates a complex, multi-factorial regulatory landscape.
 
 4. PROSPECTIVE HYPOTHESIS
 The data supports a model where {disease_choice} progression is significantly 
 mediated by {top_mech} failure, leading to a cascade of secondary 
 mitochondrial and synaptic deficits.
 
 This resource guide serves as a foundational reference for computational hypothesis generation, validation, and extension of the {disease_choice} metabolic mechanisms.
 --------------------------------------------------
 END OF REPORT
 """
         st.markdown("### Preview")
         st.info(manuscript_text)
 
         st.download_button(
             label="📥 Download Summary as .txt",
             data=manuscript_text,
             file_name=f"{disease_choice.replace(' ', '_')}_Summary.txt",
             mime="text/plain",
             use_container_width=True
         )
 
     st.markdown("---")
     st.subheader("📚 Research Bibliography")
     st.markdown(f"1. Disease Pathway: {pathway_id} | 2. KEGG API | 3. Fisher's Exact Test Analysis")
 
-st.sidebar.caption("Data: KEGG API | System: Streamlit")
+def build_digital_twin_dataframe(dataframe, ui_key_prefix="dt"):
+    stage_multiplier = {
+        "Pre-symptomatic": 0.9,
+        "Early disease": 1.0,
+        "Advanced disease": 1.2,
+    }
+
+    role_map = {
+        "🔋 Mitochondrial Dysfunction": "mitochondrial",
+        "📦 Proteostasis / PSMC": "proteostasis",
+        "♻️ Autophagy": "autophagy",
+        "💀 Apoptosis": "apoptosis",
+        "🧠 Synaptic / Excitotoxicity": "synaptic",
+    }
+
+    c_stage, c_mito, c_prot = st.columns(3)
+    with c_stage:
+        disease_stage = st.selectbox(
+            "Progression stage",
+            ["Pre-symptomatic", "Early disease", "Advanced disease"],
+            index=1,
+            key=f"{ui_key_prefix}_stage"
+        )
+    with c_mito:
+        mito_stress = st.slider("Mitochondrial stress", 0, 100, 55, key=f"{ui_key_prefix}_mito")
+    with c_prot:
+        proteostasis_load = st.slider("Proteostasis load", 0, 100, 60, key=f"{ui_key_prefix}_proteostasis")
+
+    c_auto, c_apop, c_syn = st.columns(3)
+    with c_auto:
+        autophagy_support = st.slider("Autophagy support", 0, 100, 50, key=f"{ui_key_prefix}_autophagy")
+    with c_apop:
+        apoptosis_pressure = st.slider("Apoptosis pressure", 0, 100, 58, key=f"{ui_key_prefix}_apoptosis")
+    with c_syn:
+        synaptic_toxicity = st.slider("Synaptic toxicity", 0, 100, 52, key=f"{ui_key_prefix}_synaptic")
+
+    pressure_values = {
+        "mitochondrial": mito_stress,
+        "proteostasis": proteostasis_load,
+        "autophagy": autophagy_support,
+        "apoptosis": apoptosis_pressure,
+        "synaptic": synaptic_toxicity,
+    }
+
+    twin_df = dataframe.copy()
+    stage_factor = stage_multiplier[disease_stage]
+
+    def role_pressure_multiplier(role):
+        pressure_key = role_map.get(role)
+        if pressure_key is None:
+            return 1.0
+
+        if pressure_key == "autophagy":
+            return max(0.65, 1.15 - (pressure_values[pressure_key] / 200))
+
+        return 0.85 + (pressure_values[pressure_key] / 100) * 0.7
+
+    twin_df["Twin_Multiplier"] = twin_df["Functional Role"].apply(role_pressure_multiplier) * stage_factor
+    twin_df["Twin_Score"] = (twin_df["Score"] * twin_df["Twin_Multiplier"]).round(2)
+    twin_df["Score_Delta"] = (twin_df["Twin_Score"] - twin_df["Score"]).round(2)
+    return twin_df, stage_factor
+
+
+if tab4 is not None:
+    with tab4:
+        st.subheader("🚀 Innovation Lab — Award-Winning Strategy Builder")
+        st.markdown("Generate a practical strategy to make this framework stand out in grants, competitions, and PhD evaluations.")
+
+        col_goal, col_horizon = st.columns(2)
+        with col_goal:
+            innovation_goal = st.selectbox(
+                "Primary objective",
+                [
+                    "Scientific Novelty",
+                    "Clinical Translation",
+                    "Open-Science Reproducibility",
+                    "AI + Systems Biology Integration"
+                ]
+            )
+        with col_horizon:
+            time_horizon = st.selectbox("Execution window", ["30 days", "90 days", "6 months"])
+
+        score_weights = {
+            "Scientific Novelty": {"impact": 0.45, "feasibility": 0.25, "rigor": 0.30},
+            "Clinical Translation": {"impact": 0.40, "feasibility": 0.35, "rigor": 0.25},
+            "Open-Science Reproducibility": {"impact": 0.25, "feasibility": 0.35, "rigor": 0.40},
+            "AI + Systems Biology Integration": {"impact": 0.40, "feasibility": 0.20, "rigor": 0.40},
+        }
+
+        candidate_innovations = [
+            {"Idea": "Digital Twin Simulator", "Impact": 95, "Feasibility": 50, "Rigor": 72, "Deliverable": "Interactive progression simulator with perturbation sliders"},
+            {"Idea": "Counterfactual Intervention Engine", "Impact": 90, "Feasibility": 60, "Rigor": 78, "Deliverable": "What-if ranking of single and combination interventions"},
+            {"Idea": "Cell-Type Aware Enrichment", "Impact": 85, "Feasibility": 70, "Rigor": 88, "Deliverable": "Neuron/Astrocyte/Microglia-specific pathway enrichment overlays"},
+            {"Idea": "Multi-Omics Concordance Index", "Impact": 88, "Feasibility": 55, "Rigor": 92, "Deliverable": "Cross-validation panel across transcriptomics/proteomics/metabolomics"},
+            {"Idea": "Evidence Strength Dashboard", "Impact": 78, "Feasibility": 92, "Rigor": 94, "Deliverable": "Confidence tiers, provenance trail, and manuscript supplement exports"},
+        ]
+
+        weights = score_weights[innovation_goal]
+        innovation_df = pd.DataFrame(candidate_innovations)
+        innovation_df["Strategic Score"] = (
+            innovation_df["Impact"] * weights["impact"]
+            + innovation_df["Feasibility"] * weights["feasibility"]
+            + innovation_df["Rigor"] * weights["rigor"]
+        ).round(1)
+        innovation_df = innovation_df.sort_values("Strategic Score", ascending=False).reset_index(drop=True)
+
+        st.markdown("### Ranked innovation roadmap")
+        st.dataframe(innovation_df[["Idea", "Strategic Score", "Deliverable"]], use_container_width=True, hide_index=True)
+
+        top_pick = innovation_df.iloc[0]
+        runner_up = innovation_df.iloc[1]
+        st.success(f"For **{disease_choice}** and objective **{innovation_goal}**, start with **{top_pick['Idea']}** (Strategic Score: {top_pick['Strategic Score']}).")
+        st.info(f"Second priority: **{runner_up['Idea']}**. Suggested horizon: **{time_horizon}** with milestone-driven delivery.")
+
+        action_plan = f"""AWARD-READY ACTION PLAN
+
+    Disease context: {disease_choice}
+    Objective: {innovation_goal}
+    Execution window: {time_horizon}
+
+    Priority 1: {top_pick['Idea']}
+    Deliverable: {top_pick['Deliverable']}
+    Strategic Score: {top_pick['Strategic Score']}
+
+    Priority 2: {runner_up['Idea']}
+    Deliverable: {runner_up['Deliverable']}
+    Strategic Score: {runner_up['Strategic Score']}
+
+    Judging narrative:
+    The framework is evolving from static pathway analysis into an explainable decision-support system that generates testable, reproducible, and translational hypotheses.
+    """
+        st.download_button("📥 Download Action Plan", action_plan, f"{disease_choice.replace(' ', '_')}_award_strategy.txt", "text/plain", use_container_width=True)
+
+with tab5:
+    st.subheader("🧪 Digital Twin Sandbox (Prototype)")
+    st.caption("Simulate disease-state shifts by perturbing key biological pressures and progression stage.")
+
+    twin_df, stage_factor = build_digital_twin_dataframe(df, ui_key_prefix="dt")
+    top_baseline = df.sort_values("Score", ascending=False).head(5)[["Symbol", "Score"]]
+    top_twin = twin_df.sort_values("Twin_Score", ascending=False).head(5)[["Symbol", "Twin_Score", "Score_Delta", "Functional Role"]]
+
+    m1, m2, m3 = st.columns(3)
+    with m1:
+        st.metric("Stage factor", f"{stage_factor:.2f}x")
+    with m2:
+        st.metric("Highest twin score", f"{top_twin.iloc[0]['Twin_Score']:.1f}")
+    with m3:
+        st.metric("Genes with upward shift", int((twin_df["Score_Delta"] > 0).sum()))
+
+    c_left, c_right = st.columns(2)
+    with c_left:
+        st.markdown("**Top 5 Baseline Priorities**")
+        st.dataframe(top_baseline, use_container_width=True, hide_index=True)
+    with c_right:
+        st.markdown("**Top 5 Digital Twin Priorities**")
+        st.dataframe(top_twin, use_container_width=True, hide_index=True)
+
+    role_shift = twin_df.groupby("Functional Role", as_index=False)[["Score", "Twin_Score"]].mean().rename(columns={"Score": "Baseline", "Twin_Score": "Digital Twin"})
+    st.markdown("**Mechanism burden shift (mean score by role)**")
+    st.bar_chart(role_shift.set_index("Functional Role")[["Baseline", "Digital Twin"]])
+
+with tab6:
+    st.subheader('🔁 Mechanism Counterfactual Engine ("What if we target X?")')
+    st.caption("Estimate which single or combination interventions most reduce a selected pathological burden and quantify compensatory pathway effects.")
+
+    st.markdown("**Scenario context for counterfactual analysis**")
+    twin_df_cf, _ = build_digital_twin_dataframe(df, ui_key_prefix="cf")
+
+    mechanism_options = {
+        "Apoptosis burden": "💀 Apoptosis",
+        "Mitochondrial dysfunction": "🔋 Mitochondrial Dysfunction",
+        "Proteostasis burden": "📦 Proteostasis / PSMC",
+        "Synaptic toxicity": "🧠 Synaptic / Excitotoxicity",
+    }
+    intervention_library = {
+        "Boost autophagy": {"♻️ Autophagy": -0.22, "💀 Apoptosis": -0.10, "📦 Proteostasis / PSMC": -0.06},
+        "Enhance mitochondrial biogenesis": {"🔋 Mitochondrial Dysfunction": -0.20, "💀 Apoptosis": -0.08, "🧠 Synaptic / Excitotoxicity": -0.05},
+        "Proteasome rescue": {"📦 Proteostasis / PSMC": -0.21, "💀 Apoptosis": -0.07, "🔋 Mitochondrial Dysfunction": -0.03},
+        "Anti-excitotoxic modulation": {"🧠 Synaptic / Excitotoxicity": -0.18, "💀 Apoptosis": -0.06},
+        "Caspase pathway inhibition": {"💀 Apoptosis": -0.24, "🧠 Synaptic / Excitotoxicity": -0.04},
+    }
+
+    controls_left, controls_right = st.columns(2)
+    with controls_left:
+        outcome_label = st.selectbox("Outcome to minimize", list(mechanism_options.keys()), key="cf_outcome")
+        max_combo = st.selectbox("Maximum intervention size", [1, 2, 3], index=2, key="cf_combo")
+    with controls_right:
+        ci_iterations = st.slider("Confidence interval iterations", 100, 2000, 600, step=100, key="cf_iter")
+        show_top_n = st.slider("Show top combinations", 3, 15, 8, key="cf_top_n")
+
+    selected_outcome = mechanism_options[outcome_label]
+    current_burden = twin_df_cf.groupby("Functional Role")["Twin_Score"].mean().to_dict()
+    base_outcome_value = current_burden.get(selected_outcome, float(twin_df_cf["Twin_Score"].mean()))
+
+    combo_candidates = []
+    for r in range(1, max_combo + 1):
+        combo_candidates.extend(combinations(list(intervention_library.keys()), r))
+
+    rng = np.random.default_rng(seed=sum(ord(c) for c in disease_choice + selected_outcome))
+    records = []
+    for combo in combo_candidates:
+        aggregate_effects = {}
+        for intervention in combo:
+            for mechanism, effect in intervention_library[intervention].items():
+                aggregate_effects[mechanism] = aggregate_effects.get(mechanism, 0) + effect
+
+        synergy_bonus = 1 + (len(combo) - 1) * 0.08
+        aggregate_effects = {k: v * synergy_bonus for k, v in aggregate_effects.items()}
+
+        predicted = {}
+        for mechanism, baseline_value in current_burden.items():
+            predicted[mechanism] = max(baseline_value * (1 + aggregate_effects.get(mechanism, 0)), 0.01)
+
+        outcome_after = predicted.get(selected_outcome, base_outcome_value)
+        delta = outcome_after - base_outcome_value
+
+        compensation = 0
+        for role in ["🔋 Mitochondrial Dysfunction", "📦 Proteostasis / PSMC", "🧠 Synaptic / Excitotoxicity", "💀 Apoptosis"]:
+            if role != selected_outcome:
+                compensation += max(0, predicted.get(role, 0) - current_burden.get(role, 0))
+
+        noise_scale = max(abs(delta) * 0.22, 0.5)
+        samples = base_outcome_value + delta + rng.normal(0, noise_scale, size=ci_iterations)
+        ci_low, ci_high = np.percentile(samples, [2.5, 97.5])
+
+        records.append({
+            "Intervention Combo": " + ".join(combo),
+            "Combo Size": len(combo),
+            "Predicted Outcome": round(outcome_after, 2),
+            "Change vs Baseline": round(delta, 2),
+            "Compensation Risk": round(compensation, 2),
+            "CI 95% (Low)": round(float(ci_low), 2),
+            "CI 95% (High)": round(float(ci_high), 2),
+        })
+
+    counterfactual_df = pd.DataFrame(records).sort_values(by=["Predicted Outcome", "Compensation Risk", "Combo Size"], ascending=[True, True, True])
+    st.markdown(f"**Baseline burden for {outcome_label}: {base_outcome_value:.2f}**")
+
+    best_row = counterfactual_df.iloc[0]
+    c_best, c_delta, c_risk = st.columns(3)
+    with c_best:
+        st.metric("Best intervention combo", best_row["Intervention Combo"])
+    with c_delta:
+        st.metric("Predicted change", f"{best_row['Change vs Baseline']:.2f}")
+    with c_risk:
+        st.metric("Compensation risk", f"{best_row['Compensation Risk']:.2f}")
+
+    st.dataframe(counterfactual_df.head(show_top_n), use_container_width=True, hide_index=True)
+    st.bar_chart(counterfactual_df.head(show_top_n).set_index("Intervention Combo")[["Predicted Outcome", "CI 95% (Low)", "CI 95% (High)"]])
 
EOF
)
