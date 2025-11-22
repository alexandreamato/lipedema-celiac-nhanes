#!/usr/bin/env python3
"""
NHANES Secondary Analysis: Lipedema Phenotype and HLA DQ2/DQ8-Related Diseases
================================================================================

OBJECTIVE:
Since NHANES lacks HLA genotyping, we construct a proxy for 'HLA-related
autoimmunity background' based on:
- Serologically defined celiac disease (near-certain HLA-DQ2/DQ8)
- Type 1 diabetes probable
- Autoimmune thyroiditis probable
- Ancestry (Non-Hispanic White and Hispanic have higher DQ2/DQ8 prevalence)

HYPOTHESIS:
Women with lipedema phenotype (high leg-to-trunk fat ratio) may have lower
prevalence of HLA DQ2/DQ8-related autoimmune conditions.

HLA PROXY GROUPS:
- Group 2: Near-certain DQ2/DQ8 (serologically confirmed celiac)
- Group 1: High probability DQ2/DQ8 (T1D probable OR autoimmune thyroiditis + high-risk ethnicity)
- Group 0: Lower probability (reference group)
"""

import pandas as pd
import numpy as np
import requests
import os
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configuration
CACHE_DIR = "nhanes_cache"
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs("figures", exist_ok=True)
os.makedirs("tables", exist_ok=True)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 300

# ============================================================================
# DATA SOURCES
# ============================================================================
CYCLES = {
    "2011-2012": {
        "DEMO": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DEMO_G.xpt",
        "BMX":  "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/BMX_G.xpt",
        "DEXA": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DXX_G.xpt",
        "CELIAC": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/TGEMA_G.xpt",
        "DIQ": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DIQ_G.xpt",
        "MCQ": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/MCQ_G.xpt",
        "THYROD": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/THYROD_G.xpt",
        "GLU": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/GLU_G.xpt",
    },
    "2013-2014": {
        "DEMO": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DEMO_H.xpt",
        "BMX":  "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/BMX_H.xpt",
        "DEXA": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DXX_H.xpt",
        "CELIAC": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/TGEMA_H.xpt",
        "DIQ": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DIQ_H.xpt",
        "MCQ": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/MCQ_H.xpt",
        "GLU": "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/GLU_H.xpt",
        # THYROD not available for 2013-2014
    }
}

def download_xpt(url, local_name):
    """Download NHANES XPT file with caching"""
    path = os.path.join(CACHE_DIR, local_name)
    if os.path.exists(path):
        try:
            return pd.read_sas(path, format="xport")
        except:
            os.remove(path)
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with open(path, "wb") as f:
            f.write(r.content)
        return pd.read_sas(path, format="xport")
    except Exception as e:
        print(f"  Warning: Could not download {local_name}: {e}")
        return None

def process_cycle(cycle, urls):
    """Process a single NHANES cycle"""
    print(f"\nProcessing {cycle}...")

    # Required files
    demo = download_xpt(urls["DEMO"], f"{cycle}_DEMO.xpt")
    bmx = download_xpt(urls["BMX"], f"{cycle}_BMX.xpt")
    dexa = download_xpt(urls["DEXA"], f"{cycle}_DEXA.xpt")
    cel = download_xpt(urls["CELIAC"], f"{cycle}_CELIAC.xpt")
    diq = download_xpt(urls["DIQ"], f"{cycle}_DIQ.xpt")
    mcq = download_xpt(urls["MCQ"], f"{cycle}_MCQ.xpt")

    # Optional files
    thyrod = download_xpt(urls.get("THYROD", ""), f"{cycle}_THYROD.xpt") if "THYROD" in urls else None
    glu = download_xpt(urls.get("GLU", ""), f"{cycle}_GLU.xpt") if "GLU" in urls else None

    if any(x is None for x in [demo, bmx, dexa, cel, diq, mcq]):
        print(f"  Missing required files for {cycle}")
        return None

    # Start with demographics
    df = demo[["SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "WTMEC2YR"]].copy()

    # Merge BMX
    df = df.merge(bmx[["SEQN", "BMXBMI", "BMXWAIST", "BMXHT"]], on="SEQN", how="left")

    # Merge DEXA - check for gynoid/android columns
    dexa_cols = ["SEQN", "DXXLLFAT", "DXXRLFAT", "DXXTRFAT"]
    if "DXXGYFAT" in dexa.columns:
        dexa_cols.append("DXXGYFAT")  # Gynoid fat mass
    if "DXXANFAT" in dexa.columns:
        dexa_cols.append("DXXANFAT")  # Android fat mass
    if "DXXGYPF" in dexa.columns:
        dexa_cols.append("DXXGYPF")   # Gynoid percent fat
    if "DXXANPF" in dexa.columns:
        dexa_cols.append("DXXANPF")   # Android percent fat
    df = df.merge(dexa[dexa_cols], on="SEQN", how="left")

    # Merge Celiac serology
    df = df.merge(cel[["SEQN", "LBXTTG", "LBXEMA"]], on="SEQN", how="left")

    # Merge Diabetes questionnaire
    diq_cols = ["SEQN", "DIQ010", "DIQ050", "DIQ070"]
    if "DID040" in diq.columns:
        diq_cols.append("DID040")  # Age at diabetes diagnosis
    df = df.merge(diq[diq_cols], on="SEQN", how="left")

    # Merge Medical Conditions (thyroid)
    mcq_cols = ["SEQN"]
    if "MCQ160M" in mcq.columns:
        mcq_cols.append("MCQ160M")  # Ever told thyroid problem
    if "MCQ170M" in mcq.columns:
        mcq_cols.append("MCQ170M")  # Still have thyroid problem
    df = df.merge(mcq[mcq_cols], on="SEQN", how="left")

    # Merge Thyroid antibodies if available
    if thyrod is not None:
        thyrod_cols = ["SEQN"]
        if "LBXTPO" in thyrod.columns:
            thyrod_cols.append("LBXTPO")  # TPO antibodies
        if "LBXATG" in thyrod.columns:
            thyrod_cols.append("LBXATG")  # Thyroglobulin antibodies
        if len(thyrod_cols) > 1:
            df = df.merge(thyrod[thyrod_cols], on="SEQN", how="left")

    df["cycle"] = cycle
    return df

# ============================================================================
# LOAD AND COMBINE DATA
# ============================================================================
print("="*70)
print("NHANES HLA PROXY ANALYSIS: Lipedema Phenotype and Autoimmunity")
print("="*70)

dfs = [process_cycle(cyc, urls) for cyc, urls in CYCLES.items()]
dfs = [d for d in dfs if d is not None]
df = pd.concat(dfs, ignore_index=True)

# Filter to adult women
df = df[(df["RIAGENDR"] == 2) & (df["RIDAGEYR"] >= 20)].copy()
print(f"\nTotal adult women: {len(df)}")

# ============================================================================
# CREATE DERIVED VARIABLES
# ============================================================================

# Leg-to-trunk fat ratio
df["leg_fat"] = df["DXXLLFAT"] + df["DXXRLFAT"]
df["trunk_fat"] = df["DXXTRFAT"]
df["leg_trunk_ratio"] = df["leg_fat"] / df["trunk_fat"]

# Gynoid percent fat (if available)
if "DXXGYPF" in df.columns:
    df["gynoid_pct_fat"] = df["DXXGYPF"]

# Android percent fat (if available)
if "DXXANPF" in df.columns:
    df["android_pct_fat"] = df["DXXANPF"]

# Android/Gynoid ratio
if "DXXGYPF" in df.columns and "DXXANPF" in df.columns:
    df["android_gynoid_ratio"] = df["DXXANPF"] / df["DXXGYPF"]

# Waist-to-height ratio
df["waist_height_ratio"] = df["BMXWAIST"] / df["BMXHT"]

# Ethnicity categories
df["ethnicity"] = df["RIDRETH1"].map({
    1: "Mexican American",
    2: "Other Hispanic",
    3: "Non-Hispanic White",
    4: "Non-Hispanic Black",
    5: "Other Race"
})
df["high_risk_ethnicity"] = df["RIDRETH1"].isin([1, 2, 3])  # Mexican, Hispanic, NH White

# ============================================================================
# CREATE HLA PROXY GROUPS
# ============================================================================
print("\n" + "="*70)
print("CREATING HLA DQ2/DQ8 PROXY GROUPS")
print("="*70)

# Initialize
df["hla_proxy"] = 0

# GROUP 2: Near-certain DQ2/DQ8 (Serologically confirmed celiac)
# TTG-IgA positive (1) AND/OR EMA-IgA positive (1)
df["celiac_serologic"] = ((df["LBXTTG"] == 1) | (df["LBXEMA"] == 1)).astype(int)
df.loc[df["celiac_serologic"] == 1, "hla_proxy"] = 2
n_celiac = df["celiac_serologic"].sum()
print(f"\nGroup 2 (Celiac serology+): {n_celiac}")

# TYPE 1 DIABETES PROXY
# Criteria: Diabetes + Insulin use + (young onset OR low BMI)
df["diabetes"] = (df["DIQ010"] == 1).astype(int)
df["insulin_use"] = (df["DIQ050"] == 1).astype(int)
df["diabetes_pills"] = (df["DIQ070"] == 1).astype(int)

# T1D probable: diabetes + insulin + (age at dx < 40 if available, or BMI < 30)
# More conservative: insulin without pills in non-obese person
df["t1d_probable"] = 0
# If we have age at diagnosis
if "DID040" in df.columns:
    # Young onset diabetes on insulin
    df.loc[(df["diabetes"] == 1) &
           (df["insulin_use"] == 1) &
           (df["DID040"] < 40), "t1d_probable"] = 1

# Also consider: on insulin, not on pills, BMI < 30
df.loc[(df["diabetes"] == 1) &
       (df["insulin_use"] == 1) &
       (df["diabetes_pills"] != 1) &
       (df["BMXBMI"] < 30), "t1d_probable"] = 1

n_t1d = df["t1d_probable"].sum()
print(f"Type 1 Diabetes probable: {n_t1d}")

# AUTOIMMUNE THYROIDITIS PROXY
# MCQ160M = 1 (thyroid problem) + elevated TPO antibodies
df["thyroid_problem"] = (df.get("MCQ160M", pd.Series([0]*len(df))) == 1).astype(int)
df["current_thyroid"] = (df.get("MCQ170M", pd.Series([0]*len(df))) == 1).astype(int)

# TPO antibodies: elevated if > 34 IU/mL (standard cutoff)
if "LBXTPO" in df.columns:
    df["tpo_elevated"] = (df["LBXTPO"] > 34).astype(int)
    df["autoimmune_thyroid"] = ((df["thyroid_problem"] == 1) | (df["current_thyroid"] == 1)) & (df["tpo_elevated"] == 1)
else:
    # Without TPO, use thyroid problem report alone (weaker proxy)
    df["autoimmune_thyroid"] = 0

df["autoimmune_thyroid"] = df["autoimmune_thyroid"].astype(int)
n_ait = df["autoimmune_thyroid"].sum()
print(f"Autoimmune thyroiditis probable: {n_ait}")

# GROUP 1: High probability DQ2/DQ8
# T1D probable OR (autoimmune thyroid + high-risk ethnicity)
df.loc[(df["hla_proxy"] == 0) &
       ((df["t1d_probable"] == 1) |
        ((df["autoimmune_thyroid"] == 1) & (df["high_risk_ethnicity"] == True))),
       "hla_proxy"] = 1

# Print summary
print(f"\nHLA Proxy Group Summary:")
print(f"  Group 0 (Low probability): {(df['hla_proxy'] == 0).sum()}")
print(f"  Group 1 (High probability): {(df['hla_proxy'] == 1).sum()}")
print(f"  Group 2 (Near-certain): {(df['hla_proxy'] == 2).sum()}")

# ============================================================================
# CREATE LIPEDEMA PHENOTYPE
# ============================================================================
# Filter to women with complete DXA data
df_analysis = df.dropna(subset=["leg_trunk_ratio"]).copy()
print(f"\nWomen with complete DXA data: {len(df_analysis)}")

# Lipedema phenotype: leg-to-trunk ratio > 90th percentile
p90 = df_analysis["leg_trunk_ratio"].quantile(0.90)
df_analysis["lipedema_phenotype"] = (df_analysis["leg_trunk_ratio"] > p90).astype(int)
print(f"Lipedema phenotype (>P90): {df_analysis['lipedema_phenotype'].sum()} ({100*df_analysis['lipedema_phenotype'].mean():.1f}%)")

# ============================================================================
# ANALYSIS: HLA PROXY PREVALENCE BY LIPEDEMA PHENOTYPE
# ============================================================================
print("\n" + "="*70)
print("MAIN ANALYSIS: HLA Proxy Prevalence by Lipedema Phenotype")
print("="*70)

# Compare HLA proxy groups between lipedema and non-lipedema
lip = df_analysis[df_analysis["lipedema_phenotype"] == 1]
nonlip = df_analysis[df_analysis["lipedema_phenotype"] == 0]

print(f"\nLipedema phenotype group: n={len(lip)}")
print(f"Non-lipedema group: n={len(nonlip)}")

# Prevalence of HLA proxy groups
results = []

# Any HLA-related autoimmunity (Group 1 or 2)
lip_hla_any = (lip["hla_proxy"] >= 1).sum()
nonlip_hla_any = (nonlip["hla_proxy"] >= 1).sum()
lip_prev = 100 * lip_hla_any / len(lip)
nonlip_prev = 100 * nonlip_hla_any / len(nonlip)

# Fisher's exact test
table = [[lip_hla_any, len(lip) - lip_hla_any],
         [nonlip_hla_any, len(nonlip) - nonlip_hla_any]]
odds_ratio, p_value = stats.fisher_exact(table)

print(f"\nAny HLA-related autoimmunity (Group 1+2):")
print(f"  Lipedema: {lip_hla_any}/{len(lip)} ({lip_prev:.2f}%)")
print(f"  Non-lipedema: {nonlip_hla_any}/{len(nonlip)} ({nonlip_prev:.2f}%)")
print(f"  Odds Ratio: {odds_ratio:.3f}")
print(f"  P-value: {p_value:.4f}")

results.append({
    "Condition": "Any HLA-related autoimmunity",
    "Lipedema n": lip_hla_any,
    "Lipedema %": lip_prev,
    "Non-lipedema n": nonlip_hla_any,
    "Non-lipedema %": nonlip_prev,
    "OR": odds_ratio,
    "P-value": p_value
})

# Near-certain HLA+ (Group 2 - celiac)
lip_celiac = (lip["hla_proxy"] == 2).sum()
nonlip_celiac = (nonlip["hla_proxy"] == 2).sum()

if lip_celiac + nonlip_celiac > 0:
    lip_prev2 = 100 * lip_celiac / len(lip)
    nonlip_prev2 = 100 * nonlip_celiac / len(nonlip)

    table2 = [[lip_celiac, len(lip) - lip_celiac],
              [nonlip_celiac, len(nonlip) - nonlip_celiac]]
    or2, p2 = stats.fisher_exact(table2)

    print(f"\nCeliac serology+ (Group 2):")
    print(f"  Lipedema: {lip_celiac}/{len(lip)} ({lip_prev2:.2f}%)")
    print(f"  Non-lipedema: {nonlip_celiac}/{len(nonlip)} ({nonlip_prev2:.2f}%)")
    print(f"  Odds Ratio: {or2:.3f}")
    print(f"  P-value: {p2:.4f}")

    results.append({
        "Condition": "Celiac serology+ (near-certain HLA)",
        "Lipedema n": lip_celiac,
        "Lipedema %": lip_prev2,
        "Non-lipedema n": nonlip_celiac,
        "Non-lipedema %": nonlip_prev2,
        "OR": or2,
        "P-value": p2
    })

# High probability HLA+ (Group 1)
lip_g1 = (lip["hla_proxy"] == 1).sum()
nonlip_g1 = (nonlip["hla_proxy"] == 1).sum()

if lip_g1 + nonlip_g1 > 0:
    lip_prev1 = 100 * lip_g1 / len(lip)
    nonlip_prev1 = 100 * nonlip_g1 / len(nonlip)

    table1 = [[lip_g1, len(lip) - lip_g1],
              [nonlip_g1, len(nonlip) - nonlip_g1]]
    or1, p1 = stats.fisher_exact(table1)

    print(f"\nT1D/Autoimmune thyroiditis (Group 1):")
    print(f"  Lipedema: {lip_g1}/{len(lip)} ({lip_prev1:.2f}%)")
    print(f"  Non-lipedema: {nonlip_g1}/{len(nonlip)} ({nonlip_prev1:.2f}%)")
    print(f"  Odds Ratio: {or1:.3f}")
    print(f"  P-value: {p1:.4f}")

    results.append({
        "Condition": "T1D/Autoimmune thyroiditis",
        "Lipedema n": lip_g1,
        "Lipedema %": lip_prev1,
        "Non-lipedema n": nonlip_g1,
        "Non-lipedema %": nonlip_prev1,
        "OR": or1,
        "P-value": p1
    })

# Individual conditions
for cond, var in [("Type 1 Diabetes (probable)", "t1d_probable"),
                   ("Autoimmune Thyroiditis", "autoimmune_thyroid"),
                   ("Self-reported Thyroid Problem", "thyroid_problem")]:
    if var in df_analysis.columns:
        lip_n = lip[var].sum()
        nonlip_n = nonlip[var].sum()

        if lip_n + nonlip_n > 0:
            lip_pct = 100 * lip_n / len(lip)
            nonlip_pct = 100 * nonlip_n / len(nonlip)

            tbl = [[lip_n, len(lip) - lip_n],
                   [nonlip_n, len(nonlip) - nonlip_n]]
            o, p = stats.fisher_exact(tbl)

            print(f"\n{cond}:")
            print(f"  Lipedema: {lip_n}/{len(lip)} ({lip_pct:.2f}%)")
            print(f"  Non-lipedema: {nonlip_n}/{len(nonlip)} ({nonlip_pct:.2f}%)")
            print(f"  Odds Ratio: {o:.3f}")
            print(f"  P-value: {p:.4f}")

            results.append({
                "Condition": cond,
                "Lipedema n": lip_n,
                "Lipedema %": lip_pct,
                "Non-lipedema n": nonlip_n,
                "Non-lipedema %": nonlip_pct,
                "OR": o,
                "P-value": p
            })

# ============================================================================
# SENSITIVITY ANALYSIS: Different Lipedema Thresholds
# ============================================================================
print("\n" + "="*70)
print("SENSITIVITY ANALYSIS: Different Lipedema Phenotype Thresholds")
print("="*70)

sensitivity_results = []
for percentile in [75, 80, 85, 90, 95]:
    threshold = df_analysis["leg_trunk_ratio"].quantile(percentile/100)
    lip_group = df_analysis["leg_trunk_ratio"] > threshold

    # Any HLA
    lip_hla = (df_analysis.loc[lip_group, "hla_proxy"] >= 1).sum()
    nonlip_hla = (df_analysis.loc[~lip_group, "hla_proxy"] >= 1).sum()

    n_lip = lip_group.sum()
    n_nonlip = (~lip_group).sum()

    if n_lip > 0 and n_nonlip > 0:
        tbl = [[lip_hla, n_lip - lip_hla],
               [nonlip_hla, n_nonlip - nonlip_hla]]
        o, p = stats.fisher_exact(tbl)

        lip_pct = 100 * lip_hla / n_lip
        nonlip_pct = 100 * nonlip_hla / n_nonlip

        print(f"\n>P{percentile} (ratio > {threshold:.3f}):")
        print(f"  Lipedema: {lip_hla}/{n_lip} ({lip_pct:.2f}%)")
        print(f"  Non-lipedema: {nonlip_hla}/{n_nonlip} ({nonlip_pct:.2f}%)")
        print(f"  OR: {o:.3f}, p = {p:.4f}")

        sensitivity_results.append({
            "Percentile": percentile,
            "Threshold": threshold,
            "Lipedema n": n_lip,
            "Lipedema HLA+ %": lip_pct,
            "Non-lipedema HLA+ %": nonlip_pct,
            "OR": o,
            "P-value": p
        })

# ============================================================================
# SAVE RESULTS
# ============================================================================
print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

# Main results table
results_df = pd.DataFrame(results)
results_df.to_csv("tables/Table_HLA_Proxy_Analysis.csv", index=False)
print("Saved: tables/Table_HLA_Proxy_Analysis.csv")

# Sensitivity analysis table
sens_df = pd.DataFrame(sensitivity_results)
sens_df.to_csv("tables/Table_HLA_Sensitivity_Analysis.csv", index=False)
print("Saved: tables/Table_HLA_Sensitivity_Analysis.csv")

# ============================================================================
# CREATE FIGURE: HLA Proxy Prevalence by Lipedema Phenotype
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Prevalence comparison
conditions = ["Any HLA-related\nautoimmunity", "Celiac\n(near-certain)", "T1D/Thyroiditis\n(high probability)"]
lip_prev_list = []
nonlip_prev_list = []

for cond in ["Any HLA-related autoimmunity", "Celiac serology+ (near-certain HLA)", "T1D/Autoimmune thyroiditis"]:
    row = [r for r in results if r["Condition"] == cond]
    if row:
        lip_prev_list.append(row[0]["Lipedema %"])
        nonlip_prev_list.append(row[0]["Non-lipedema %"])
    else:
        lip_prev_list.append(0)
        nonlip_prev_list.append(0)

x = np.arange(len(conditions))
width = 0.35

bars1 = axes[0].bar(x - width/2, nonlip_prev_list, width, label='Non-lipedema', color='gray', alpha=0.7)
bars2 = axes[0].bar(x + width/2, lip_prev_list, width, label='Lipedema phenotype', color='#4CAF50', alpha=0.7)

axes[0].set_ylabel('Prevalence (%)')
axes[0].set_title('A. HLA Proxy Prevalence by Lipedema Phenotype')
axes[0].set_xticks(x)
axes[0].set_xticklabels(conditions, fontsize=9)
axes[0].legend()

# Add value labels
for bar in bars1:
    height = bar.get_height()
    if height > 0:
        axes[0].annotate(f'{height:.2f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)

for bar in bars2:
    height = bar.get_height()
    if height > 0:
        axes[0].annotate(f'{height:.2f}%',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)

# Panel B: Sensitivity analysis (OR across thresholds)
if len(sensitivity_results) > 0:
    percentiles = [r["Percentile"] for r in sensitivity_results]
    ors = [r["OR"] for r in sensitivity_results]

    axes[1].plot(percentiles, ors, 'o-', color='#2196F3', linewidth=2, markersize=8)
    axes[1].axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('Lipedema Phenotype Threshold (Percentile)')
    axes[1].set_ylabel('Odds Ratio')
    axes[1].set_title('B. Sensitivity Analysis: OR Across Thresholds')

    # Add p-value annotations
    for i, r in enumerate(sensitivity_results):
        p_text = f"p={r['P-value']:.3f}" if r['P-value'] >= 0.001 else "p<0.001"
        axes[1].annotate(p_text,
                        xy=(r["Percentile"], r["OR"]),
                        xytext=(0, 10), textcoords="offset points",
                        ha='center', fontsize=8)

plt.tight_layout()
plt.savefig("figures/Figure_HLA_Proxy_Analysis.png", dpi=300, bbox_inches='tight')
print("Saved: figures/Figure_HLA_Proxy_Analysis.png")
plt.close()

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print("""
HLA PROXY CONSTRUCTION:
- Group 2 (near-certain DQ2/DQ8): Serologically confirmed celiac disease
- Group 1 (high probability): T1D probable OR autoimmune thyroiditis (with TPO+)
- Group 0 (reference): Lower probability

KEY FINDINGS:
""")

for r in results[:3]:
    direction = "LOWER" if r["OR"] < 1 else "HIGHER" if r["OR"] > 1 else "SIMILAR"
    sig = "*" if r["P-value"] < 0.05 else ""
    print(f"- {r['Condition']}: {direction} in lipedema phenotype")
    print(f"  OR={r['OR']:.3f}, p={r['P-value']:.4f}{sig}")

print("""
INTERPRETATION:
Since HLA-DQ2/DQ8 genotyping is not available in NHANES, we constructed a
proxy for 'HLA-related autoimmunity background' based on serologically
defined celiac disease, type 1 diabetes, and autoimmune thyroid disease.

Self-reported celiac disease without serological support was NOT used alone,
given the high likelihood of misclassification with non-celiac gluten
sensitivity (NCGS).
""")

# Save analysis dataset
df_analysis.to_csv("nhanes_hla_proxy_dataset.csv", index=False)
print("\nSaved: nhanes_hla_proxy_dataset.csv")

print("\n" + "="*70)
print("ANALYSIS COMPLETE")
print("="*70)
