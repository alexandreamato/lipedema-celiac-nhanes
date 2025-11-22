# Gynoid Fat Distribution and Celiac Disease: NHANES Analysis

This repository contains the analytic code and supplementary materials for the research article:

**"Gynoid Fat Distribution and Celiac Disease: An Exploratory Analysis of NHANES 2011-2014"**

## Overview

This study investigates the relationship between gynoid fat distribution (lipedema-compatible phenotype) and celiac disease prevalence using data from the National Health and Nutrition Examination Survey (NHANES) 2011-2014.

The analysis employs complex survey design weighting to provide nationally representative estimates for the US adult female population.

**Key findings:**
- Women with celiac disease had significantly lower gynoid region percent fat (39.5% vs 42.6%, p=0.0007)
- This association persisted when restricted to overweight/obese women, ruling out leanness bias
- Women with lipedema phenotype showed favorable immunometabolic profiles (lower inflammation, lower insulin resistance)

## Data Source

This analysis uses publicly available data from the [National Health and Nutrition Examination Survey (NHANES)](https://www.cdc.gov/nchs/nhanes/) 2011-2014 cycles. The analysis script automatically downloads required data files from the CDC website.

**NHANES modules used:**
- Demographics (DEMO)
- Body Measures (BMX)
- Dual-Energy X-ray Absorptiometry (DXX)
- Tissue Transglutaminase & Endomysial Antibody (TGEMA)
- Complete Blood Count (CBC)
- Plasma Fasting Glucose & Insulin (GLU)
- Diabetes (DIQ)
- Thyroid Profile (THYROD)

## How to Run the Analysis

The script is designed to automatically download the required raw `.xpt` (SAS transport) files directly from the CDC/NHANES website, ensuring reproducibility without hosting restricted data.

### Prerequisites

Python 3.8+ installed.

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/alexandreamato/lipedema-celiac-nhanes.git
   cd lipedema-celiac-nhanes
   ```

2. Install required packages:
   ```bash
   pip install pandas numpy requests statsmodels scipy matplotlib seaborn
   ```

3. Run the main analysis:
   ```bash
   python nhanes_weighted_final_analysis.py
   ```

4. For HLA proxy sensitivity analysis:
   ```bash
   python nhanes_hla_proxy_analysis.py
   ```

This will:
1. Download NHANES data files from CDC (cached locally in `nhanes_cache/`)
2. Process and merge datasets
3. Calculate lipedema phenotype proxies
4. Perform survey-weighted statistical analyses
5. Generate figures and tables

## Repository Structure

```
.
├── nhanes_weighted_final_analysis.py   # Main analysis script
├── nhanes_hla_proxy_analysis.py        # HLA proxy sensitivity analysis
├── figures/                            # Generated figures
│   ├── Figure1_Gynoid_Fat.png
│   ├── Figure2_Multiple_Proxies.png
│   ├── Figure3_Celiac_Distribution.png
│   ├── Figure4_Immunometabolic.png
│   ├── FigureS1_Quartiles.png
│   └── FigureS2_Leanness_Bias.png
└── README.md
```

## Key Variables

### Celiac Disease Definition
Strict dual-positive criterion:
- tTG-IgA positive (LBXTTG = 1) **AND**
- EMA-IgA positive (LBXEMA = 1)

### Lipedema Phenotype Proxies
1. **Gynoid region percent fat** (primary) - from DXA
2. **Leg-to-trunk fat ratio** - leg fat mass / trunk fat mass
3. **Lipedema phenotype** - leg-to-trunk ratio > 90th percentile

### Statistical Approach
- All analyses incorporate NHANES complex survey design weights (WTMEC2YR adjusted for combined cycles)
- Survey-weighted t-tests for continuous variables
- Fisher's exact test for categorical variables (small cell counts)
- Sensitivity analyses stratified by BMI category

## Citation

If you use this code or data structure in your research, please cite the original paper:

> [Citation will be added upon publication]

## Contact

For questions regarding the code or methodology:

- **Author:** Alexandre Campos Moraes Amato
- **Institution:** Amato - Instituto de Medicina Avancada

---

*Disclaimer: This repository is for research reproducibility purposes. The raw data is property of the CDC/NCHS.*
