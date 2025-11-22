# Gynoid Fat Distribution and Celiac Disease

Exploratory analysis investigating the relationship between gynoid fat distribution (lipedema-compatible phenotype) and celiac disease prevalence using NHANES 2011-2014 data.

## Citation

If you use this code, please cite our manuscript:
> [Citation will be added upon publication]

## Overview

This repository contains the analysis code for our study examining whether women with higher gynoid (lower-body) fat distribution have different celiac disease prevalence. We use validated body composition measures from dual-energy X-ray absorptiometry (DXA) as proxies for lipedema phenotype.

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

## Requirements

```bash
pip install pandas numpy requests statsmodels scipy matplotlib seaborn
```

## Usage

Run the main analysis:

```bash
python nhanes_weighted_final_analysis.py
```

This will:
1. Download NHANES data files from CDC (cached locally in `nhanes_cache/`)
2. Process and merge datasets
3. Calculate lipedema phenotype proxies
4. Perform survey-weighted statistical analyses
5. Generate figures and tables

For HLA proxy sensitivity analysis:

```bash
python nhanes_hla_proxy_analysis.py
```

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
├── tables/                             # Generated tables (CSV)
├── manuscript_final.md                 # Methods and results
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
- All analyses incorporate NHANES complex survey design weights
- Survey-weighted t-tests for continuous variables
- Fisher's exact test for categorical variables (small cell counts)
- Sensitivity analyses stratified by BMI category

## License

This code is released under the MIT License. NHANES data is in the public domain.

## Contact

For questions about this analysis, please open an issue on this repository.
