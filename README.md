# Lipedema Phenotype vs. Celiac Disease: NHANES Analysis

This repository contains the analytic code, data dictionary, and supplementary materials for the research article:

**"The Immunological Shield Hypothesis: Phenotypic Divergence and Immunometabolic Differences Between Lipedema and Celiac Disease in NHANES"**

## 📄 Overview

This study investigates the relationship between the lipedema phenotype (defined by DXA-derived leg-to-trunk fat ratio) and celiac disease autoimmunity using data from the National Health and Nutrition Examination Survey (NHANES) 2011–2014.

The analysis employs complex survey design weighting to provide nationally representative estimates for the US adult female population.

## 📂 Repository Contents

* **`analysis_script.py`**: The complete Python script used to download NHANES data, process variables, apply survey weights, and perform statistical testing (T-tests, Fisher's Exact, Regression models).
* **`data_dictionary.md`**: Definitions of all variables created or used in the study (e.g., how "Lipedema Phenotype" or "Celiac Disease" were operationally defined).
* **`requirements.txt`**: List of Python libraries required to run the script.
* **`Supplementary_Material.pdf`**: The supplementary figures and tables referenced in the manuscript (Sensitivity Analysis & Secondary Analysis).

## 🚀 How to Run the Analysis

The script is designed to automatically download the required raw `.xpt` (SAS transport) files directly from the CDC/NHANES website, ensuring reproducibility without hosting restricted data.

### Prerequisites

You need Python 3.8+ installed.

### Installation

1.  Clone this repository:
    ```bash
    git clone https://github.com/alexandreamato/lipedema-celiac-nhanes.git
    cd lipedema-celiac-nhanes
    ```

2.  Install required packages:
    ```bash
    pip install pandas numpy scipy statsmodels matplotlib requests
    ```

3.  Run the script:
    ```bash
    python analysis_script.py
    ```

## 📊 Methodology Notes

* **Data Source:** [NHANES 2011-2014](https://wwwn.cdc.gov/nchs/nhanes/) (CDC/NCHS).
* **Weighting:** All primary analyses utilize `WTMEC2YR` (Mobile Examination Center weights) adjusted for the combined 4-year cycle.
* **Statistical Software:** Analysis performed using Python (Pandas, SciPy, Statsmodels).

## 📝 Citation

If you use this code or data structure in your research, please cite the original paper:

> [Insert Citation Here once published: Author Name, Title, Journal, Year, DOI]

## 📬 Contact

For questions regarding the code or methodology:

* **Author:** Alexandre Campos Moraes Amato
* **Institution:** Amato - Instituto de Medicina Avançada
 

---

*Disclaimer: This repository is for research reproducibility purposes. The raw data is property of the CDC/NCHS.*
