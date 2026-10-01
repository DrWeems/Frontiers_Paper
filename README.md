# Frontiers_Paper
# HepatoSwitch: A Multi-Platform Transcriptomic Framework for Hepatic Diabesity

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational workflow, data harmonization logic, and statistical analysis scripts for the study:

> **"A Multi-Platform Transcriptomic Framework for Molecular Staging and Regulatory Prioritization in Hepatic Diabesity."**  
> *Ebony Weems-Oluremi et al.*  
> Submitted to *Frontiers in Endocrinology*.

---

## Overview

HepatoSwitch integrates transcriptomic profiling across multiple discovery cohorts ($n=103$) and an independent validation cohort ($n=21$) to characterize the metabolic and transcriptional continuum from lean liver to obese non-diabetic and obese type 2 diabetic (T2D) states.

The pipeline automates:
1. **GEO Data Ingestion**: Direct retrieval and extraction of multi-platform microarray and RNA-seq accessions.
2. **Harmonization & Batch Correction**: Inter-platform quantile standardization and Empirical Bayes batch harmonization via ComBat.
3. **Differential Expression Analysis**: Multi-contrast statistical evaluation (Mann-Whitney U with Benjamini-Hochberg FDR correction).
4. **Pathway Enrichment**: Pre-ranked Gene Set Enrichment Analysis (GSEA) and Over-Representation Analysis (ORA) against MSigDB Hallmark and KEGG pathways, with an automatic offline hypergeometric fallback.
5. **Molecular Switch Prioritization**: Non-parametric progression-axis correlation analysis prioritizing key regulatory nodes (notably *SIRT1* suppression and *TGFB1* activation).
6. **Machine Learning Classifier**: Leave-One-Out Cross-Validation (LOOCV) Random Forest staging model evaluated against the independent validation cohort ($n=21$).
7. **Publication Figures**: Automated generation of 6 publication-ready figures (300 DPI) adhering strictly to journal styling guidelines.

---

## Repository Structure

```text
.
├── config.py                     # Global paths, thresholds, gene panels, color schemes
├── geo_ingest.py                 # Download and parsing of GEO datasets (GEOparse)
├── preprocess.py                 # Multi-platform normalization and ComBat batch correction
├── differential_expression.py    # Non-parametric DEG testing across 3 disease contrasts
├── pathway_enrichment.py         # GSEA and Enrichr ORA (with offline fallback)
├── network_switches.py           # Molecular switch prioritization along the disease axis
├── ml_classifier.py              # Random Forest classifier with LOOCV & external validation
├── figures.py                    # Generation of Figures 1–6 (300 DPI)
├── simulate.py                   # Verified synthetic data generator for testing/CI
├── Frontiers_Pipeline_Corrected.ipynb  # Primary interactive pipeline orchestrator
├── requirements.txt              # Pinned Python package dependencies
└── README.md                     # Documentation and reproduction guide
```

---

## Datasets

The framework integrates four NCBI Gene Expression Omnibus (GEO) datasets:

| Dataset Accession | Platform | Platform ID | Role in Study | Sample Count ($n$) | Phenotypes |
| :--- | :--- | :--- | :--- | :---: | :--- |
| **GSE121344** | RNA-seq | GPL20301 | Discovery Pool | 12 | Lean, Obese ND, Obese T2D |
| **GSE15653** | Microarray | GPL96 | Discovery Pool | 18 | Lean, Obese ND, Obese T2D |
| **GSE48452** | Microarray | GPL11532 | Discovery Pool | 73 | Control, Healthy Obese, Steatosis, NASH |
| **GSE64998** | Microarray | GPL11532 | Independent Validation | 21 | Non-obese, Obese ND, Obese T2D |

---

## Installation & Setup

```bash
git clone https://github.com/EbonyWeems/HepatoSwitch-Transcriptomics.git
cd HepatoSwitch-Transcriptomics

python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

pip install --upgrade pip
pip install -r requirements.txt
```

---

## Running the Pipeline

Launch the orchestrator notebook:
```bash
jupyter notebook Frontiers_Pipeline_Corrected.ipynb
```
Select **Kernel -> Restart & Run All**.

---

## License

This project is licensed under the MIT License.
