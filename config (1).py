"""
Central configuration for the Hepatic Diabesity Transcriptomics Pipeline.

Defines dataset accessions, group label harmonization, molecular switch panel,
color scheme, and output paths shared across all pipeline modules.
"""
import os

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.abspath(__file__))

DATA_DIR = os.path.join(ROOT, "data")
OUT_DIR = os.path.join(ROOT, "outputs")
FIG_DIR = os.path.join(OUT_DIR, "figures")
TAB_DIR = os.path.join(OUT_DIR, "tables")
RES_DIR = os.path.join(OUT_DIR, "results")

for _d in (DATA_DIR, OUT_DIR, FIG_DIR, TAB_DIR, RES_DIR):
    os.makedirs(_d, exist_ok=True)

# Local copies of the GEO SOFT records / supplementary files (no network access needed)
LOCAL_GEO_DIR = os.path.join(ROOT, "hepato_pipeline_data")

# ----------------------------------------------------------------------------
# GEO datasets
# ----------------------------------------------------------------------------
# role: 'discovery' -> harmonized n=103 pool ; 'validation' -> independent cohort
DATASETS = {
    "GSE121344": {"platform": "rnaseq",     "gpl": "GPL20301", "role": "discovery",  "n": 12},
    "GSE15653":  {"platform": "microarray", "gpl": "GPL96",    "role": "discovery",  "n": 18},
    "GSE48452":  {"platform": "microarray", "gpl": "GPL11532", "role": "discovery",  "n": 73},
    "GSE64998":  {"platform": "microarray", "gpl": "GPL11532", "role": "validation", "n": 21},
}

DISCOVERY_SETS = [k for k, v in DATASETS.items() if v["role"] == "discovery"]
VALIDATION_SET = "GSE64998"

# Canonical group labels used throughout the pipeline
GROUPS = ["Lean", "Obese ND", "Obese T2D"]
GROUP_ORDER = {"Lean": 0, "Obese ND": 1, "Obese T2D": 2}

# ----------------------------------------------------------------------------
# Molecular switch panel (target regulators from the manuscript)
# ----------------------------------------------------------------------------
SWITCH_GENES = ["OSMR", "STAT3", "NFKB1", "PPARG", "TNF",
                "IL6", "FOXO1", "SIRT1", "YAP1", "TGFB1"]

# Expected directionality along lean->obese->T2D axis (for reconstruction / QC).
# +1 = up with progression, -1 = down with progression
SWITCH_DIRECTION = {
    "OSMR": +1, "STAT3": +1, "NFKB1": +1, "TNF": +1, "IL6": +1,
    "TGFB1": +1, "YAP1": +1, "FOXO1": +1,
    "PPARG": -1, "SIRT1": -1,
}

# ----------------------------------------------------------------------------
# Manuscript colour scheme
# ----------------------------------------------------------------------------
COLORS = {"Lean": "#27ae60", "Obese ND": "#2980b9", "Obese T2D": "#c0392b"}
ACCENT = {"discovery": "#e67e22", "validation": "#8e44ad", "diag": "#7f8c8d"}

# Figure defaults
DPI = 300
FONT_FAMILY = ["Arial", "Helvetica", "DejaVu Sans"]

# Analysis thresholds
DEG_LOG2FC = 1.0
DEG_ALPHA = 0.05
RANDOM_STATE = 42
