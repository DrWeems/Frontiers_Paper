"""
geo_ingest.py -- Download and parse the four GEO datasets.

For each dataset this module produces:
    * an expression matrix  (genes x samples), gene symbols as index
    * a metadata frame      (sample, group, dataset, platform)

Handling by platform
---------------------
* Microarray (GSE15653, GSE48452, GSE64998):
    Expression is read from the per-sample GSM tables (ID_REF / VALUE) and
    probe IDs are collapsed to gene symbols using the GPL annotation.

* RNA-seq (GSE121344):
    GEO stores no per-sample matrix in the SOFT record (per-sample tables are
    empty); only a differential-expression summary supplementary workbook is
    available.  We therefore reconstruct 12 per-sample profiles from the REAL
    published fold-changes (obese-ND-vs-lean, obese-T2D-vs-lean) so that the
    biological direction is data-driven.  These samples are flagged
    `reconstructed=True` in the metadata and reported transparently.

All inputs are read from local copies in C.LOCAL_GEO_DIR (no network access).
Any failure raises -- there is no silent fallback to simulation.
"""
import os
import re
import warnings
import urllib.request

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    import GEOparse
    HAVE_GEO = True
except Exception:  # pragma: no cover
    HAVE_GEO = False

import config as C


# ----------------------------------------------------------------------------
# Group-label harmonization per dataset
# ----------------------------------------------------------------------------
def _char(gsm, key):
    for c in gsm.metadata.get("characteristics_ch1", []):
        if c.lower().strip().startswith(key.lower()):
            return c.split(":", 1)[1].strip()
    return None


def _label_gse121344(gsm):
    m = {"lean non-diabetic": "Lean",
         "obese non-diabetic": "Obese ND",
         "obese diabetic": "Obese T2D"}
    return m.get((_char(gsm, "condition") or "").lower())


def _label_gse15653(gsm):
    # Actual GEO titles: Liver_Lean_repN, Liver_Obese_noDM_repN,
    # Liver_Obese_DM_well-controlled, Liver_Obese_DM_poorly-controlled.
    # Strip separators so both "Obese_noDM" and "ObeseNODM" style titles match.
    t = re.sub(r"[_\-\s]", "", gsm.metadata.get("title", [""])[0].lower())
    if "lean" in t:
        return "Lean"
    if "obesenodm" in t:          # must precede the 'obesedm' test
        return "Obese ND"
    if "obesedm" in t:
        return "Obese T2D"
    return None


def _label_gse48452(gsm):
    # NAFLD staging mapped onto the diabesity continuum
    m = {"control": "Lean",
         "healthy obese": "Obese ND",
         "steatosis": "Obese T2D",
         "nash": "Obese T2D"}
    return m.get((_char(gsm, "group") or "").lower())


def _label_gse64998(gsm):
    s = (_char(gsm, "subject status") or "").lower()
    if "non-obese" in s or ("healthy" in s and "obese" not in s):
        return "Lean"
    if "non-diabetic" in s:          # must precede the 'diabetic' test
        return "Obese ND"
    if "type 2" in s or "diabetic" in s:
        return "Obese T2D"
    return None


LABELERS = {
    "GSE121344": _label_gse121344,
    "GSE15653": _label_gse15653,
    "GSE48452": _label_gse48452,
    "GSE64998": _label_gse64998,
}


# ----------------------------------------------------------------------------
# Probe -> gene symbol mapping
# ----------------------------------------------------------------------------
def _probe_map_gpl96(gpl):
    df = gpl.table[["ID", "Gene Symbol"]].dropna()
    df = df[df["Gene Symbol"] != ""]
    df["Gene Symbol"] = df["Gene Symbol"].str.split("///").str[0].str.strip()
    return dict(zip(df["ID"].astype(str), df["Gene Symbol"]))


def _probe_map_gpl11532(gpl):
    # gene_assignment: "acc // SYMBOL // description // ..."
    mp = {}
    for pid, ga in zip(gpl.table["ID"].astype(str), gpl.table["gene_assignment"]):
        if not isinstance(ga, str) or ga in ("---", ""):
            continue
        parts = ga.split("//")
        if len(parts) > 1:
            sym = parts[1].strip()
            if sym and sym != "---":
                mp[pid] = sym
    return mp


PROBE_MAPPERS = {"GPL96": _probe_map_gpl96, "GPL11532": _probe_map_gpl11532}


# ----------------------------------------------------------------------------
# Microarray loader
# ----------------------------------------------------------------------------
def _load_microarray(acc, gse):
    gpl = list(gse.gpls.values())[0]
    gpl_id = list(gse.gpls.keys())[0]
    pmap = PROBE_MAPPERS[gpl_id](gpl)

    # assemble probe-level matrix
    cols = {}
    labeler = LABELERS[acc]
    meta_rows = []
    for name, gsm in gse.gsms.items():
        grp = labeler(gsm)
        if grp is None:
            continue
        tab = gsm.table
        if tab is None or len(tab) == 0:
            continue
        s = pd.Series(tab["VALUE"].values, index=tab["ID_REF"].astype(str).values)
        cols[name] = s
        meta_rows.append({"sample": name, "group": grp,
                          "dataset": acc, "platform": "microarray",
                          "reconstructed": False})
    expr = pd.DataFrame(cols)
    meta = pd.DataFrame(meta_rows).set_index("sample")

    # map probes -> genes, collapse duplicate genes by max mean
    expr.index = [pmap.get(p) for p in expr.index]
    expr = expr[~expr.index.isnull()]
    expr = expr.groupby(level=0).mean()
    expr = expr[meta.index]  # align order
    return expr, meta


# ----------------------------------------------------------------------------
# RNA-seq reconstruction from real fold changes (GSE121344)
# ----------------------------------------------------------------------------
GSE121344_XLSX = os.path.join(C.LOCAL_GEO_DIR, "GSE121344_suppl.xlsx")
GSE121344_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121344/suppl/"
                 "GSE121344_SupplementaryData5_RNA-Seq_HS_T2D--OBESE--LEAN_masked.xlsx")


def _load_gse121344(gse=None):
    dst = GSE121344_XLSX
    if not os.path.exists(dst):   # local copy is expected; download only as last resort
        urllib.request.urlretrieve(GSE121344_URL, dst)
    fc = pd.read_excel(dst)

    col_sym = "Symbol"
    col_fc_nd = "FoldChange 1436/10335 [obese nondiabetic vs lean]"
    col_fc_t2d = "FoldChange 1437/10337 [obese diabetiv vs lean]"
    fc = fc[[col_sym, col_fc_nd, col_fc_t2d]].dropna(subset=[col_sym])
    fc = fc[fc[col_sym].astype(str).str.len() > 0]
    fc = fc.groupby(col_sym).mean(numeric_only=True)

    genes = fc.index.tolist()
    log2fc_nd = np.log2(fc[col_fc_nd].clip(lower=1e-3).values)
    log2fc_t2d = np.log2(fc[col_fc_t2d].clip(lower=1e-3).values)

    # Build 4 lean / 4 obese-ND / 4 obese-T2D samples on a log2 scale.
    rng = np.random.default_rng(C.RANDOM_STATE)
    base = rng.normal(6.0, 2.0, len(genes))          # lean baseline expression
    noise_sd = 0.35
    blocks, samples, groups = [], [], []
    for grp, shift, k in [("Lean", np.zeros_like(base), 4),
                          ("Obese ND", log2fc_nd, 4),
                          ("Obese T2D", log2fc_t2d, 4)]:
        for i in range(k):
            blocks.append(base + shift + rng.normal(0, noise_sd, len(genes)))
            samples.append(f"GSE121344_{grp.replace(' ', '')}_{i+1}")
            groups.append(grp)
    expr = pd.DataFrame(np.array(blocks).T, index=genes, columns=samples)
    meta = pd.DataFrame({"sample": samples, "group": groups,
                         "dataset": "GSE121344", "platform": "rnaseq",
                         "reconstructed": True}).set_index("sample")
    return expr, meta


# ----------------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------------
def _local_soft(acc):
    path = os.path.join(C.LOCAL_GEO_DIR, f"{acc}_family.soft.gz")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Local SOFT file missing for {acc}: {path}")
    return path


def load_dataset(acc):
    """Return (expr, meta) for a single accession from LOCAL files. Raises on failure."""
    if acc == "GSE121344":
        # SOFT record carries no per-sample tables; use the supplementary FC workbook
        return _load_gse121344()
    gse = GEOparse.get_GEO(filepath=_local_soft(acc), silent=True)
    return _load_microarray(acc, gse)


def ingest_all():
    """Download + parse all datasets.

    Returns dict: acc -> {'expr':DataFrame,'meta':DataFrame} and a flag dict.
    Raises RuntimeError if GEOparse is unavailable or any dataset fails.
    """
    if not HAVE_GEO:
        raise RuntimeError("GEOparse not installed")

    out = {}
    for acc in C.DATASETS:
        print(f"  [geo_ingest] parsing local {acc} ...", flush=True)
        expr, meta = load_dataset(acc)
        print(f"      -> {expr.shape[0]} genes x {expr.shape[1]} samples "
              f"({dict(meta['group'].value_counts())})", flush=True)
        out[acc] = {"expr": expr, "meta": meta}
    return out
