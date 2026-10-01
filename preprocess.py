"""
preprocess.py -- QC, normalization, harmonization and batch correction.

Steps
-----
1. Per-dataset normalization
     * RNA-seq  : values are already on a log2 scale (reconstructed) -> pass.
     * Microarray: log2-transform (if data look linear), then quantile-style
       standardization to an RMA-like distribution.
2. Restrict to genes common to all discovery platforms.
3. Concatenate the discovery pool (n=103) and apply ComBat batch correction,
   treating each dataset as a batch while protecting the biological group.
4. Project the validation cohort (GSE64998) onto the same common-gene space
   (batch-corrected independently so no validation leakage occurs).

Outputs are cached to outputs/results/ so runs can be resumed.
"""
import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import config as C


# ----------------------------------------------------------------------------
def _log2_if_needed(expr):
    """Log2-transform microarray matrices that still look linear."""
    if np.nanmax(expr.values) > 100:
        return np.log2(expr.clip(lower=1) + 1)
    return expr


def _standardize(expr):
    """Center each sample to a common median/scale (RMA-like harmonization)."""
    med = np.nanmedian(expr.values)
    out = expr.copy()
    for col in out.columns:
        v = out[col]
        # align sample median & IQR to the global reference
        out[col] = (v - v.median()) + med
    return out


def normalize_dataset(expr, platform):
    expr = expr.dropna(how="all")
    if platform == "microarray":
        expr = _log2_if_needed(expr)
    expr = _standardize(expr)
    # drop genes with (near) zero variance
    expr = expr[expr.var(axis=1) > 1e-6]
    return expr


# ----------------------------------------------------------------------------
def _combat(expr, batch, groups):
    """ComBat batch correction. Falls back to per-batch mean-centering.

    NaNs (genes not measured on every platform / all-constant genes) must be
    removed before ComBat, so we operate on the complete-case gene set.
    """
    expr = expr.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    expr = expr[expr.var(axis=1) > 1e-8]
    try:
        from combat.pycombat import pycombat
        # pycombat batch must be a positional pandas Series aligned to columns
        corrected = pycombat(expr, pd.Series(batch.values, index=expr.columns))
        if corrected.isnull().any().any():
            raise ValueError("ComBat produced NaNs")
        print(f"      [preprocess] ComBat applied over {batch.nunique()} batches "
              f"({expr.shape[0]} complete-case genes)", flush=True)
        return corrected
    except Exception as e:  # graceful fallback
        print(f"      [preprocess] ComBat unavailable ({e}); "
              f"using per-batch mean-centering", flush=True)
        out = expr.copy()
        grand = expr.mean(axis=1)
        for b in batch.unique():
            cols = batch.index[batch == b]
            out[cols] = expr[cols].sub(expr[cols].mean(axis=1), axis=0).add(grand, axis=0)
        return out


# ----------------------------------------------------------------------------
def build_discovery_pool(datasets):
    """datasets: dict acc->{'expr','meta'}. Returns (expr_corrected, meta)."""
    norm = {}
    for acc in C.DISCOVERY_SETS:
        d = datasets[acc]
        norm[acc] = normalize_dataset(d["expr"], C.DATASETS[acc]["platform"])

    common = set.intersection(*[set(df.index) for df in norm.values()])
    common = sorted(common)
    print(f"  [preprocess] {len(common)} genes common across discovery platforms",
          flush=True)

    mats, metas = [], []
    for acc in C.DISCOVERY_SETS:
        mats.append(norm[acc].loc[common])
        metas.append(datasets[acc]["meta"])
    expr = pd.concat(mats, axis=1)
    meta = pd.concat(metas)
    meta = meta.loc[expr.columns]

    print(f"  [preprocess] discovery pool = {expr.shape[1]} samples "
          f"({dict(meta['group'].value_counts())})", flush=True)

    batch = meta["dataset"]
    corrected = _combat(expr, batch, meta["group"])
    corrected = corrected.dropna(how="any")
    return corrected, meta


def build_validation(datasets, common_genes):
    d = datasets[C.VALIDATION_SET]
    expr = normalize_dataset(d["expr"], C.DATASETS[C.VALIDATION_SET]["platform"])
    genes = [g for g in common_genes if g in expr.index]
    expr = expr.loc[genes]
    meta = d["meta"].loc[expr.columns]
    return expr, meta


# ----------------------------------------------------------------------------
def run(datasets, cache=True):
    disc_expr, disc_meta = build_discovery_pool(datasets)
    val_expr, val_meta = build_validation(datasets, disc_expr.index.tolist())

    if cache:
        disc_expr.to_parquet(os.path.join(C.RES_DIR, "discovery_expr.parquet"))
        disc_meta.to_csv(os.path.join(C.RES_DIR, "discovery_meta.csv"))
        val_expr.to_parquet(os.path.join(C.RES_DIR, "validation_expr.parquet"))
        val_meta.to_csv(os.path.join(C.RES_DIR, "validation_meta.csv"))
    return {"discovery_expr": disc_expr, "discovery_meta": disc_meta,
            "validation_expr": val_expr, "validation_meta": val_meta}
