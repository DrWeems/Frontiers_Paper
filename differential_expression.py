"""
differential_expression.py -- DEG analysis across the three group contrasts.

Contrasts:
    Lean vs Obese ND
    Obese ND vs Obese T2D
    Lean vs Obese T2D

For each gene we compute log2 fold-change (groupB - groupA on the log2 scale),
a Mann-Whitney U p-value (robust, non-parametric) and Benjamini-Hochberg FDR.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import config as C


CONTRASTS = [
    ("Lean", "Obese ND"),
    ("Obese ND", "Obese T2D"),
    ("Lean", "Obese T2D"),
]


def _one_contrast(expr, meta, a, b):
    sa = meta.index[meta["group"] == a]
    sb = meta.index[meta["group"] == b]
    Xa = expr[sa].values
    Xb = expr[sb].values

    log2fc = Xb.mean(axis=1) - Xa.mean(axis=1)
    pvals = np.ones(expr.shape[0])
    for i in range(expr.shape[0]):
        try:
            _, p = stats.mannwhitneyu(Xa[i], Xb[i], alternative="two-sided")
        except ValueError:
            p = 1.0
        pvals[i] = p
    padj = multipletests(pvals, method="fdr_bh")[1]

    df = pd.DataFrame({
        "gene": expr.index,
        "contrast": f"{a} vs {b}",
        "log2FC": log2fc,
        "pvalue": pvals,
        "padj": padj,
    })
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["significant"] = (df["padj"] < C.DEG_ALPHA) & (df["log2FC"].abs() >= C.DEG_LOG2FC)
    return df


def run(expr, meta, cache=True):
    frames = [_one_contrast(expr, meta, a, b) for a, b in CONTRASTS]
    degs = pd.concat(frames, ignore_index=True)
    degs = degs.sort_values(["contrast", "padj"])
    n_sig = degs.groupby("contrast")["significant"].sum().to_dict()
    print(f"  [DEG] significant genes per contrast: {n_sig}", flush=True)
    if cache:
        degs.to_csv(os.path.join(C.TAB_DIR, "table_degs.csv"), index=False)
    return degs
