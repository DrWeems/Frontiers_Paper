"""
network_switches.py -- Molecular switch prioritization (Table II).

For each target regulator we correlate its harmonized expression with the
lean -> obese -> T2D progression axis (0/1/2) across the discovery pool.
Both Spearman and Pearson coefficients are reported; genes are ranked by the
strength of the Spearman association and FDR-corrected.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import config as C


def run(expr, meta, cache=True):
    axis = meta["group"].map(C.GROUP_ORDER).astype(float)
    rows = []
    for gene in C.SWITCH_GENES:
        if gene not in expr.index:
            rows.append({"gene": gene, "present": False})
            continue
        v = expr.loc[gene, meta.index].astype(float).values
        rho, p_s = stats.spearmanr(axis.values, v)
        r, p_p = stats.pearsonr(axis.values, v)
        rows.append({
            "gene": gene, "present": True,
            "spearman_r": rho, "spearman_p": p_s,
            "pearson_r": r, "pearson_p": p_p,
            "expected_direction": C.SWITCH_DIRECTION.get(gene, np.nan),
        })
    tab = pd.DataFrame(rows)

    present = tab["present"]
    tab.loc[present, "FDR"] = multipletests(
        tab.loc[present, "spearman_p"].values, method="fdr_bh")[1]
    tab["concordant"] = np.sign(tab["spearman_r"]) == tab["expected_direction"]
    tab = tab.sort_values("spearman_p", na_position="last").reset_index(drop=True)

    n_sig = int((tab["FDR"] < 0.05).sum())
    print(f"  [switches] {n_sig}/{present.sum()} switches significant "
          f"(FDR<0.05) along progression axis", flush=True)
    if cache:
        tab.to_csv(os.path.join(C.TAB_DIR, "table_switches.csv"), index=False)
    return tab
