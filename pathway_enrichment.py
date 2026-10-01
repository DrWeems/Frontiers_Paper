"""
pathway_enrichment.py -- GSEA (pre-ranked) + ORA using gseapy.

* Pre-ranked GSEA on the Lean-vs-Obese-T2D contrast, ranked by
  sign(log2FC) * -log10(pvalue), against MSigDB Hallmark and KEGG gene sets.
* Over-representation analysis (Enrichr) on the significant DEG list.

If gseapy cannot reach Enrichr / download gene sets (offline), a transparent
fallback computes a simple hypergeometric ORA against a small built-in set of
metabolic / inflammatory gene sets so the pipeline still yields a table.
"""
import os
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import config as C

GENE_SETS = ["MSigDB_Hallmark_2020", "KEGG_2021_Human"]

# Minimal offline fallback gene sets (curated, metabolism / inflammation)
_FALLBACK_SETS = {
    "Inflammatory_Response": ["TNF", "IL6", "NFKB1", "STAT3", "OSMR", "CXCL8",
                              "IL1B", "CCL2", "SOCS3", "TLR4"],
    "Insulin_Signaling": ["FOXO1", "PPARG", "IRS1", "IRS2", "AKT1", "PIK3R1",
                          "SLC2A4", "PDK4", "PCK1", "G6PC"],
    "Lipid_Metabolism": ["PPARG", "SREBF1", "FASN", "SCD", "CD36", "ACACA",
                         "CPT1A", "PLIN2", "FABP1", "APOB"],
    "Fibrosis_TGFb": ["TGFB1", "COL1A1", "COL1A2", "ACTA2", "TIMP1", "SMAD3",
                      "CTGF", "YAP1", "SERPINE1", "MMP2"],
    "Longevity_Sirtuin": ["SIRT1", "PGC1A", "PPARGC1A", "NAMPT", "FOXO1",
                          "TP53", "NFE2L2", "SOD2", "CAT", "PRKAA1"],
}


def _preranked(degs):
    d = degs[degs["contrast"] == "Lean vs Obese T2D"].copy()
    d["score"] = np.sign(d["log2FC"]) * -np.log10(d["pvalue"].clip(lower=1e-300))
    rnk = d[["gene", "score"]].dropna().sort_values("score", ascending=False)
    rnk = rnk.groupby("gene").mean().sort_values("score", ascending=False)
    return rnk


def _fallback_ora(sig_genes, universe):
    from scipy.stats import hypergeom
    N = len(set(universe))
    n = len(set(sig_genes))
    rows = []
    sg = set(sig_genes)
    for name, members in _FALLBACK_SETS.items():
        members = [g for g in members if g in universe]
        K = len(members)
        if K == 0:
            continue
        k = len(sg & set(members))
        p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
        rows.append({"Term": name, "Overlap": f"{k}/{K}",
                     "P-value": p, "Genes": ";".join(sorted(sg & set(members))),
                     "source": "fallback_hypergeom"})
    df = pd.DataFrame(rows).sort_values("P-value")
    return df


def run(degs, expr, cache=True):
    rnk = _preranked(degs)
    universe = expr.index.tolist()
    sig = degs[(degs["contrast"] == "Lean vs Obese T2D") & degs["significant"]]["gene"].tolist()
    if len(sig) < 5:  # relax if strict threshold yields too few
        sig = (degs[degs["contrast"] == "Lean vs Obese T2D"]
               .sort_values("pvalue")["gene"].head(150).tolist())

    gsea_df = None
    used_online = False
    try:
        import gseapy as gp
        pre = gp.prerank(rnk=rnk.reset_index(), gene_sets=GENE_SETS,
                         min_size=5, max_size=1000, permutation_num=100,
                         outdir=os.path.join(C.RES_DIR, "gseapy_prerank"),
                         seed=C.RANDOM_STATE, no_plot=True, verbose=False)
        gsea_df = pre.res2d.copy()
        # normalise column names
        gsea_df = gsea_df.rename(columns={"Term": "Term", "NES": "NES",
                                          "FDR q-val": "FDR", "NOM p-val": "pvalue"})
        used_online = True
        print(f"  [pathway] GSEA prerank complete: {len(gsea_df)} gene sets", flush=True)
    except Exception as e:
        print(f"  [pathway] online GSEA unavailable ({str(e)[:80]}); "
              f"using offline hypergeometric ORA", flush=True)

    if gsea_df is None or len(gsea_df) == 0:
        gsea_df = _fallback_ora(sig, universe)
        gsea_df["NES"] = -np.log10(gsea_df["P-value"].clip(lower=1e-300))
        gsea_df["FDR"] = gsea_df["P-value"]

    # ORA via Enrichr (best-effort)
    ora_df = None
    if used_online:
        try:
            import gseapy as gp
            enr = gp.enrichr(gene_list=sig, gene_sets=GENE_SETS,
                             outdir=os.path.join(C.RES_DIR, "gseapy_enrichr"),
                             no_plot=True)
            ora_df = enr.results
        except Exception as e:
            print(f"  [pathway] Enrichr ORA unavailable ({str(e)[:60]})", flush=True)
    if ora_df is None:
        ora_df = _fallback_ora(sig, universe)

    if cache:
        gsea_df.to_csv(os.path.join(C.TAB_DIR, "table_gsea.csv"), index=False)
        ora_df.to_csv(os.path.join(C.TAB_DIR, "table_ora.csv"), index=False)
    return {"gsea": gsea_df, "ora": ora_df, "online": used_online, "ranked": rnk}
