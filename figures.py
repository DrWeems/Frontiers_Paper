"""
figures.py -- All publication-ready figures (300 DPI, manuscript colour scheme).

Figures
-------
fig_consort.png : dataset harmonization / flow diagram
fig_pca.png     : PCA of harmonized n=103 pool coloured by group
fig_volcano.png : Lean vs Obese T2D volcano, switch genes highlighted
fig_heatmap.png : top-50 DEG heatmap across all discovery samples
fig_gsea.png    : top enriched pathways (NES bar / enrichment)
fig_roc.png     : discovery + validation ROC curves
"""
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")
import config as C

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = C.FONT_FAMILY
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["figure.dpi"] = 110


def _save(fig, name):
    path = os.path.join(C.FIG_DIR, name)
    fig.savefig(path, dpi=C.DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"      saved {name}", flush=True)
    return path


# ----------------------------------------------------------------------------
def fig_consort(meta, val_meta, simulated=False):
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.axis("off")
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)

    def box(x, y, w, h, text, color="#ecf0f1", ec="#34495e"):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle="round,pad=0.08", fc=color, ec=ec, lw=1.4))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=9, wrap=True)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                     arrowstyle="-|>", mutation_scale=14, lw=1.3, color="#34495e"))

    # top: four datasets
    box(0.2, 8.3, 2.1, 1.2, "GSE121344\nRNA-seq\nn=12", "#fdf2e9")
    box(2.6, 8.3, 2.1, 1.2, "GSE15653\nMicroarray\nn=18", "#eaf2f8")
    box(5.0, 8.3, 2.1, 1.2, "GSE48452\nMicroarray\nn=73", "#eaf2f8")
    box(7.5, 8.3, 2.3, 1.2, "GSE64998\nMicroarray\nn=21\n(validation)", "#f5eef8")

    # normalization
    box(1.3, 6.3, 4.6, 1.0, "Per-platform normalization\n(log2 / RMA-like harmonization)", "#eafaf1")
    for x in (1.25, 3.65, 6.05):
        arrow(x, 8.3, 3.6, 7.35)

    # combat
    box(1.3, 4.4, 4.6, 1.0, "Common genes + ComBat batch correction", "#eafaf1")
    arrow(3.6, 6.3, 3.6, 5.45)

    # discovery pool
    nl = int((meta['group'] == 'Lean').sum())
    nn = int((meta['group'] == 'Obese ND').sum())
    nt = int((meta['group'] == 'Obese T2D').sum())
    box(1.0, 2.3, 5.2, 1.4,
        f"Harmonized discovery pool  n={len(meta)}\n"
        f"Lean={nl}   Obese ND={nn}   Obese T2D={nt}", "#d5f5e3")
    arrow(3.6, 4.4, 3.6, 3.75)

    # analyses
    box(0.6, 0.3, 5.9, 1.4,
        "DEG - GSEA/ORA - Switch prioritization\nRandom Forest staging (LOOCV)", "#fcf3cf")
    arrow(3.5, 2.3, 3.5, 1.75)

    # validation arm
    vl = int((val_meta['group'] == 'Lean').sum())
    vn = int((val_meta['group'] == 'Obese ND').sum())
    vt = int((val_meta['group'] == 'Obese T2D').sum())
    box(7.3, 2.3, 2.6, 1.4,
        f"Independent validation\nn={len(val_meta)}\nLean={vl} ND={vn} T2D={vt}", "#ebdef0")
    arrow(8.6, 8.3, 8.6, 3.75)
    box(7.3, 0.3, 2.6, 1.4, "External AUROC\n+ concordance", "#fcf3cf")
    arrow(8.6, 2.3, 8.6, 1.75)
    arrow(6.5, 1.0, 7.3, 1.0)

    title = "Figure 1. Multi-platform dataset harmonization and analysis workflow"
    if simulated:
        title += "  [SIMULATED DATA]"
    ax.set_title(title, fontsize=11, fontweight="bold")
    return _save(fig, "fig_consort.png")


# ----------------------------------------------------------------------------
def fig_pca(expr, meta, simulated=False):
    X = expr[meta.index].T.values
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    pca = PCA(n_components=2, random_state=C.RANDOM_STATE)
    pcs = pca.fit_transform(X)
    ev = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(7, 6))
    for grp in C.GROUPS:
        idx = (meta["group"] == grp).values
        ax.scatter(pcs[idx, 0], pcs[idx, 1], c=C.COLORS[grp], label=grp,
                   s=70, edgecolors="white", linewidths=0.8, alpha=0.9)
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}% variance) - progression axis", fontsize=11)
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}% variance)", fontsize=11)
    t = f"Figure 2. PCA of harmonized discovery pool (n={len(meta)})"
    if simulated:
        t += "  [SIMULATED]"
    ax.set_title(t, fontsize=11, fontweight="bold")
    ax.legend(title="Group", frameon=False)
    ax.grid(True, linestyle="--", alpha=0.35)
    return _save(fig, "fig_pca.png")


# ----------------------------------------------------------------------------
def fig_volcano(degs, simulated=False):
    d = degs[degs["contrast"] == "Lean vs Obese T2D"].copy()
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    sig = d["significant"]
    ax.scatter(d.loc[~sig, "log2FC"], d.loc[~sig, "nlp"], s=10,
               c="#bdc3c7", alpha=0.5, label="NS")
    up = sig & (d["log2FC"] > 0)
    dn = sig & (d["log2FC"] < 0)
    ax.scatter(d.loc[up, "log2FC"], d.loc[up, "nlp"], s=14, c="#c0392b", alpha=0.7, label="Up in T2D")
    ax.scatter(d.loc[dn, "log2FC"], d.loc[dn, "nlp"], s=14, c="#2980b9", alpha=0.7, label="Down in T2D")

    # highlight switch genes
    sw = d[d["gene"].isin(C.SWITCH_GENES)]
    ax.scatter(sw["log2FC"], sw["nlp"], s=90, facecolors="none",
               edgecolors="#f39c12", linewidths=1.8, label="Switch gene", zorder=5)
    try:
        from adjustText import adjust_text
        texts = [ax.text(r["log2FC"], r["nlp"], r["gene"], fontsize=8, fontweight="bold")
                 for _, r in sw.iterrows()]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#7f8c8d", lw=0.5))
    except Exception:
        for _, r in sw.iterrows():
            ax.annotate(r["gene"], (r["log2FC"], r["nlp"]), fontsize=8, fontweight="bold")

    ax.axhline(-np.log10(C.DEG_ALPHA), ls="--", c="#7f8c8d", lw=0.8)
    ax.axvline(C.DEG_LOG2FC, ls="--", c="#7f8c8d", lw=0.8)
    ax.axvline(-C.DEG_LOG2FC, ls="--", c="#7f8c8d", lw=0.8)
    ax.set_xlabel("log2 fold-change (Obese T2D vs Lean)", fontsize=11)
    ax.set_ylabel("-log10 adjusted p-value", fontsize=11)
    t = "Figure 3. Differential expression: Lean vs Obese T2D"
    if simulated:
        t += "  [SIMULATED]"
    ax.set_title(t, fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    return _save(fig, "fig_volcano.png")


# ----------------------------------------------------------------------------
def fig_heatmap(expr, meta, degs, simulated=False):
    top = (degs[degs["contrast"] == "Lean vs Obese T2D"]
           .sort_values("padj").head(50)["gene"].tolist())
    top = [g for g in top if g in expr.index]
    order = meta.sort_values("group").index
    M = expr.loc[top, order]
    Z = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1) + 1e-9, axis=0)
    Z = Z.clip(-3, 3)

    fig, (axc, axh) = plt.subplots(2, 1, figsize=(11, 9),
                                   gridspec_kw={"height_ratios": [0.5, 12], "hspace": 0.02})
    # group colour bar
    grp = meta.loc[order, "group"]
    cbar = np.array([[list(int(C.COLORS[g].lstrip("#")[i:i+2], 16) / 255 for i in (0, 2, 4))
                      for g in grp]])
    axc.imshow(cbar, aspect="auto")
    axc.set_xticks([]); axc.set_yticks([])
    axc.set_title("Figure 4. Top-50 DEG heatmap (Lean vs Obese T2D) across discovery samples"
                  + ("  [SIMULATED]" if simulated else ""),
                  fontsize=11, fontweight="bold")

    im = axh.imshow(Z.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    axh.set_yticks(range(len(top)))
    axh.set_yticklabels(top, fontsize=6)
    axh.set_xticks([])
    axh.set_xlabel(f"Samples (n={len(order)}), grouped Lean -> Obese ND -> Obese T2D", fontsize=10)
    cb = fig.colorbar(im, ax=axh, fraction=0.02, pad=0.01)
    cb.set_label("row z-score", fontsize=9)

    handles = [plt.Line2D([0], [0], marker="s", ls="", markersize=9,
               markerfacecolor=C.COLORS[g], label=g) for g in C.GROUPS]
    axh.legend(handles=handles, title="Group", bbox_to_anchor=(1.13, 1),
               loc="upper left", frameon=False, fontsize=8)
    return _save(fig, "fig_heatmap.png")


# ----------------------------------------------------------------------------
def fig_gsea(gsea, simulated=False):
    df = gsea.copy()
    # normalise NES / term columns
    nes_col = "NES" if "NES" in df.columns else df.columns[-1]
    df[nes_col] = pd.to_numeric(df[nes_col], errors="coerce")
    df = df.dropna(subset=[nes_col])
    df["absnes"] = df[nes_col].abs()
    df = df.sort_values("absnes", ascending=False).head(10)
    df = df.sort_values(nes_col)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    colors = ["#c0392b" if v > 0 else "#2980b9" for v in df[nes_col]]
    terms = [str(t)[:45] for t in df["Term"]]
    ax.barh(terms, df[nes_col], color=colors, edgecolor="white")
    ax.axvline(0, color="#34495e", lw=0.8)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
    t = "Figure 5. Top enriched pathways (Lean vs Obese T2D)"
    if simulated:
        t += "  [SIMULATED]"
    ax.set_title(t, fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8)
    return _save(fig, "fig_gsea.png")


# ----------------------------------------------------------------------------
def fig_roc(ml, simulated=False):
    d, v = ml["discovery"], ml["validation"]
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot(d["fpr"], d["tpr"], color=C.ACCENT["discovery"], lw=2.4,
            label=f"Discovery (n={d['n']}, AUC={d['auc']:.2f})")
    if not np.isnan(v["auc"]):
        ax.plot(v["fpr"], v["tpr"], color=C.ACCENT["validation"], lw=2.4,
                label=f"Validation (n={v['n']}, AUC={v['auc']:.2f})")
    ax.plot([0, 1], [0, 1], ls="--", color=C.ACCENT["diag"], lw=1)
    ax.set_xlabel("False positive rate", fontsize=11)
    ax.set_ylabel("True positive rate", fontsize=11)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
    t = "Figure 6. Random Forest staging performance (Obese T2D vs others)"
    if simulated:
        t += "  [SIMULATED]"
    ax.set_title(t, fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(True, alpha=0.3)
    return _save(fig, "fig_roc.png")


def run_all(pre, degs, gsea, ml, simulated=False):
    print("  [figures] generating publication figures ...", flush=True)
    paths = []
    paths.append(fig_consort(pre["discovery_meta"], pre["validation_meta"], simulated))
    paths.append(fig_pca(pre["discovery_expr"], pre["discovery_meta"], simulated))
    paths.append(fig_volcano(degs, simulated))
    paths.append(fig_heatmap(pre["discovery_expr"], pre["discovery_meta"], degs, simulated))
    paths.append(fig_gsea(gsea, simulated))
    paths.append(fig_roc(ml, simulated))
    return paths
