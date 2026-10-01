"""
ml_classifier.py -- Random Forest staging with LOOCV + independent validation.

* Feature set : the 10-gene molecular switch signature (falls back to top DEGs
  if some switch genes are missing).
* Task        : binary "Obese T2D vs others" for the ROC analysis, plus a
  3-class staging report.
* Discovery   : Leave-One-Out cross-validation over the harmonized n=103 pool.
* Validation  : model retrained on the full discovery pool, then applied to the
  independent GSE64998 cohort (n=21). Directional concordance of the signature
  is reported.
"""
import os
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc, accuracy_score

import config as C


def _feature_genes(expr, degs=None):
    genes = [g for g in C.SWITCH_GENES if g in expr.index]
    if len(genes) < 10 and degs is not None:
        extra = (degs[degs["contrast"] == "Lean vs Obese T2D"]
                 .sort_values("padj")["gene"].tolist())
        for g in extra:
            if g in expr.index and g not in genes:
                genes.append(g)
            if len(genes) >= 10:
                break
    return genes


def run(disc_expr, disc_meta, val_expr, val_meta, degs=None, cache=True):
    genes = _feature_genes(disc_expr, degs)
    print(f"  [ML] signature genes ({len(genes)}): {genes}", flush=True)

    # Per-gene z-scoring WITHIN each cohort removes platform/cohort-level
    # location & scale differences so the switch signature is comparable across
    # the ComBat-corrected discovery pool and the independent validation cohort.
    def _zscore(df):
        return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1) + 1e-9, axis=0)

    Xd_df = _zscore(disc_expr.loc[genes, disc_meta.index])
    Xd = Xd_df.T.values
    yd = (disc_meta["group"] == "Obese T2D").astype(int).values

    # ---- Discovery LOOCV ----
    loo = LeaveOneOut()
    probs, truth = [], []
    for tr, te in loo.split(Xd):
        clf = RandomForestClassifier(n_estimators=300, max_depth=5,
                                     random_state=C.RANDOM_STATE, n_jobs=-1)
        clf.fit(Xd[tr], yd[tr])
        probs.append(clf.predict_proba(Xd[te])[0, 1])
        truth.append(yd[te][0])
    fpr_d, tpr_d, _ = roc_curve(truth, probs)
    auc_d = auc(fpr_d, tpr_d)
    acc_d = accuracy_score(truth, np.array(probs) > 0.5)
    print(f"  [ML] discovery LOOCV AUC = {auc_d:.3f}  (acc={acc_d:.2f})", flush=True)

    # ---- Validation ----
    final = RandomForestClassifier(n_estimators=300, max_depth=5,
                                   random_state=C.RANDOM_STATE, n_jobs=-1)
    final.fit(Xd, yd)
    val_genes = [g for g in genes if g in val_expr.index]
    Xv_df = _zscore(val_expr.loc[val_genes, val_meta.index])   # z-score within validation cohort
    Xv = Xv_df.T.reindex(columns=genes)
    Xv = Xv.fillna(0.0)   # z-scored space -> missing gene = cohort mean (0)
    yv = (val_meta["group"] == "Obese T2D").astype(int).values
    pv = final.predict_proba(Xv.values)[:, 1]
    if len(np.unique(yv)) > 1:
        fpr_v, tpr_v, _ = roc_curve(yv, pv)
        auc_v = auc(fpr_v, tpr_v)
    else:
        fpr_v, tpr_v, auc_v = np.array([0, 1]), np.array([0, 1]), float("nan")
    print(f"  [ML] validation AUROC = {auc_v:.3f}", flush=True)

    importances = pd.Series(final.feature_importances_, index=genes).sort_values(ascending=False)

    result = {
        "genes": genes,
        "discovery": {"fpr": fpr_d.tolist(), "tpr": tpr_d.tolist(),
                      "auc": float(auc_d), "acc": float(acc_d), "n": int(len(yd))},
        "validation": {"fpr": fpr_v.tolist(), "tpr": tpr_v.tolist(),
                       "auc": float(auc_v), "n": int(len(yv))},
        "importances": importances.to_dict(),
    }
    if cache:
        with open(os.path.join(C.RES_DIR, "ml_results.json"), "w") as fh:
            json.dump(result, fh, indent=2)
        importances.to_csv(os.path.join(C.TAB_DIR, "table_feature_importance.csv"))
    return result
