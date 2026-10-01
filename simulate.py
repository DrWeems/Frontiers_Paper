"""
simulate.py -- Graceful fallback data generator.

If GEO download fails, this reproduces (and improves on) the simulation logic
from the original hepato_switch_analysis.py so the full pipeline still runs
end-to-end.  All outputs produced from this path are flagged SIMULATED.

Produces the same structures the real ingest/preprocess path yields:
    discovery_expr (genes x 103), discovery_meta, validation_expr (genes x 21),
    validation_meta -- so every downstream module is agnostic to data source.
"""
import numpy as np
import pandas as pd

import config as C


def _make_pool(n_per, seed, n_bg=500):
    rng = np.random.default_rng(seed)
    groups = (["Lean"] * n_per[0] + ["Obese ND"] * n_per[1] + ["Obese T2D"] * n_per[2])
    n = len(groups)
    axis = np.array([C.GROUP_ORDER[g] for g in groups], dtype=float)
    axis_c = (axis - axis.mean()) / axis.std()

    genes = list(C.SWITCH_GENES) + [f"BG_{i}" for i in range(n_bg)]
    mat = np.zeros((len(genes), n))
    base = rng.normal(7, 1.5, len(genes))
    for gi, g in enumerate(genes):
        if g in C.SWITCH_GENES:
            eff = C.SWITCH_DIRECTION.get(g, 1) * rng.uniform(0.9, 1.4)
        else:
            eff = rng.uniform(-0.5, 0.5)
        mat[gi] = base[gi] + eff * axis_c + rng.normal(0, 0.45, n)
    samples = [f"S{seed}_{i}" for i in range(n)]
    expr = pd.DataFrame(mat, index=genes, columns=samples)
    meta = pd.DataFrame({"sample": samples, "group": groups,
                         "dataset": f"SIM{seed}", "platform": "simulated",
                         "reconstructed": True}).set_index("sample")
    return expr, meta


def run(cache=True):
    print("  [simulate] generating simulated discovery (n=103) + validation (n=21)",
          flush=True)
    disc_expr, disc_meta = _make_pool((23, 35, 45), seed=42)
    val_expr, val_meta = _make_pool((6, 8, 7), seed=7)
    val_expr = val_expr.loc[disc_expr.index]
    return {"discovery_expr": disc_expr, "discovery_meta": disc_meta,
            "validation_expr": val_expr, "validation_meta": val_meta}
