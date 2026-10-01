#!/usr/bin/env python3
"""
pipeline.py -- HepatoSwitch Hepatic Diabesity Transcriptomics Pipeline

Executable orchestrator for the HepatoSwitch modular workflow supporting
the Frontiers in Endocrinology manuscript:
"An Integrative Multi-Platform Transcriptomic Framework for Molecular
Staging and Identification of Regulatory Drivers of the Hepatic Molecular
Switch in Diabesity."

Place this file in the same directory as:
    config.py, geo_ingest.py, preprocess.py, differential_expression.py,
    pathway_enrichment.py, network_switches.py, ml_classifier.py,
    figures.py, simulate.py

Usage
-----
    python pipeline.py            # full run (downloads GEO data)
    python pipeline.py --resume   # reuse cached harmonized matrices
    python pipeline.py --simulate # force simulated-data fallback

Stages
------
    1. GEO ingestion            (GSE121344, GSE15653, GSE48452, GSE64998)
    2. Preprocessing            (normalization, harmonization, ComBat)
    3. Differential expression  (3 ordered contrasts, BH-FDR)
    4. Pathway enrichment       (GSEA + ORA, MSigDB Hallmark + KEGG)
    5. Molecular-switch prioritization  (Spearman/Pearson, Table II)
    6. Random Forest staging    (LOOCV, discovery) + independent validation
    7. Publication figures (6, 300 DPI) + results tables + run_summary.json

If live GEO retrieval or preprocessing fails at any point the pipeline
switches to a clearly-flagged SIMULATED data path so it always produces a
complete set of outputs. Every figure title is annotated [SIMULATED] in
that mode.
"""

import argparse
import os
import traceback
import warnings

import pandas as pd

import config as C
import differential_expression
import figures
import geo_ingest
import ml_classifier
import network_switches
import pathway_enrichment
import preprocess
import simulate

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(text):
    print("\n" + "=" * 74)
    print(f"  {text}")
    print("=" * 74, flush=True)


def load_cached_data():
    """Load cached harmonized matrices and report whether they are simulated."""
    pre = {
        "discovery_expr": pd.read_parquet(
            os.path.join(C.RES_DIR, "discovery_expr.parquet")
        ),
        "discovery_meta": pd.read_csv(
            os.path.join(C.RES_DIR, "discovery_meta.csv"), index_col=0
        ),
        "validation_expr": pd.read_parquet(
            os.path.join(C.RES_DIR, "validation_expr.parquet")
        ),
        "validation_meta": pd.read_csv(
            os.path.join(C.RES_DIR, "validation_meta.csv"), index_col=0
        ),
    }
    simulated = (
        "platform" in pre["discovery_meta"].columns
        and pre["discovery_meta"]["platform"].eq("simulated").any()
    )
    return pre, simulated


def get_data(args):
    """Return the preprocessed data dictionary and simulation-status flag."""
    cache_files = [
        "discovery_expr.parquet",
        "discovery_meta.csv",
        "validation_expr.parquet",
        "validation_meta.csv",
    ]
    cache_ready = all(
        os.path.exists(os.path.join(C.RES_DIR, f)) for f in cache_files
    )

    if args.simulate:
        banner("SIMULATION MODE (forced)")
        return simulate.run(), True

    if args.resume and cache_ready:
        banner("STAGES 1-2  Resuming from cached harmonized matrices")
        return load_cached_data()

    try:
        banner("STAGE 1  GEO data ingestion")
        datasets = geo_ingest.ingest_all()

        banner("STAGE 2  Preprocessing, harmonization, and batch correction")
        return preprocess.run(datasets), False

    except Exception as error:
        print("\n[WARNING] Live GEO ingestion/preprocessing failed.")
        print(f"          {type(error).__name__}: {str(error)[:300]}")
        print("          Switching to clearly-labelled simulation fallback.")
        traceback.print_exc(limit=2)
        return simulate.run(), True


def validate_alignment(expr, meta, cohort_name):
    """Ensure expression columns and metadata indices describe identical samples."""
    missing = set(expr.columns) - set(meta.index)
    if missing:
        raise ValueError(
            f"{cohort_name} metadata missing samples: {sorted(missing)[:10]}"
        )
    meta = meta.loc[expr.columns].copy()
    if not expr.columns.equals(meta.index):
        raise ValueError(
            f"{cohort_name} expression and metadata are not aligned."
        )
    if "group" not in meta.columns:
        raise ValueError(
            f"{cohort_name} metadata do not contain a 'group' column."
        )
    return meta


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="HepatoSwitch hepatic diabesity transcriptomics pipeline"
    )
    parser.add_argument(
        "--simulate", action="store_true", help="force simulated data"
    )
    parser.add_argument(
        "--resume", action="store_true", help="reuse cached harmonized matrices"
    )
    # parse_known_args avoids crashes when run inside Jupyter (extra kernel args)
    args, _ = parser.parse_known_args()

    banner("HEPATIC DIABESITY TRANSCRIPTOMICS PIPELINE")
    print(f"  Output directory: {C.OUT_DIR}")

    # ------------------------------------------------------------------
    # STAGES 1-2: data ingestion and preprocessing
    # ------------------------------------------------------------------
    pre, simulated = get_data(args)
    data_mode = "SIMULATED" if simulated else "REAL GEO DATA"
    print(f"\n  >>> DATA MODE: {data_mode} <<<")

    disc_expr = pre["discovery_expr"]
    disc_meta = validate_alignment(disc_expr, pre["discovery_meta"], "Discovery cohort")
    val_expr  = pre["validation_expr"]
    val_meta  = validate_alignment(val_expr,  pre["validation_meta"], "Validation cohort")

    # Write validated metadata back so downstream stages use aligned copies
    pre["discovery_meta"]  = disc_meta
    pre["validation_meta"] = val_meta

    print(
        f"  Discovery pool : {disc_expr.shape[0]} genes x "
        f"{disc_expr.shape[1]} samples"
    )
    print(f"  Validation cohort: {val_expr.shape[1]} samples")
    print("\n  Discovery cohort composition:")
    print(disc_meta["group"].value_counts().reindex(C.GROUPS, fill_value=0))
    print("\n  Validation cohort composition:")
    print(val_meta["group"].value_counts().reindex(C.GROUPS, fill_value=0))

    # ------------------------------------------------------------------
    # STAGE 3: differential expression
    # ------------------------------------------------------------------
    banner("STAGE 3  Differential expression (3 contrasts)")
    degs = differential_expression.run(disc_expr, disc_meta)
    target = degs[degs["contrast"] == "Lean vs Obese T2D"]
    print(
        f"  Lean vs Obese T2D: {len(target)} genes tested, "
        f"{int(target['significant'].sum())} significant (adj. p < {C.DEG_ALPHA}, "
        f"|log2FC| > {C.DEG_LOG2FC})"
    )

    # ------------------------------------------------------------------
    # STAGE 4: pathway enrichment
    # ------------------------------------------------------------------
    banner("STAGE 4  Pathway enrichment (GSEA + ORA)")
    enrichment = pathway_enrichment.run(degs, disc_expr)

    # ------------------------------------------------------------------
    # STAGE 5: molecular-switch prioritization
    # ------------------------------------------------------------------
    banner("STAGE 5  Molecular-switch prioritization")
    switches = network_switches.run(disc_expr, disc_meta)
    display_cols = [
        col for col in
        ["gene", "present", "spearman_r", "spearman_p", "FDR", "concordant"]
        if col in switches.columns
    ]
    print(switches[display_cols].to_string(index=False))

    # ------------------------------------------------------------------
    # STAGE 6: machine-learning staging + independent validation
    # ------------------------------------------------------------------
    banner("STAGE 6  Molecular-staging classifier (LOOCV + validation)")
    ml_results = ml_classifier.run(
        disc_expr, disc_meta, val_expr, val_meta, degs
    )

    # ------------------------------------------------------------------
    # STAGE 7: publication figures and tables
    # ------------------------------------------------------------------
    banner("STAGE 7  Publication figures and tables")
    figure_paths = figures.run_all(
        pre, degs, enrichment["gsea"], ml_results, simulated
    )

    # ------------------------------------------------------------------
    # Write machine-readable run summary
    # ------------------------------------------------------------------
    summary = {
        "data_mode":        data_mode,
        "discovery_auc":    ml_results["discovery"]["auc"],
        "validation_auroc": ml_results["validation"]["auc"],
        "n_discovery":      int(disc_meta.shape[0]),
        "n_validation":     int(val_meta.shape[0]),
        "gsea_online":      enrichment.get("online", False),
        "figures":          [os.path.basename(f) for f in figure_paths],
    }
    summary_path = os.path.join(C.RES_DIR, "run_summary.json")
    pd.Series(summary).to_json(summary_path, indent=2)

    # ------------------------------------------------------------------
    # Final report
    # ------------------------------------------------------------------
    banner("PIPELINE COMPLETE")
    print(f"  Data mode              : {data_mode}")
    print(f"  Discovery AUC (LOOCV)  : {ml_results['discovery']['auc']:.3f}")
    print(f"  Validation AUROC       : {ml_results['validation']['auc']:.3f}")
    print(f"  Figures ({len(figure_paths)})           -> {C.FIG_DIR}")
    print(f"  Tables                 -> {C.TAB_DIR}")
    for f in sorted(os.listdir(C.TAB_DIR)):
        print(f"      - {f}")
    print(f"  Run summary            -> {summary_path}")
    print("\n  Done.")


if __name__ == "__main__":
    main()
