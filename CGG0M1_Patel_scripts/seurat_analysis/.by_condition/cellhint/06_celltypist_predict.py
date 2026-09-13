#!/usr/bin/env python3
# Separate HNOCA/Braun CellTypist models (avoids reconciling their differing label
# taxonomies). Idempotent: skips training/prediction if the .pkl/CSV already exists.

import celltypist
import scanpy as sc
import pandas as pd
from pathlib import Path

CONDITIONS = ["NR1", "CASPR2", "Control", "NoTreated"]
ROOT = Path("/well/bsg/projects/CGG0M1_Patel/seurat_analysis")
REF_DIR = ROOT / "reference"
H5AD_DIR = ROOT / ".by_condition" / "cellhint" / "h5ads"
OUT_DIR = ROOT / ".by_condition" / "cellhint" / "celltypist_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MODELS = {
    "hnoca": {
        "ref_h5ad": REF_DIR / "ref_hnoca_subsampled.h5ad",
        "label_col": "annot_level_2",
        "model_path": OUT_DIR / "model_hnoca.pkl",
    },
    "braun": {
        "ref_h5ad": REF_DIR / "ref_braun_subsampled.h5ad",
        "label_col": "CellClass",
        "model_path": OUT_DIR / "model_braun.pkl",
    },
}


def normalize(adata):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def train_model(name, cfg):
    if cfg["model_path"].exists():
        print(f"{name}: model already exists at {cfg['model_path']}, skipping training")
        return celltypist.models.Model.load(str(cfg["model_path"]))

    print(f"{name}: training on {cfg['ref_h5ad']}")
    ref = sc.read_h5ad(cfg["ref_h5ad"])
    ref = normalize(ref)
    model = celltypist.train(ref, labels=cfg["label_col"], n_jobs=-1, feature_selection=True)
    model.write(cfg["model_path"])
    print(f"{name}: model saved to {cfg['model_path']}")
    return model


def predict_condition(cond, models):
    out_csv = OUT_DIR / f"06_{cond}_celltypist.csv"
    if out_csv.exists():
        print(f"{cond}: {out_csv} already exists, skipping prediction")
        return pd.read_csv(out_csv)

    print(f"{cond}: predicting")
    adata = sc.read_h5ad(H5AD_DIR / f"01_{cond}.h5ad")
    adata = normalize(adata)

    result = pd.DataFrame({"barcode": adata.obs_names})
    for name, model in models.items():
        pred = celltypist.annotate(adata, model=model)
        result[f"celltypist_{name}_label"] = pred.predicted_labels["predicted_labels"].values
        result[f"celltypist_{name}_conf"] = pred.probability_matrix.max(axis=1).values

    result.to_csv(out_csv, index=False)
    print(f"{cond}: wrote {out_csv}")
    return result


def main():
    models = {name: train_model(name, cfg) for name, cfg in MODELS.items()}

    all_results = []
    for cond in CONDITIONS:
        res = predict_condition(cond, models)
        all_results.append(res)

    combined = pd.concat(all_results, ignore_index=True)
    combined_path = OUT_DIR / "06_celltypist_predictions.csv"
    combined.to_csv(combined_path, index=False)
    print("Done. Combined output:", combined_path)


if __name__ == "__main__":
    main()
