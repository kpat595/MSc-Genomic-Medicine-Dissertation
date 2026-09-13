#!/usr/bin/env python3

import scanpy as sc
import cellhint
import pandas as pd
from pathlib import Path

CONDITIONS = ["NR1", "CASPR2", "Control", "NoTreated"]
ROOT = Path("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/cellhint")
TABLE_ROOT = Path("/well/bsg/projects/CGG0M1_Patel/seurat_analysis/.by_condition/tables_by_condition")
OUT_DIR = ROOT / "harmonisation_output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def best_resolution(cond):
    cmp = pd.read_csv(TABLE_ROOT / cond / "05_bootstrap_stability_resolution_comparison.tsv", sep="\t")
    return cmp.loc[cmp["pct_cells_stable_median"].idxmax(), "resolution"]

adatas = {}
for cond in CONDITIONS:
    res = str(best_resolution(cond)).replace(".", "_")
    h5ad_path = ROOT / "h5ads" / f"01_{cond}.h5ad"
    adata = sc.read_h5ad(h5ad_path)
    adata.obs["condition"] = cond
    adata.obs["cellhint_cluster"] = adata.obs[f"cluster_sct_res_{res}"]
    adatas[cond] = adata
    print(f"{cond}: using resolution {res}")

adata = sc.concat(adatas, join="outer", index_unique="-")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="condition")
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=30, use_highly_variable=True)

alignment = cellhint.harmonize(
    adata,
    dataset="condition",
    cell_type="cellhint_cluster",
)

alignment.relation.to_csv(OUT_DIR / "02_harmonisation_relation.csv")
pd.DataFrame(alignment.groups).to_csv(OUT_DIR / "02_harmonisation_groups.csv")  # groups is a numpy array, not a DataFrame
alignment.reannotation.to_csv(OUT_DIR / "02_harmonisation_reannotation.csv")

cellhint.treeplot(alignment, save=str(OUT_DIR / "02_treeplot.png"))
cellhint.sankeyplot(alignment, save=str(OUT_DIR / "02_sankeyplot.html"), show=False)  # show=True (the default) hangs on a headless node

adata.write_h5ad(OUT_DIR / "02_harmonised.h5ad")

print("Done. Outputs in", OUT_DIR)
