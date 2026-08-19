#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File:    Elvin_Tcf15_Manuscript
Author:  David Rodriguez Morales

"""
import os
import pathlib

from tqdm import tqdm
from prelude_py import ad, plt, sc, do, dc
from dotools_py.pl._plot_utils import spine_format

# ! Figure 1e,f
# region Expression of TCF5 in AMI Spatial Transcriptomics
input_path = pathlib.Path("/mnt/davidr/scStorage/DavidR/BioData/AMI_Kuppe_Spatial_2022/Objects")
out_path = pathlib.Path("/mnt/davidr/ZMM_shared/DavidR/Requests/Elvin/Tcf15_Spatial")
os.makedirs(out_path, exist_ok=True)

files = [f for f in os.listdir(input_path) if "snRNA" not in f]
skip_regions = ["BZ", "IZ_BZ", "RZ_BZ", "RZ_FZ"]

# * Figure 1e
for f in files:
    adata = ad.read_h5ad(out_path / f)
    name = f.split(".h5ad")[0].strip()
    adata.var_names = adata.var["feature_name"].values
    del adata.raw
    library = list(adata.uns["spatial"].keys())
    library.remove("is_single")
    if len(library) != 1:
        raise ValueError
    fig, axs = plt.subplots(1, 1, figsize=(15, 6))
    ax = sc.pl.spatial(
        adata, color="TCF15", library_id=library[0], size=1.5, ax=axs, show=False, cmap="magma", vmax=2)
    spine_format(ax[0], "SP")
    plt.savefig(out_path / f"SP_Visium_{name}_TCF15.pdf", bbox_inches="tight")


database = {}
for f in tqdm(files):
    if any(pattern in f for pattern in skip_regions):
        continue

    adata = do.io.read_h5ad(input_path / f)
    adata.var_names = adata.var["feature_name"].values
    adata.var_names_make_unique()

    try:
        adata.layers["counts"] = adata.raw.X.copy()
        del adata.raw
    except AttributeError as e:
        adata.layers["counts"] = adata.X.copy()
    name = f.split(".h5ad")[0].strip()
    database[name] = adata

adata_concat = ad.concat(database.values(), keys=database.keys(), label="name", fill_value=0, index_unique="-")
adata_concat.obs["region"] = adata_concat.obs["name"].map({
    'IZ_BZ_P2': "IZ_BZ", 'RZ_BZ_P2': "RZ_BZ", 'IZ_P3': "IZ", 'RZ_P6': "RZ",
    'RZ_GT_P2': "RZ", 'control_P1': "Control", 'IZ_P15': "IZ",
    'control_P7': "Control", 'RZ_P6_2': "RZ", 'FZ_GT_P4': "FZ",
    'FZ_P20': "FZ", 'control_P17': "Control", 'control_P8': "Control",
    'GT_IZ_P13': "IZ", 'RZ_BZ_P3': "RZ_BZ", 'RZ_FZ_P5': "RZ_FZ",
    'RZ_P9': "RZ", 'IZ_P10': "IZ", 'RZ_P11': "RZ", 'RZ_P3': "RZ",
    'FZ_P14': "FZ", 'GT_IZ_P9': "IZ", 'GT_IZ_P9_rep2': "IZ",
    'IZ_P16': "IZ", 'RZ_BZ_P12': "RZ_BZ", 'FZ_P18': "FZ",
    'FZ_GT_P19': "FZ", 'GT_IZ_P15': "IZ"
})
adata_concat.X = adata_concat.layers["counts"].copy()
do.pp.log_normalize(adata_concat, target_sum=10_000)
adata_concat = adata_concat[adata_concat.obs.region.isin(["Control", "RZ", "IZ", "FZ"])].copy()

do.utility.free_memory()

sc.pp.filter_genes(adata_concat, min_cells=5)
pdata = dc.pp.pseudobulk(adata_concat, sample_col="donor_id", groups_col="region", layer="counts")
pdata.layers["counts"] = pdata.X.copy()
dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000)

tester = do.tl.DGEAnalysis(adata=pdata, groupby="region", batch_key="donor_id", is_pseudobulk=True)
tester.edger(design="~region+donor_id", reference="Control", groups=["IZ", "RZ", "FZ"])
table = tester.get_dge["EdgeR"]
print(table[table.GeneName =="TCF15"])
#       GeneName  statistic    log2fc      pval      padj group
# 4254     TCF15   9.973417 -1.189486  0.004577  0.026414    IZ
# 16258    TCF15   0.238949 -0.201091  0.629825  0.950816    RZ
# 10035    TCF15   2.137417 -0.602142  0.157944  0.386456    FZ

adata_clean = adata_concat[adata_concat.obs.region.isin(["Control", "IZ", "RZ", "FZ"])]

# * Figure 1f
do.pl.barplot(
    adata_clean, 'region', 'TCF15', batch_key='donor_id', xticks_order=["Control", "RZ", "IZ", "FZ"],
    reference="Control", groups=["RZ", "IZ", "FZ"], groups_pvals=[0.9508, 0.026414, 0.386456], figsize=(3, 5),
    line_offset=0.1, txt_size=10, xticks_rotation=45, path=out_path,
    filename="Barplot_TCF15_HumanRegions_Statistics.pdf"
)

# endregion

