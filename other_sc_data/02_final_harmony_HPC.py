import os
# os.chdir('/rd2/user/xiacr/sle/other_sc_data/')
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import importlib

import glob
import anndata as ad
import collections

import scanpy.external as sce
from mycolorpy import colorlist as mcp

import utils
importlib.reload(utils)

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3) 
sc.logging.print_header() 
sc.settings.set_figure_params(dpi=160, facecolor='white', fontsize=10) sc._settings.ScanpyConfig.n_jobs = 36

nc_sle_flare = sc.read_h5ad('./h5ad/nc_19_sle_flare.h5ad') 
nc_sle= sc.read_h5ad('./h5ad/nc_19_sle.h5ad')

nc_sle_flare.obs['dataset'] = 'nc_19'
nc_sle_flare.obs['group'] = 'sle_flare'
nc_sle.obs['dataset'] = 'nc_19'
nc_sle.obs['group'] = 'sle'
nc_sle.obs.loc[nc_sle.obs[(nc_sle.obs['index'].isin(['sample12','sample18','sample8','sample22','sample30','sample34','sample45','sample26','sample35','sample27','sample31','sample39']))].index.tolist(),'group'] = 'hc'
nc_sle_flare.obs['sample'] =  nc_sle_flare.obs['dataset'].astype(str) + '_' + nc_sle_flare.obs['value'].astype(str)
nc_sle.obs['sample'] = nc_sle.obs['dataset'].astype(str) + '_' + nc_sle.obs['index'].astype(str)

core_data =  sc.read_h5ad('./h5ad/04-pbmc_all_anno_modify_meta.h5ad')

core_data.obs['dataset'] = 'core_dataset'
core_data.obs['group'] = core_data.obs['group'].astype(str)
core_data.obs.loc[core_data.obs[(core_data.obs['treatment'].isin(['untreated']))].index.tolist(),'group'] = 'sle_flare'
core_data.obs.loc[core_data.obs[(core_data.obs['treatment'].isin(['treated']))].index.tolist(),'group'] = 'sle_treated'
core_data.obs.loc[core_data.obs[(core_data.obs['treatment'].isin(['HC']))].index.tolist(),'group'] = 'hc'
core_data.obs['sample'] = core_data.obs['dataset'].astype(str) + '_' + core_data.obs['orig.ident'].astype(str)

nbt_17= sc.read_h5ad('./h5ad/nbt_17_meta.h5ad')
pnas_19 = sc.read_h5ad('./h5ad/pnas_19_meta.h5ad')
ni_20 =  sc.read_h5ad('./h5ad/ni_20_meta.h5ad')
ebi_21 = sc.read_h5ad('./h5ad/ebi_21_meta.h5ad')

pnas_19.obs['group'] = 'sle' 
nbt_17.obs['dataset'] = 'nbt_17'
pnas_19.obs['dataset'] = 'pnas_19'
ni_20.obs['dataset'] = 'ni_20'
ebi_21.obs['dataset'] = 'ebi_21'

core_data.obs['dataset'] = core_data.obs['dataset'].astype(str)
nc_sle_flare.obs['dataset'] = nc_sle_flare.obs['dataset'].astype(str)
nc_sle.obs['dataset'] = nc_sle.obs['dataset'].astype(str)

adata = core_data.concatenate([nc_sle,nc_sle_flare, ebi_21,ni_20,pnas_19,nbt_17],join="outer")

adata.write_h5ad('./output/01-pbmc_concat_raw.h5ad')
del nc_sle_flare; del nc_sle; del core_data; del ebi_21; del ni_20; del pnas_19; del nbt_17

# -----------------------------------------------------------
print('step2')

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=20)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.highly_variable_genes(adata, flavor ='seurat_v3',n_top_genes =2500, batch_key ='sample') # must specify n_top_genes when using seurat v3

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata # You can get back an AnnData of the object in .raw by calling .raw.to_adata().
adata = adata[:, adata.var.highly_variable]

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

sce.pp.harmony_integrate(adata, 'sample',theta=1) # theta is the key parameter and I find set 1 works well

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50,use_rep='X_pca_harmony' ) sc.tl.leiden(adata,resolution = 0.8)  
sc.tl.umap(adata)

adata.write_h5ad('./output/02-pbmc_concat_harm_anno.h5ad')