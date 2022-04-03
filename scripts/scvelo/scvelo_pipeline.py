import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
# import glob
# import anndata as ad
# import collections
# from tqdm import tqdm
import argparse

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo') 

parser = argparse.ArgumentParser(description="scvelo for special cell type")
parser.add_argument('--name',type=str, required=True,
                   help='cell type, will be prefix of out')
parser.add_argument('--meta',type=str, required=True,
                   help='seurat meta of the cell used in scVelo')
parser.add_argument('-n','--jobs', default=32,type=int,
                   help='(optional) max cores for computing in some steps')
parser.add_argument('-g','--group', default='subtype',type=str,
                   help='(optional) varaible in seurat meta to group the cell subtype')
parser.add_argument('-s','--save', action='store_true',
                   help='(optional) if save the adata')
results = parser.parse_args()
name = results.name
meta = results.meta
n_jobs = results.jobs
group = results.group

out_pic1 = name + '_Velocity_pic.pdf'
out_pic2 = name + '_Velocity_Paga_pic.pdf'
out_pic3 = name + '_Velocity_grid_pic.pdf'

# os.chdir('/rd2/user/xiacr/sle/')
# results_file = 'output_file/scvelo/' 
adata =sc.read_h5ad('/rd2/user/xiacr/sle/final/scvelo/all_cell_exclude_berry.h5ad')
cell_meta = pd.read_csv(meta)
cell_meta['1'], cell_meta['2'] = cell_meta['Row.names'].str.split('_', 1).str
cell_meta['barcode'] = cell_meta['1']+ '_' + cell_meta['orig.ident']

# check the cell meta 
if not set(['UMAP_1','UMAP_1','barcode']).issubset(cell_meta.columns):
    # print()
    raise RuntimeError('ERROR: seurat cell meta is not valid, please check it!')

adata_sub = adata[adata.obs['barcode'].isin(cell_meta['barcode'] )]
adata_sub.obs = adata_sub.obs.rename_axis("CellID").reset_index()
adata_sub.obs = adata_sub.obs.merge(cell_meta, how='left', left_on='barcode', right_on='barcode')
adata_sub.obs = adata_sub.obs.set_index('CellID')

# scv.pp.filter_genes(adata_sub, min_shared_counts=20)
# scv.pp.normalize_per_cell(adata_sub)
# scv.pp.filter_genes_dispersion(adata_sub, n_top_genes=2000)
# scv.pp.log1p(adata_sub)

scv.pp.filter_and_normalize(adata_sub, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_sub, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata_sub,groupby='batch')
scv.tl.velocity_graph(adata_sub, n_jobs =n_jobs)

adata_sub.obsm['X_harmony'] =  adata_sub.obs[['UMAP_1','UMAP_2']].to_numpy()
scv.pl.velocity_embedding_stream(adata_sub, basis='harmony',color=group,
                                 save= out_pic1)
scv.pl.velocity_embedding_grid(adata_sub, dpi=120,color=group, 
                               basis='harmony',scale=0.05, arrow_length =3, save= out_pic3)

adata_sub.uns['neighbors']['distances'] = adata_sub.obsp['distances']
adata_sub.uns['neighbors']['connectivities'] = adata_sub.obsp['connectivities']

scv.tl.paga(adata_sub, groups=group)
# df = scv.get_df(adata_sub, 'paga/transitions_confidence', precision=2).T
# df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata_sub, basis='harmony', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save=out_pic2)

if  results.save:
    adata_sub.write_h5ad(name+'_scVelo.h5ad')
