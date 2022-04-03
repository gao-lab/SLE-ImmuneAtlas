import pandas as pd 
import numpy as np
import scanpy as sc

marker_list  =  ['IL7R','CCR7','TCF7','LEF1',  # Naive T
                 'CD3D','CD3E','CD3G',         # T cell marker 
                 'PRF1','EOMES',               # activated T cell marker also NK
                 'S100A4',                     # Memory CD4+ T also 'IL7R+'
                 'CD8A','CD8B', 'CD4',         # CD4 and CD8+ T
                 'GNLY','GZMK','GZMA',         # T and NK
                 'NKG7','KLRC1','KLRD1','KLRF1',# NK
                 'CD68','S100A9','CSF1R','LYZ',# Mono/Macro
                 'FCGR3A','MS4A7',             # FCGR3A+ Mono
                 'CD14',                       # CD14+ Mono
                 'MS4A1','CD79A',              # B
                 'JCHAIN','IGKC','MZB1',       # plasma
                 'CST3','IRF7','ITGAX','CD86', # DC
                 'LILRA4','CLEC4C',            # pDCs
                 'CD1C','FCER1A',              # mDCs also 'LYZ'
                 'PPBP', 'GP9','PF4',          # Platelet
                 'CPA3',                       # MAST
                 'GATA1',                      # red blood cell
                 'PTPRC',                      # immune
                 'PCNA',                       # profiling
                 'CD34'                        # hspc
                 # 'AIDIA'                     # regrouping 
                ]


IFN_genes = ["ABCE1","ADAR","BST2","CACTIN","CDC37","CNOT7","DCST1","EGR1","FADD","GBP2","HLA-A","HLA-B","HLA-C",
             "HLA-E","HLA-F","HLA-G","HLA-H","HSP90AB1","IFI27","IFI35","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2",
             "IFITM3","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7",
             "IFNA8","IFNAR1","IFNAR2","IFNB1","IKBKE","IP6K2","IRAK1","IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8",
             "IRF9","ISG15","ISG20","JAK1","LSM14A","MAVS","METTL3","MIR21","MMP12","MUL1","MX1","MX2","MYD88","NLRC5",
             "OAS1","OAS2","OAS3","OASL","PSMB8","PTPN1","PTPN11","PTPN2","PTPN6","RNASEL","RSAD2","SAMHD1","SETD2","SHFL",
             "SHMT2","SP100","STAT1","STAT2","TBK1","TREX1","TRIM56","TRIM6","TTLL12","TYK2","UBE2K","USP18","WNT5A","XAF1",
             "YTHDF2","YTHDF3","ZBP1"]


###############################################################
#
#  Funtions
#
###############################################################

def do_harmony(adata, n_top_genes = 2500,batch = 'sample', theta = 1, resolution = 0.8,hvg_batch = True,max_iter_harmony = 20):
    '''
    do harmony in scanpy
    '''
    import scanpy.external as sce
    if hvg_batch:
        sc.pp.highly_variable_genes(adata, flavor ='seurat_v3',n_top_genes =n_top_genes, batch_key = batch) # must specify n_top_genes when using seurat v3
    else:
        sc.pp.highly_variable_genes(adata, flavor ='seurat_v3',n_top_genes =n_top_genes)
        

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    sce.pp.harmony_integrate(adata, batch,theta=theta,max_iter_harmony = max_iter_harmony)
    
    # default parameters
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50,use_rep='X_pca_harmony' )
    sc.tl.leiden(adata,resolution = resolution)

    sc.tl.umap(adata)
    
    return adata
