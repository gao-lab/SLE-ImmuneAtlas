
####---------------------------source----------------------------------####
b_marker_list.traditional <- c("CD79A",'CD19','CD74',
                   'MME','IGHD','IGHM','CD24','CD38',
                   'CD27','MS4A1','CD69','CD83','CD86',
                   # 'CD14','CD5','CD23',
                   'CD2','CD3E','CD3D','CD3G','CD8A','CD8B','CD4',
                   'SDC1','IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','JCHAIN','IGKC','MZB1','XBP1','BLIMP1',
                   # "PCNA","MKI67","SERPINA9","AICDA","CD44","PTPRC"
                   "NONE"
)
marker_list <- c('IL7R','CCR7','TCF7','LEF1',  # Naive T
                 'CD3D','CD3E','CD3G',         # T cell marker 
                 'PRF1',                       # activated T cell marker also NK
                 'S100A4',                     # Memory CD4+ T also 'IL7R+'
                 'CD8A','CD8B', 'CD4',         # CD4 and CD8+ T
                 'GNLY','GZMK','GZMA',         # T and NK
                 'NKG7','KLRC1','KLRD1','KLRF1',# NK
                 'CD68','S100A9','CSF1R','LYZ',# Mono/Macro
                 'FCGR3A','MS4A7','CD16',      # FCGR3A+ Mono
                 'CD14',                       # CD14+ Mono
                 'MS4A1','CD79A',              # B
                 'JCHAIN','IGKC','MZB1',       # plasma
                 'CST3','IRF7','ITGAX','CD86', # DC
                 'LILRA4','CLEC4C',            # pDCs
                 'CD1C','FCER1A',               # mDCs also 'LYZ'
                 'PPBP', 'GP9','PF4',          # Platelet
                 'CPA3',                       # MAST
                 'GATA1',                      # red blood cell
                 'PTPRC',                      # immune
                 # 'MK167','PCNA',               # profiling
                 # 'AIDIA'                       # regrouping 
                 "NONE"
)
tlb_marker_list <- c('IL7R','SELL','LEF1','TCF7','MAL','CCR7','CCR6',
                     'GNLY','NKG7','CCL4','CCL5','GZMA','GZMB','GZMH','CD3D','CD3E',
                     'CD79A','CD7A','SDC1','CD138','IGHG1','IGHG2','IGHG3','IGHA1','IGHA2','IGHA3',
                     'IGHM','IGHD','CD27','MS4A1','JCHAIN','IGKC','MZB1')
tlb_marker_list.T <- c('CD4','IL7R','SELL','LEF1','TCF7','MAL','CCR7','CCR6',
                       'GNLY','NKG7','CCL4','CCL5','GZMA','GZMB','GZMH')

plasma_marker <- c('SDC1','MZB1','PRDM1','SLAMF7','ABCB1','XBP1','BCMA','CD38','CXCR4',
                   'CD27','CD19','CD79A',
                   'IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','JCHAIN','IGKC')

t_cell_marker <- c('CCR7','LEF1','SELL','TCF7','CD27','CD28','S1PR1', # naive marker 
                   'ITGAE',  # TRM marker 
                   'CCR4','CCR5', # TCM marker
                   'NCAM1','FCGR3A','NKG7','KLRC1','KLRD1','KLRF1', # NK,
                   'CD44','CD69', # activate T 
                   'FOXP3','IL2RA',  # Treg
                   'PDCD1','LAYN','HAVCR2', # exhausted T
                   'BCL6','ICOS','CD40LG',  # Tfh
                   'CD1D'  # NKT
)

TLB_maual_marker <- c('IL7R','CCR7','TCF7',
                      'NKG7','GNLY','GZMA',
                      'MZB1','CST3','IGJ')

stem_cell <- c('CD34','AVP')
# for plot 
feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo")
gc_b_marker <- c('RGS13','MEF2B','SEMA4A','TCL1A','HMGN1','HMGB2','STMN1','TOP2A','CD40')
