library(cowplot)
library(harmony)
library(Seurat)
library(future)

setwd('/data/sle')
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./final//seurat/t_cell/03-NK_cell_raw_harm.rdata')


DimPlot(nk_filter, label = T)
VlnPlot(nk_filter, features = feats)
DotPlot2(nk_filter, marker_list = marker_list)
DotPlot2(nk_filter, marker_list = t_cd8_all_marker)
VlnPlot(nk_filter, features = c('FCGR3A','NCAM1','XCL1','IFIT3','TNF','CXCR4','CD3D','KLRF1','KLRB1'),stack = T)

nk_filter <- subset(nk_filter, idents =c(9,10), invert =T)

#---------------------------------- Anno ---------------------------------------
nk_filter$subtype <- 'NK'
nk_filter$subtype[which(nk_filter$seurat_clusters %in% c(2,4,5))] <- 'NK.cd56+'
DimPlot(nk_filter, label = T, group.by = 'subtype',cols = get_color(2,set_len =4), pt.size = 0.3) + NoAxes()

#---------------------------------- Save ---------------------------------------
save(nk_filter,file = './final/seurat/t_cell/04-NK_Tcell_filter_anno.rdata')
