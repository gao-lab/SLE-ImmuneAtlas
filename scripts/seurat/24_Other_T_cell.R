library(future)

setwd('/data/sle')
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./final//seurat/t_cell/03-other_Tcell_raw_harm.rdata')

DotPlot2(other_T_filter, marker_list = marker_list)
DotPlot2(other_T_filter, marker_list = t_rare_marker)
DimPlot(other_T_filter ,label = T)
VlnPlot(other_T_filter, features = t_rare_marker[-c(6,11)], stack = T)

#---------------------------------- Anno ---------------------------------------
other_T_filter$subtype <- 'T.yd'
other_T_filter$subtype[which(other_T_filter$seurat_clusters %in% c(2))] <- 'T.MAIT'
DimPlot(other_T_filter, label = T, group.by = 'subtype',cols = get_color(2,set_len =4), pt.size = 0.3) + NoAxes()

#---------------------------------- Save ---------------------------------------
save(other_T_filter,file = './final/seurat/t_cell/04-Other_Tcell_filter_anno.rdata')
