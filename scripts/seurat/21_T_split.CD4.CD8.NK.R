library(cowplot)
library(harmony)
library(RISC)
library(future)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')
# load('./final/seurat/t_cell/03_t_cell_anno_harm.rdata')
load('./final/seurat/t_cell/02-t_cell_raw_harm.rdata')
#---------------------------- Split and Harmony --------------------------------
DimPlot(tcell_harm, group.by = 'subtype',label = T) + DimPlot(tcell_harm, label = T)
# 11,15,19 is doublet; 13,18 is proliferation
cd4 <- subset(tcell_harm, idents = c(1,2,6,7,12,21))
cd8 <- subset(tcell_harm, idents = c(0,3,4,9,13,18))

back_run(do_harmony, out_name = 'cd4_harm', job_name = 'cd4', from_begin = T,
         seu_obj = cd4, harmony_slot = 'orig.ident',max.iter = 30, res = c(0.8))
back_run(do_harmony, out_name = 'cd8_harm', job_name = 'cd8', from_begin = T,
         seu_obj = cd8, harmony_slot = 'orig.ident',max.iter = 30, res = c(0.8))
rm(cd4, cd8)

DimPlot(cd4_harm, label = T) + DimPlot(cd8_harm, label = T)

#----------------------------Vis Marker Genes ----------------------------------
# -------- CD4 --------  
# cluster 10 and 11 is CD8 T 
DimPlot(cd4_harm, label = T) + DimPlot(cd4_harm, group.by = 'group')
DoHeatmap(subset(cd4_harm, downsample = 2000), features = c('CD3D','CD4','CD8A','CD8B'), label = F)
DotPlot(cd4_harm, features = c('CD3D','CD4','CD8A','CD8B'))
VlnPlot(cd4_harm, features = c('CD3D','CD4','CD8A','CD8B'), stack = T)
DotPlot2(cd4_harm, marker_list  = marker_list)
VlnPlot(cd4_harm, features = feats,stack = F)

# -------- CD8 --------  
# cluster 11 is nk cell, cluster 8 is other T, cluster 4 9 may be doublet  
DimPlot(cd8_harm, label = T) + DimPlot(cd8_harm, group.by = 'group')
# DimPlot(cd8_harm, label = T) + DimPlot(cd8_harm, group.by = 'subtype')
DoHeatmap(subset(cd8_harm, downsample = 2000), features = c('CD3D','CD4','CD8A','CD8B'), label = F)
FeaturePlot(cd8_harm, features = c('CD3D','CD4','CD8A','CD8B'), ncol = 2)
VlnPlot(cd8_harm, features = c('CD3D','CD4','CD8A','CD8B'), stack = T)
VlnPlot(cd8_harm, features = feats,stack = F)
DotPlot(cd8_harm, features = t_rare_marker)
DotPlot2(cd8_harm, marker_list = marker_list)

#----------------------------- Ratio of cluster --------------------------------
# Cluster by group
cd8_harm@meta.data %>% group_by(orig.ident, seurat_clusters) %>% 
  summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
  mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_harm@meta.data[, c(1,4,5,6)]  %>%  distinct()) %>%
  ggpubr::ggboxplot(x = 'group',y = 'Ratio', color = 'group',
                    palette = c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters, scales = "free") + stat_compare_means(method = 't.test')

#------------------------------ Re-Subset ALL T cells --------------------------
# CD 4
CD4_barcode <- cd4_harm@meta.data %>% 
  filter(!seurat_clusters %in% c(10, 11)) %>% row.names()
CD4_rescue <- cd8_harm@meta.data %>% 
  filter(seurat_clusters %in% c(8)) %>% row.names()
CD4_barcode_all <- c(CD4_barcode, CD4_rescue)

# CD 8
CD8_barcode <- cd8_harm@meta.data %>% 
  filter(!seurat_clusters %in% c(8)) %>% row.names()
CD8_rescue <- cd4_harm@meta.data %>% 
    filter(seurat_clusters %in% c(10, 11)) %>% row.names()
CD8_barcode_all <- c(CD8_barcode, CD8_rescue)

#------------------------------- Save Files ------------------------------------
# subset via barcodes
cd4_filter <- subset(tcell_harm, cells = CD4_barcode_all)
cd8_filter <- subset(tcell_harm, cells = CD8_barcode_all)
back_run(do_harmony,out_name = 'cd4_filter', job_name = 'cd4_filter', from_begin = T,
         seu_obj = cd4_filter, harmony_slot = 'orig.ident', max.iter = 30, res = c(0.8))
back_run(do_harmony,out_name = 'cd8_filter', job_name = 'cd8_filter',from_begin = T,
         seu_obj = cd8_filter, harmony_slot = 'orig.ident', max.iter = 30, res = c(0.8))
cd4_filter <- subset(cd4_filter, idents = c(11, 12), invert = T)
save(cd4_filter,file = 'final/addtional/seurat/t_cell/03-CD4_Tcell_filter_harm.rdata')
save(cd8_filter,file = 'final/addtional/seurat/t_cell/03-CD8_Tcell_filter_harm.rdata')
rm(cd4_harm, cd8_harm)
