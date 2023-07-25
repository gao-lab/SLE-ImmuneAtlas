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
# 11,15,19 is doublet; 13,18 is prolifing 
cd8 <- subset(tcell_harm,idents = c(0,3,4,6,14,17) )
cd4 <- subset(tcell_harm,idents = c(1,2,8,9,12))
nk <- subset(tcell_harm,idents = c(5,16))
other_T <- subset(tcell_harm,idents = c(7,10))
prolife_T <- subset(tcell_harm,idents =c(13,18))

back_run(do_harmony,out_name = 'cd4_harm', job_name = 'cd4', from_begin = T,
         seu_obj = cd4, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'cd8_harm', job_name = 'cd8',from_begin = T,
         seu_obj = cd8, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'nk_harm', job_name = 'nk',from_begin = T,
         seu_obj = nk, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'other_T_harm', job_name = 'other_T',from_begin = T,
         seu_obj = other_T, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'prolife_T_harm', job_name = 'prolife_T',from_begin = T,
         seu_obj = prolife_T, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
rm(cd4,cd8,nk, other_T,prolife_T)

DimPlot(other_T_harm, label = T) + DimPlot(nk_harm, label = T)
DimPlot(cd4_harm, label = T) + DimPlot(cd8_harm, label = T)
DimPlot(other_T_harm, cells.highlight = tcell_harm@meta.data %>% filter(subtype == 'T.yd') %>% rownames())
DimPlot(other_T_harm, cells.highlight = tcell_harm@meta.data %>% filter(subtype == 'MAIT') %>% rownames())

#----------------------------Vis Marker Genes ----------------------------------
# -------- CD4 --------  
# cluster 10 and 11 is CD8 T 
DimPlot(cd4_harm, label = T) + DimPlot(cd4_harm, group.by = 'group')
DoHeatmap(subset(cd4_harm, downsample = 2000), features = c('CD3D','CD4','CD8A','CD8B'), label = F)
DotPlot(cd4_harm, features = c('CD3D','CD4','CD8A','CD8B'))
VlnPlot(cd4_harm, features = c('CD3D','CD4','CD8A','CD8B'),stack = T)
DotPlot2(cd4_harm, marker_list  = marker_list)
VlnPlot(cd4_harm, features = feats,stack = F)

# -------- CD8 --------  
# cluster 11 is nk cell, cluster 8 is other T, cluster 4 9 may be doublet  
DimPlot(cd8_harm, label = T) + DimPlot(cd8_harm, group.by = 'group')
# DimPlot(cd8_harm, label = T) + DimPlot(cd8_harm, group.by = 'subtype')
DoHeatmap(subset(cd8_harm, downsample = 2000), features = c('CD3D','CD4','CD8A','CD8B'), label = F)
FeaturePlot(cd8_harm, features = c('CD3D','CD4','CD8A','CD8B'), ncol = 2)
VlnPlot(cd8_harm, features = c('CD3D','CD4','CD8A','CD8B'),stack = T)
VlnPlot(cd8_harm, features = feats,stack = F)
DotPlot(cd8_harm, features = t_rare_marker)
DotPlot2(cd8_harm, marker_list  = marker_list)

# -------- NK --------  
# cluster 6 is CD8 T cell, cluster 8 is doublet 
DimPlot(nk_harm, label = T) + DimPlot(nk_harm, group.by = 'group')
VlnPlot(nk_harm, features = feats,stack = F)
DotPlot2(nk_harm, marker_list  = marker_list)

# -------- other --------  
# cluster 6 is CD8 T cell, cluster 8 is doublet 
DimPlot(other_T_harm, label = T) + DimPlot(other_T_harm, group.by = 'group')
DotPlot2(other_T_harm, marker_list  = marker_list)
DotPlot2(other_T_harm, marker_list  = t_rare_marker)
VlnPlot(other_T_harm, features = feats,stack = F, ncol = 2)

#----------------------------- Ratio of cluster --------------------------------
# Cluster by group
cd8_harm@meta.data  %>% group_by(orig.ident,seurat_clusters) %>% 
  summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
  mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_harm@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free")+ stat_compare_means(method = 't.test')

#------------------------------- Finder DEGs -----------------------------------
plan("multiprocess", workers = 16)
options(future.globals.maxSize= 120 * 1000 * 1024 ^2 ) 
marker_cd8_c4 <- FindMarkers(cd8_harm, ident.1 = 4, only.pos = T, min.pct = 0.2 ,logfc.threshold = 0.25)
marker_cd8_c8 <- FindMarkers(cd8_harm, ident.1 = 8, only.pos = T, min.pct = 0.2 ,logfc.threshold = 0.25)
marker_cd8_c9<- FindMarkers(cd8_harm, ident.1 = 9, only.pos = T, min.pct = 0.2,logfc.threshold = 0.25)
marker_cd8_c11<- FindMarkers(cd8_harm, ident.1 = 11, only.pos = T, min.pct = 0.2,logfc.threshold = 0.25)
back_run(FindMarkers, out_name = 'marker_cd8_c4',job_name = 'marker_cd8_c4',
         cd8_harm, ident.1 = 4, only.pos = T, min.pct = 0.2 ,logfc.threshold = 0.25)
back_run(FindMarkers, out_name = 'marker_cd8_c8',job_name = 'marker_cd8_c8',
         cd8_harm, ident.1 = 8, only.pos = T, min.pct = 0.2 ,logfc.threshold = 0.25)
back_run(FindMarkers, out_name = 'marker_cd8_c9',job_name = 'marker_cd8_c9',
         cd8_harm, ident.1 = 9, only.pos = T, min.pct = 0.2,logfc.threshold = 0.25)
back_run(FindMarkers, out_name = 'marker_cd8_c11',job_name = 'marker_cd8_c11',
         cd8_harm, ident.1 = 11, only.pos = T, min.pct = 0.2,logfc.threshold = 0.25)

#------------------------------ Re-Subset ALL T cells --------------------------
# CD 4
CD4_barcode <- cd4_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(!seurat_clusters %in% c(10,11)) %>% row.names()
CD4_rescue <- cd8_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(seurat_clusters %in% c(5)) %>% row.names()
CD4_barcode_all <- c(CD4_barcode,CD4_rescue)

# CD 8
CD8_rescue1 <- cd4_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(seurat_clusters %in% c(10,11)) %>% row.names()
CD8_rescue2 <- nk_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(seurat_clusters %in% c(6)) %>% row.names()
CD8_barcode <- cd8_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(!seurat_clusters %in% c(4,5,8,9,11)) %>% row.names()
CD8_barcode_all <- c(CD8_barcode,CD8_rescue1,CD8_rescue2)

# NK
NK_rescue <- cd8_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(seurat_clusters %in% c(11)) %>% row.names() 
nk_barcode <- nk_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
  filter(!seurat_clusters %in% c(6,8)) %>% row.names()
nk_barcode_all <- c(nk_barcode,NK_rescue)

# Other T cell
# other_rescue <- cd8_harm@meta.data %>% filter(scrublet_doublet == 'False') %>% 
#   filter(seurat_clusters %in% c(8)) %>% row.names()
# other_barcode_all <- c(Cells(other_T_harm), other_rescue)

#------------------------------- Save Files ------------------------------------
# subset via barcode
cd4_filter <- subset(tcell_harm, cells = CD4_barcode_all)
cd8_filter <- subset(tcell_harm, cells = CD8_barcode_all)
nk_filter <- subset(tcell_harm, cells = nk_barcode_all)
other_T_filter <- other_T_harm
prolife_T_harm <- subset(prolife_T_harm, idents = 0 ,invert =T )
back_run(do_harmony,out_name = 'cd4_filter', job_name = 'cd4_filter', from_begin = T,
         seu_obj = cd4_filter, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'cd8_filter', job_name = 'cd8_filter',from_begin = T,
         seu_obj = cd8_filter, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'nk_filter', job_name = 'nk_filter',from_begin = T,
         seu_obj = nk_filter, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
back_run(do_harmony,out_name = 'other_T_filter', job_name = 'other_T_filter',from_begin = T,
         seu_obj = other_T_filter, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1,0.8))
# rm(cd4,cd8,nk, other_T)

save(cd4_filter,file =  'final/seurat/t_cell/03-CD4_Tcell_raw_harm.rdata')
save(cd8_filter,file =  'final/seurat/t_cell/03-CD8_Tcell_raw_harm.rdata')
save(nk_filter,file =  'final/seurat/t_cell/03-NK_cell_raw_harm.rdata')
save(other_T_filter,file =  'final/seurat/t_cell/03-other_Tcell_raw_harm.rdata')
save(prolife_T_harm, file = 'final/seurat/t_cell/03-prolife_Tcell_raw_harm.rdata')
rm(cd4_filter_harm,cd4_harm,cd8_harm,nk_harm,other_T_harm)
