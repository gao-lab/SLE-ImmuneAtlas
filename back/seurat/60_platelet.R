library(cowplot)
library(harmony)
# library(RISC)
library(Seurat)
library(tidyverse)
library(ggpubr)

setwd('/data/sle')
output_path <- './output_file/seurat/platelet/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/platelet/platelet_subset_01.rdata')

#------------------------------ Re Analysis ------------------------------------
platelet <- do_seurat(platelet)
DimPlot(platelet, label = T) + DimPlot(platelet, group.by = 'orig.ident')


#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
platelet_harmony <- platelet %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

# save(platelet_harmony, file = paste0(output_path, 'platelet_harmony_14iter.rdata'))

# DimPlot(platelet_harmony, label = T) + DimPlot(platelet_harmony, group.by = 'orig.ident')
plot_scdata(platelet_harmony, color_by = "treatment", pal_setup = 'Set2') + 
  plot_scdata(platelet_harmony, color_by = "subtype", pal_setup = 'Set1')


#----------------------------- Vis Marker Genes --------------------------------
VlnPlot(platelet_harmony, features = feats)

FeaturePlot(platelet_harmony, features = c('CD58','CD69','IFITM1'), pt.size = 0.7,
            cols = c('grey','red'), ncol = 3) &NoAxes()

DotPlot(platelet_harmony, features =c('CD58','CD69','IFITM1','PRKRA'))
DotPlot(platelet_harmony, features =c('CD58','CD69','IFITM1','PRKRA') , 
        group.by = 'group')
DotPlot2(platelet_harmony, marker_list = c('CD58','CD69','IFITM1'), group.by = 'subtype')
VlnPlot(platelet_harmony, features =c('CD58','CD69','IFITM1'), stack = T , group.by = 'subtype')
marker_platelet_slepah <- FindMarkers(platelet_harmony, group.by = 'group', 
                                      ident.1 = 'SLE_pah')
marker_platelet_slepah <- marker_platelet_slepah%>% filter(p_val_adj<0.05)

# by cell type annotation result 
DotPlot(platelet_harmony, features = c('CD58','CD69','IFITM1','PRKRA'),
        group.by = 'subtype')

#---------------------------- Anno the Sub type --------------------------------
platelet_harmony$subtype <- 'Platelet'
platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c(6))] <- 'Platelet.CD58+CD69+'
platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c(2,4))] <- 'Platelet.IFITM1+'
# platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c())] <- ''

table(platelet_harmony$subtype)
DimPlot(platelet_harmony, group.by = 'subtype', label = T) + NoLegend()


#-------------------------- Compare at Sample Level ----------------------------
# plate to all pbmc

# By Group
platelet_harmony@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free") 

# By Treatment
platelet_harmony@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'treatment',
                    palette =c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  facet_wrap(~seurat_clusters,scales = "free") 


#-------------------------------- Save Files -----------------------------------
save(platelet_harmony, file = paste0(output_path,'platelet_harmony_anno.rdata'))
save(platelet_harmony, file = './output_file/seurat/platelet/final_platelet.rdata')
