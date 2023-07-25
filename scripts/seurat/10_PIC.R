# PIC cells
library(cowplot)
library(harmony)
# library(RISC)
library(Seurat)
library(tidyverse)
library(ggpubr)

setwd('/data/sle')
output_path <- './output_file/seurat/pic/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/pic/pic_subset_01.rdata')
load('./output_file/seurat/mono_dc/pic_mono_dc_harmony.rdata')
load('./output_file/seurat/b_cell/pic_bcell_harmony_11iter.rdata')

#------------------------------ Re Analysis ------------------------------------
pic_all <- merge(pic, c(mono_dc_pic, bcell_pic))
rm(pic,mono_dc_pic,bcell_pic)
pic_all <- do_seurat(pic_all)
DimPlot(pic_all, label = T) + DimPlot(pic_all, group.by = 'orig.ident')
save(pic_all, file = paste0(output_path,'pic_all_recluster.rdata'))

#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
pic_all_harmony <- pic_all %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

save(pic_all_harmony, file = paste0(output_path, 'pic_all_harmony_14iter.rdata'))

DimPlot(pic_all_harmony, label = T) + DimPlot(pic_all_harmony, group.by = 'orig.ident')
