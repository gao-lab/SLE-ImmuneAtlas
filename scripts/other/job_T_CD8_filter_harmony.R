library(tidyverse)
library(harmony)
library(Seurat)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/t_cell/CD8_tcell_filter_raw.rdata')

#------------------------------ Re Analysis ------------------------------------
t.cd8.filter <- do_seurat(t.cd8.filter)


#--------------------------- Harmony Batch Remove ------------------------------
# run Harmony ------
t.cd8.filter <- t.cd8.filter %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)%>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

#------------------------------- Save Results-----------------------------------
save(t.cd8.filter, file = paste0(output_path,'cd8_tcell_filter_harmony_01.rdata'))
