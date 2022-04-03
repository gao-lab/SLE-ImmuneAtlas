library(tidyverse)
library(harmony)
library(Seurat)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/t_cell/CD4_tcell_filter_raw.rdata')

#------------------------------ Re Analysis ------------------------------------
t.cd4.filter <- do_seurat(t.cd4.filter)


#--------------------------- Harmony Batch Remove ------------------------------
# run Harmony ------
t.cd4.filter <- t.cd4.filter %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)%>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

#------------------------------- Save Results-----------------------------------
save(t.cd4.filter, file = paste0(output_path,'CD4_tcell_filter_harmony_01.rdata'))
