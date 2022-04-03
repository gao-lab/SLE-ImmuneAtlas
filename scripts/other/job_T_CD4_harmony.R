library(tidyverse)
library(harmony)
library(Seurat)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/t_cell/tcell_harmony_anno.rdata')
#------------------------------ Re Analysis ------------------------------------
t.cd4 <- subset(tcell, idents = c(1,2,6,7,10,15))
# t.cd8 <- subset(tcell, idents = c(0,3,5,9,11,12,14,16,18,20,21))

rm(tcell)
t.cd4 <- do_seurat(t.cd4)
# t.cd8 <- do_seurat(t.cd8)


#--------------------------- Harmony Batch Remove ------------------------------
# run Harmony ------
t.cd4 <- t.cd4 %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)%>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

# t.cd8 <- t.cd8 %>% 
#   RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)%>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.8) %>% 
#   identity()


#------------------------------- Save Results-----------------------------------
save(t.cd4, file = paste0(output_path,'CD4_tcell_harmony_01.rdata'))
# save(t.cd8, file = paste0(output_path,'CD8_tcell_harmony_01.rdata'))
