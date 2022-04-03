# RISC on all T cell and NK cell in SLE data set 

library(RISC)
# library(Seurat)
# library(future)
# library(stringr)
# library(tidyverse)

#------------------------------ Input ------------------------------------------
setwd('/data/sle/')
load('./tmp/RISC_input_t_cell_list.rdata')
output_path <- './output_file/RISC/'


#------------------------------ Code ------------------------------------------
risc_list <- list()
for(i in c(1:length(t_cell_list))){
  data = t_cell_list[[i]]@assays$RNA@counts
  coldat <- t_cell_list[[i]]@meta.data
  rawdata <- data.frame(Symbol = rownames(data), RNA = "Gene Expression", row.names = rownames(data))
  risc_list[[i]] <- readscdata(count = data, cell = coldat, gene = rawdata)
}

for(i in risc_list){
  print(FilterPlot(i))
}

process0 <- function(obj0){
  # Filter cells and genes
  obj0 = scFilter(obj0, min.UMI = 500, max.UMI = 40000, min.gene = 200, min.cell = 3)
  # Normalize the raw counts
  obj0 = scNormalize(obj0)
  # Find highly variable genes
  obj0 = scDisperse(obj0)
  # print(length(obj0@vargene))
  return(obj0)
}

# options(repr.plot.width=6, repr.plot.height=5)
for(i in c(1:length(risc_list))){
  risc_list[[i]] <- process0(risc_list[[i]])
}

var0 = risc_list[[1]]@rowdata$Symbol
for(i in c(1:length(risc_list)) ){
  var0 = intersect(var0,risc_list[[i]]@rowdata$Symbol )
}

InPlot(risc_list,var.gene = var0, Std.cut = 0.99, ncore = 32, minPC = 16, nPC = 40)

data0 = scMultiIntegrate(
  objects = risc_list, eigens = 20, add.Id = NULL, var.gene = var0,
  method = "RPCI", align = 'OLS', npc = 50, adjust = TRUE,
  ncore = 32, do.fast = "AUTO"
)

data0 = scUMAP(data0, npc = 18, use = 'PLS')
DimPlot(data0, slot = "cell.umap", colFactor = "cell_type")

data0 = scCluster(data0, slot = "cell.pls", neighbor = 3, method = "louvain", npc = 20)
print(table(data0@coldata$Cluster))


#------------------------------ Output -----------------------------------------
rm(risc_list)
gc()
save.image(file = paste0(output_path,'t_cell_RISC.rdata'))
