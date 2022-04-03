
library(shiny)
library(Seurat)
setwd('/data/sle/vdj/VDJview/')

# Inside the docker we may face port problem, so set the launch.broswer False
# shiny::runApp('./VDJView', launch.browser = F)
# open the ui.R and click runapp in rstudio

# GW h5ad to mtx dir
library(DropletUtils)
gw_seu <- Read10X_h5('/data/sle/data/10x_rna/GW_rna.h5')
gw_seu <- CreateSeuratObject(gw_seu)
write10xCounts(x = gw_seu@assays$RNA@counts, path = './test_data/GW_10X_mtx',version = '3')
rm(gw_seu)
gw_bcr <- read.csv('./test_data/GW_10X_bcr/filtered_contig_annotations.csv')

gw_b_seu <- subset(gw_seu, cells = intersect(Cells(gw_seu), gw_bcr$barcode))
gw_b_seu <- CreateSeuratObject(gw_b_seu@assays$RNA@counts, min.cells = 10)

b_cell_df <- as.data.frame(gw_b_seu@assays$RNA@counts,keep.rownames = T )
fwrite(b_cell_df,file = './test_data/GW_10X_mtx/genes.fpkm_table',sep = '\t',row.names =T,col.names = T )

ls()


