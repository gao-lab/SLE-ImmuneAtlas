# try to use cca to integrate 6 smart-seq data set
# 2 of these data set do not contain TLB cells
library(Seurat)
library(future)

source('/data/teamwork/sle/function_R/markers.R')

#--------------------------------------------------------------------------------------
# only for B cell
#--------------------------------------------------------------------------------------

setwd('/data/teamwork/sle/data/xiehe/tlb_scan/batch_correct/02_6_smartseq/')
b_cell_GSE146771 <- readRDS('../../smartseq2/GSE146771/b_cell_GSE146771.Rds')
b_cell_GSE115978 <- readRDS('../../smartseq2/GSE115978/b_cell_GSE115978.Rds')
b_cell_GSE103322 <- readRDS('../../smartseq2/GSE103322/b_cell_GSE103322.Rds')
b_cell_GSE120575 <- readRDS('../../smartseq2/GSE120575/b_cell_GSE120575.Rds')
b_cell_GSE72056  <- readRDS('../../smartseq2/GSE72056/b_cell_GSE72056.Rds')
b_cell_GSE140228 <- readRDS('../../smartseq2/GSE140228/b_cell_GSE140228.Rds')

b_cell_GSE146771$dataset <- 'GSE146771'
b_cell_GSE115978$dataset <- 'GSE115978'
b_cell_GSE103322$dataset <- 'GSE103322'
b_cell_GSE120575$dataset <- 'GSE120575'
b_cell_GSE72056$dataset <- 'GSE72056'
b_cell_GSE140228$dataset <- 'GSE140228'

seu_list <- c(b_cell_GSE146771,b_cell_GSE115978,b_cell_GSE103322,
              b_cell_GSE120575,b_cell_GSE72056,b_cell_GSE140228)

# cca align cells
plan("multiprocess", workers = 16); options(future.globals.maxSize= 150 * 1000 * 1024 ^2 ) 
b_cell_anchors <- FindIntegrationAnchors(seu_list)
plan('sequential')
b_integrated <- IntegrateData(anchorset = b_cell_anchors)
DefaultAssay(b_integrated) <- "integrated"

b_integrated <- ScaleData(b_integrated, verbose = FALSE)
b_integrated <- RunPCA(b_integrated, npcs = 30, verbose = FALSE)
b_integrated <- RunUMAP(b_integrated, reduction = "pca", dims = 1:30,verbose = F)
b_integrated <- FindNeighbors(b_integrated, dims = 1:30,verbose = F)
b_integrated <- FindClusters(b_integrated, resolution = 0.8,verbose = F)

DimPlot(b_integrated, label = T)
DimPlot(b_integrated, group.by = 'dataset')

# downstream analysis
DefaultAssay(b_integrated) <- "RNA"
FeaturePlot(b_integrated, features = c('CD79A','CD3D','MS4A1','XBP1'))
VlnPlot(b_integrated,features = c('CD79A','CD3D','MS4A1','CD79B','PAX5','XBP1','MZB1','SLAMF7','SDC1'), stack = T)
VlnPlot(b_integrated,features = plasma_marker, stack = T)
VlnPlot(b_integrated,features = feats, stack = T)


# find markers 
plan("multiprocess", workers = 16); options(future.globals.maxSize= 150 * 1000 * 1024 ^2 )
marker_binte_c9 <- FindMarkers(b_integrated,ident.1 = 9, min.pct = 0.25, only.pos = T)
marker_binte_c11 <- FindMarkers(b_integrated,ident.1 = 11, min.pct = 0.25, only.pos = T)
plan('sequential')
head(marker_binte_c9[order(marker_binte_c9$avg_log2FC),],20)
head(marker_binte_c11[order(marker_binte_c11$avg_log2FC),],20)

# cell type annotation 


table(b_integrated$integrated_snn_res.0.8, b_integrated$dataset)
#--------------------------------------------------------------------------------------
# tlb analysis
#--------------------------------------------------------------------------------------
tlb <- subset(b_integrated, idents = 7)


#--------------------------------------------------------------------------------------
# single R confirm
#--------------------------------------------------------------------------------------
ref <- BlueprintEncodeData()
pred.b_integrated<- SingleR(b_integrated@assays$RNA@counts,ref = ref, labels = ref$label.fine,
                            BPPARAM=MulticoreParam(8))
plotScoreHeatmap(pred.b_integrated)
tab.b_integrated <- table(Assigned=pred.b_integrated$pruned.labels, Cluster=b_integrated$integrated_snn_res.0.8)
pheatmap(log2(tab.b_integrated +10), color=colorRampPalette(c("white", "blue"))(101))
