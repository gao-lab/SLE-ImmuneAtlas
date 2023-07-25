
library(monocle)
library(magrittr)
library(Seurat)

setwd('/data/sle')
output_path <- './output_file/seurat/platelet/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/platelet/platelet_harmony_anno.rdata')

#---------------------------- Monocle Analysis----------------------------------
Idents(platelet_harmony) <- 'subtype'
marker_platelet_harmony_all_cluster <- FindAllMarkers(platelet_harmony,only.pos = F,
                                                      min.pct = 0.2, min.diff.pct = 0.2)
marker_platelet_harmony_all_cluster %<>% filter(p_val_adj <0.05)
sel.gene <- unique(marker_platelet_harmony_all_cluster$gene)

platelet_harmony.cds <- as.CellDataSet(platelet_harmony,assay = '')
platelet_harmony.cds <- estimateSizeFactors(platelet_harmony.cds)
platelet_harmony.cds <- estimateDispersions(platelet_harmony.cds)

platelet_harmony.cds <- monocle::setOrderingFilter(platelet_harmony.cds, sel.gene)
platelet_harmony.cds <- monocle::reduceDimension(platelet_harmony.cds, method = 'DDRTree')
platelet_harmony.cds <- monocle::orderCells(platelet_harmony.cds)


## ordering cells by assigning root nodes
GM_state <- function(cds){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds$subtype)[,"Platelet"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
platelet_harmony.cds <- orderCells(platelet_harmony.cds, 
                                   root_state =  GM_state(platelet_harmony.cds))

# vis
plot_cell_trajectory(platelet_harmony.cds, color_by = "Pseudotime")+ 
  scale_color_gradient2 (low="#FFFFFF", mid="#DCDCDC", high="#FF4500")+ 
  theme(legend.position = "right")

plot_cell_trajectory(platelet_harmony.cds, color_by = "subtype") + facet_wrap(~group)
# plot_cell_trajectory(platelet_harmony.cds, color_by = "group") + facet_wrap(~subtype)
plot_cell_trajectory(platelet_harmony.cds, color_by = "treatment")  + facet_wrap(~treatment)
#---------------------------- CCA before monocle -------------------------------




