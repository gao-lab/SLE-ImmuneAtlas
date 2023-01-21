library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)

# -------------------------- follow the tutouial -------------------------------
# copy the tutotial from https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/monocle2.html
# monocle2 will run very slow in a big single cell dataset .So I downsample the 
# dataset , and ensure every subtype has same cell number

Idents(bcell_filter) <- 'group'
bcell_sle <- subset(bcell_filter, idents = 'SLE')
bcell_hc <- subset(bcell_filter, idents = 'HC')
Idents(bcell_hc) <- 'subtype'
Idents(bcell_filter) <- 'subtype'
cds <- subset(bcell_filter, downsample = 200)

cds <- as.CellDataSet(cds)
cds <- estimateSizeFactors(cds)
# cds <- estimateDispersions(cds) # long time about 10 min for 25k cell
back_run(estimateSizeFactors,out_name = 'cds',job_name = 'estimateSizeFactors',cds)

cds <- detectGenes(cds,min_expr = 0.1)

# deg <- marker_all_bcell_filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10)
# deg <- marker_all_bcell_filter %>% filter(avg_log2FC > 0.5 )
deg <- unique(marker_all_bcell_filter$gene)
cds <- monocle::setOrderingFilter(cds, deg)

## dimension reduciton
# cds <- monocle::reduceDimension(cds, method = 'DDRTree')
back_run(reduceDimension, out_name = 'cds',job_name = 'DDRTree',
         cds, method = 'DDRTree')

## ordering cells(use the follow version)
# cds <- monocle::orderCells(cds) # very slow , use backend mode
back_run(orderCells,out_name = 'cds',job_name = 'monocle::orderCells',cds,num_paths =1)
plot_cell_trajectory(cds, color_by = "subtype")
plot_cell_trajectory(cds, color_by = "subtype") + facet_wrap(~subtype)

## ordering cells by assigning root nodes
GM_state <- function(cds){
    if (length(unique(cds$State)) > 1){
        T0_counts <- table(cds$State, cds$subtype)[,"B.naive"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
cds <- monocle::orderCells(cds, root_state =  GM_state(cds))
monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")

#-----------------------------  vis the genes ----------------------------------

my_genes <- c("IGHD", "IGHG1", "VPREB3", "FCRL5", "ISG15")
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "subtype")


#-----------------------------  branch analysis ---------------------------------
BEAM_res <- BEAM(cds, branch_point = 1, cores = 24)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res <- BEAM_res %>% filter(!str_detect(gene_short_name, "^RP[LS]")) %>% 
    filter(!str_detect(gene_short_name, "^MT")) %>% 
    filter(!str_detect(gene_short_name, "^IG[KL]")) %>% 
    filter(!str_detect(gene_short_name, "^IGHV"))
plot_genes_branched_heatmap(cds[c(row.names(subset(BEAM_res,
                                                  qval < 1e-25)),'CD27','IGHD','IGHM'),],
                            branch_labels = c("B.mem.CD27-", "B.memory"),
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T)
plot_genes_branched_pseudotime(cds[c('CD27','CD79A','IGHG2','IGHA1','IGHD','CD38','FCRL5','CXCR3','CXCR4'),],
                               branch_point = 1,branch_labels = c("B.mem.CD27-", "B.memory"),
                               color_by = "subtype",
                               ncol = 3)
plot_genes_branched_pseudotime(cds[c('CXCR5','CD27','ISG15','TLR7','TLR9')],
                               branch_point = 1,branch_labels = c("B.mem.CD27-", "B.memory"),
                               color_by = "subtype")
plot_genes_branched_pseudotime(cds[c('CXCR5','CD27','ISG15','TLR7','TLR9')],
                               branch_point = 1,branch_labels = c("B.mem.CD27-", "B.memory"),
                               color_by = "group")
