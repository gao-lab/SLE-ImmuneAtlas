setwd('/data/sle/')
source('./scripts/function_R/utils.R')

library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)


ligand_target_matrix = readRDS('./scripts/NicheNet/source/ligand_target_matrix.rds')
lr_network = readRDS('./scripts/NicheNet/source/lr_network.rds')
weighted_networks = readRDS('./scripts/NicheNet/source/weighted_networks.rds')

load('./output_file/seurat/all_pbmc/modify_subtype_final_concat_pbmc_without_pSS.rdata')

Idents(pbmc_final) <- 'subtype'
pbmc_final_downsample <- subset(pbmc_final, downsample =1000)
save(pbmc_final_downsample, file = './output_file/seurat/all_pbmc/downsample_final_concat_pbmc_without_pSS.rdata')
load('./output_file/seurat/all_pbmc/downsample_final_concat_pbmc_without_pSS.rdata')

all_cell_type <- names(table(pbmc_final_downsample$subtype))
Idents(pbmc_final_downsample) <- 'subtype'
nichenet_output.pDC = nichenet_seuratobj_aggregate(
    seurat_obj = pbmc_final_downsample,
    receiver = "pDC",sender = all_cell_type,
    condition_colname = "disease", condition_oi = "SLE", condition_reference = "HC",
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, organism = "human",lfc_cutoff = 0.2)
# back_run(func = nichenet_seuratobj_aggregate, out_name = 'nichenet_output.pDC',job_name = 'nichenet_output.pDC',
#          seurat_obj = pbmc_final_downsample, 
#          receiver = "B.IFN-response",sender = names(table(pbmc_final_downsample$subtype)), 
#          condition_colname = "disease", condition_oi = "SLE", condition_reference = "HC", 
#          ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
#          weighted_networks = weighted_networks, organism = "human")
nichenet_output.CD8_Temra_active = nichenet_seuratobj_aggregate(
    seurat_obj = pbmc_final_downsample,
    receiver = "CD8_Temra_active",sender = all_cell_type,
    condition_colname = "disease", condition_oi = "SLE", condition_reference = "HC",
    ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
    weighted_networks = weighted_networks, organism = "human",lfc_cutoff = 0.2)

nichenet_output.CD8_Temra_active$ligand_expression_dotplot + scale_y_discrete(guide = guide_axis(angle = 30))
nichenet_output.CD8_Temra_active$ligand_target_heatmap
nichenet_output.CD8_Temra_active$ligand_differential_expression_heatmap
nichenet_output.CD8_Temra_active$ligand_activity_target_heatmap # overall about 3 pics before
nichenet_output.CD8_Temra_active$ligand_receptor_heatmap
nichenet_output.CD8_Temra_active$ligand_receptor_heatmap_bonafide

Idents(pbmc_final_downsample) <- 'compare_meta'
DotPlot(pbmc_final_downsample %>% subset(idents = 'B cells'), features = nichenet_output.pDC$top_targets %>% rev(), 
        group.by = 'subtype', split.by = "disease",cols = "RdYlBu") + RotatedAxis() & coord_flip() 

DotPlot(pbmc_final_downsample %>% subset(idents = c('B cells','CD4 T cells','CD8 T cells')), 
        features = 'LTB', split.by = 'disease', group.by = 'compare_meta', cols = 'RdYlBu')
DotPlot(pbmc_final_downsample %>% subset(idents = c('CD8 T cells','CD4 T cells')), 
        features = c('LTB','LTA'), split.by = 'disease', group.by = 'compare_meta', cols = 'RdYlBu')
DotPlot(pbmc_final_downsample, features= c('LTB','LTA'),cols = 'RdYlBu',split.by = 'disease')

