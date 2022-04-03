library(glue, lib.loc = "/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)

setwd('/data/sle/scripts/NicheNet/')
source('../function_R/utils.R')
load('../../final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')

Idents(pbmc_all) <- 'subtype'
pbmc_all_downsample <- subset(pbmc_all, downsample =800)
save(pbmc_all_downsample, file = 'pbmc_all_downsample.rdata')

ligand_target_matrix <- readRDS('./source/ligand_target_matrix.rds')
lr_network <- readRDS('./source/lr_network.rds')
weighted_networks <- readRDS('./source/weighted_networks.rds')

sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

# ----------------------------- Define a receiver(target) ----------------------
## receiver
Idents(pbmc_all_downsample) <- 'main_type'
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = pbmc_all_downsample, 
  receiver = "Plasma", 
  condition_colname = "group", condition_oi = "SLE", condition_reference = "HC", 
  # sender = c("NK","cDC", "T.prolife", "pDC", "Plasma", "Mono.CD16"), 
  sender = 'all',
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, organism = "human",lfc_cutoff = 0.2,cutoff_visualization = 0.1)
nichenet_output$top_ligands
nichenet_output$ligand_activities
nichenet_output$ligand_differential_expression_heatmap

nichenet_output2 = nichenet_seuratobj_aggregate(
  seurat_obj = pbmc_all_downsample, 
  receiver = "Bcell", 
  condition_colname = "group", condition_oi = "SLE", condition_reference = "HC", 
  # sender = c("NK","cDC", "T.prolife", "pDC", "Plasma", "Mono.CD16"), 
  sender = 'undefined',
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, organism = "human",lfc_cutoff = 0.2)


# ------------------------------ build network ---------------------------------
ligands_all = "CD40LG" # this can be a list of multiple ligands if required
targets_all = c("CCR2","CD38")


active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix,
                                                     ligands_all = ligands_all,
                                                     targets_all = targets_all, 
                                                     weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>%
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(
  signaling_graph_list = active_signaling_network_min_max, 
  ligands_all = ligands_all, 
  targets_all = targets_all, 
  sig_color = "indianred", 
  gr_color = "steelblue")
DiagrammeR::render_graph(graph_min_max, layout = "tree")
