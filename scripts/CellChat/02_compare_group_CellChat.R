library(CellChat)
library(patchwork)
library(Seurat)
setwd('/data/sle')
mc <- getOption("mc.cores", 36)
CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

# -------------------------- Prepare the Data ----------------------------------
data.dir <- './final/CellChat/comparison'
dir.create(data.dir)
setwd(data.dir)
load('../../seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')
group_list <- SplitObject(pbmc_all, split.by = 'treatment')
rm(pbmc_all);gc()
group_list <- lapply(group_list, function(x) 
  createCellChat(object = x, group.by = "main_type", assay = "RNA"))
# This step need about 30min and 100G MEM
options(future.globals.maxSize= 8*1024^3)
group_list <- mclapply(group_list, process_cellchat, mc.cores = mc)
# not need anymore because we add this step to process_cellchat()
group_list <- mclapply(group_list, function(x)
  netAnalysis_computeCentrality(x, slot.name = "netP"), mc.cores = mc)

cellchat <- mergeCellChat(group_list, add.names = names(group_list))
rm(group_list);gc()

process_cellchat <- function(cellchat_obj, workers = 8){
  cellchat_obj@DB <- CellChatDB.use
  cellchat_obj <- subsetData(cellchat_obj)
  future::plan("multiprocess", workers = workers) 
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  # cellchat_obj <- projectData(cellchat_obj, PPI.human)
  cellchat_obj <- computeCommunProb(cellchat_obj)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
}

# ------------------------------------ Plot ------------------------------------
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight",size.text =10)
gg1 + gg2
gg2 + scale_fill_manual(values = c('#DA9494','#9FB1D4','#B4D493')) 

# circle plot 
par(mfrow = c(1,3), xpd=TRUE)
netVisual_diffInteraction(cellchat,comparison = c(1,3), weight.scale = T,arrow.size = 0.4,
                          title.name = 'Untreated compare HC')
netVisual_diffInteraction(cellchat,comparison = c(2,1), weight.scale = T,arrow.size = 0.4,
                          title.name = 'Treated compare Untreated')
SnetVisual_diffInteraction(cellchat,comparison = c(1,2), weight.scale = T,arrow.size = 0.4)
netVisual_diffInteraction(cellchat,comparison = c(2,3), weight.scale = T,arrow.size = 0.4)

# heatmap
gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 3))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1, 3))
gg1 + gg2

# compare signal pathway
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasma", signaling.exclude = "MIF")

# see variable pathway
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = T,comparison = c(1, 3),
        color.use = c('#DA9494','#B4D493'), font.size = 10)
netAnalysis_contribution(group_list[[1]], signaling = 'CCL',return.data = T, font.size = 10)
tmp <- netAnalysis_contribution(group_list[[1]], signaling = 'CCL',return.data = T, font.size = 10)
ggbarplot(tmp$LR.contribution,x = 'name',y = 'contribution', fill = 'name', color = 'name', palette = c( "#5D9DCC","#E7525A"))
netAnalysis_contribution(group_list[[2]], signaling = 'CCL',return.data = T, font.size = 10)

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5:11),  comparison = c(1, 3), angle.x = 45)

features.name = 'untreated';pos.dataset = 'untreated'
plan('sequential')
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "untreated",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "untreated",ligand.logFC = -0.1, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), 
                        comparison = c(1, 3),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in Untreated"))


################################################################################
#
# Show specific signal pathway
#
################################################################################
pathways.show <- c("CCL") 
# debugonce(getMaxWeight)
weight.max <- getMaxWeight(group_list[1:2], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(group_list[1:2])) {
  netVisual_aggregate(group_list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(group_list)[i]))
}




################################################################################
#
# Show general signal pathway
#
################################################################################
library(ComplexHeatmap)
pathway.union <- Reduce(union,list(group_list[[1]]@netP$pathways, group_list[[2]]@netP$pathways, group_list[[3]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(group_list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(group_list)[1], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(group_list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(group_list)[2], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(group_list[[3]], pattern = "outgoing", signaling = pathway.union, title = names(group_list)[3], width = 5, height = 6)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(group_list[[1]], pattern = "incoming", signaling = pathway.union, title = names(group_list)[1], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(group_list[[2]], pattern = "incoming", signaling = pathway.union, title = names(group_list)[2], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(group_list[[3]], pattern = "incoming", signaling = pathway.union, title = names(group_list)[3], width = 5, height = 6)
draw(ht1 + ht2 + ht3, ht_gap = unit(0.1, "cm"))
