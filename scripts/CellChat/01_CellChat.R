library(CellChat)

load('./final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')

cellchat <- createCellChat(object = pbmc_all, group.by = "main_type", assay = "RNA")
CellChatDB.use <- CellChatDB
CellChatDB <- CellChatDB.human 
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 16) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

future::plan("sequential")
# ---------------------- Part 2 ------------------------------------------------
# if use truncatedMean parameter, very slow(longer than 1h)
cellchat <- computeCommunProb(cellchat,type ='truncatedMean',trim = 0.1 )
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# future::plan("multiprocess", workers = 8)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat) 

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
heatmap(mat)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   arrow.width = 0.1,arrow.size = 0.05, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# ----------------------- show specific pathway  -------------------------------
pathways.show <- c("MHC-I") # CXCL  BAFF IFN-II MHC-I
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,  arrow.width = 0.1,arrow.size = 0.1,layout = 'chord')
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netAnalysis_contribution(cellchat, signaling = pathways.show)

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object


plotGeneExpression(cellchat, signaling = "MCH-I")


# ----------------------------- Graph theory -----------------------------------
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
gg1 + gg2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

# global (NMF)-----
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_dot(cellchat, pattern = "outgoing")

selectK(cellchat, pattern = "ingoing")

# classify --------
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = 'uwot')
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)



#----------------------------- Save the File ------------------------------------
save(cellchat, file = './final/CellChat/01-all_pbmc_merge_CellChat.rdata')







