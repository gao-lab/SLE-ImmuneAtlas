library(ArchR)
set.seed(1)

setwd('/data/sle/atac/working/')
addArchRThreads(threads = 32) 
addArchRGenome("hg38")
source('/data/sle/scripts/function_R/utils.R')

# ------------------------------ Create Arrow Files ----------------------------
input_files <-c('./data/atac_WYF.fragments.tsv.gz','./data/atac_ZS.fragments.tsv.gz')
names(input_files) <- c('wyf','zs')

ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# ---------------------- Remove Doublets and Create Object ---------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "results/",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)


# ------------------------ Reduce Dimension and Cluster ------------------------
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
# proj <- addIterativeLSI(ArchRProj = proj,
#                               useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, 
#                               clusterParams = list(resolution = c(0.2),  sampleCells = 10000, n.start = 10), 
#                               varFeatures = 25000,  dimsToUse = 1:30,force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
proj <- addImputeWeights(proj)


p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")


# --------------------------Use Harmony Remove Batch ---------------------------
## Harmony remove batch effect
proj_harmony<- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI",
                         name = "Harmony",groupBy = "Sample")

proj_harmony <- addClusters(input = proj_harmony, reducedDims = "Harmony",
                          method = "Seurat", name = "Clusters",resolution = 0.5, force = T )

proj_harmony <- addUMAP(ArchRProj = proj_harmony, reducedDims = "IterativeLSI", 
                      name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)
proj_harmony <- addUMAP(ArchRProj = proj_harmony, reducedDims = "Harmony", 
                      name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
proj_harmony <- addImputeWeights(proj_harmony)

### plot the UMAP results
p1_harm <- plotEmbedding(proj_harmony, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2_harm <- plotEmbedding(proj_harmony, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p1_harm, p2_harm, type = "h")
rm(list = c('p1_harm','p2_harm'))

# ------------------------ Vis Gene Scores and Tracks --------------------------
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A", #TCells
  "KLRC1", "KLRD1", "KLRF1", #NK
  
)

p_marker <- plotEmbedding(proj_harmony, colorBy = "GeneScoreMatrix", 
  name = setdiff(marker_list,c('CD16','IGKC','NONE')), 
  embedding = "UMAPHarmony",imputeWeights = getImputeWeights(proj_harmony))

# multi plots
p_marker_all <- lapply(p_marker, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 8),p_marker_all))

# vis the tracks
p_track <- plotBrowserTrack(proj_harmony, groupBy = "Clusters", 
                      geneSymbol = setdiff(marker_list,c('CD16','IGKC','NONE')), 
                      upstream = 50000,downstream = 50000)
grid::grid.newpage()
grid::grid.draw(p$CD3D)

# ---------------------Assigning Clusters with Gene Scores----------------------

# ------------------------ Save the ArchR Project-------------------------------
save(proj, file = './results/Project_save/pbmc_cluster_without_anno.rdata')
save(proj_harmony, file = './results/Project_save/pbmc_cluster_harmony_without_anno.rdata')
