
library(slingshot)
library(dyndimred)

DimPlot(bcell_filter,label = T) + DimPlot(bcell_filter,label = T, group.by = 'subtype')
sds <- slingshot(Embeddings(bcell_filter, "umap"), clusterLabels = bcell_filter$seurat_clusters, 
                 start.clus = 0, end.clus = 7,stretch = 0)

#-------------------------------- Fisrt for naive ------------------------------
# MDS, we use the landmark mds because it is much faster and memory efficient
naive_b <- subset(bcell_filter, idents =c(0,1,2,5))
naive_b <- do_harmony(seu_obj = naive_b,from_begin = T,max.iter = 30)
Idents(naive_b) <- 'disease'
naive_b_HC <- subset(naive_b, idents = 'HC')
naive_b_SLE <- subset(naive_b, idents = 'SLE')
back_run(func = do_harmony,out_name = 'naive_b_HC',job_name = 'naive_b_HC',
         seu_obj = naive_b_HC,from_begin = T,max.iter = 30)
back_run(func = do_harmony,out_name = 'naive_b_SLE',job_name = 'naive_b_SLE',
         seu_obj = naive_b_SLE,from_begin = T,max.iter = 30)
back_run(func = do_harmony,out_name = 'naive_b',job_name = 'naive_b',
         seu_obj = naive_b,from_begin = T,max.iter = 30)

DimPlot(naive_b, group.by = 'subtype') + DimPlot(naive_b, label =T)

sds <- slingshot(Embeddings(naive_b, "umap"), clusterLabels = naive_b$subtype, 
                 start.clus = 'B.transition',end.clus = c('B.IFN-response','B.naive'), stretch = 2 )
cell_colors <- cell_pal(naive_b$subtype, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(naive_b$seurat_clusters, hue_pal())
plot(Embeddings(naive_b, "umap"), col = cell_colors, pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd = 2, type = 'lineages', col = 'black')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sds@assays@data$pseudotime, breaks=100)]

plot(Embeddings(naive_b, "umap"), col = cell_colors, pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd=2, col='black')


# -----------------------------Failure try -------------------------------------
# NOTE follow https://github.com/NBISweden/single-cell_sib_scilifelab/blob/master/session-trajectories/3_slingshot.md
# naive_b <- NormalizeData(object = naive_b)
# naive_b <- FindVariableFeatures(object = naive_b,nfeatures = 2000)
# scale_all =F
# if(scale_all ){
#     all.genes <- rownames(naive_b)
#     naive_b <- ScaleData(object = naive_b, features = all.genes)
# }else{naive_b <- ScaleData(object = naive_b)}
# lmds <- dyndimred::dimred_landmark_mds(naive_b@reductions$harmony@cell.embeddings)
# colnames(lmds) <- paste0("lmds_", seq_len(ncol(lmds)))
# lmds_object <- CreateDimReducObject(lmds, key = "lmds_", assay = "scale_data")
# naive_b@reductions$lmds <- lmds_object
# DimPlot(naive_b, reduction = "lmds",pt.size = 0.5, label = TRUE, repel = TRUE,group.by = 'subtype')

# dimred <- naive_b@reductions$umap@cell.embeddings
# clustering <- naive_b@meta.data$subtype
# lineages <- slingshot::getLineages(dimred, clustering,start.clus= 'B.naive')

# plot.new()
# plot(dimred, col = RColorBrewer::brewer.pal(9,"Set1")[clustering], asp = 1, pch = 16)
# lines(lineages, lwd = 3, col = 'black')