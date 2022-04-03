
load('./output_file/seurat/b_cell/bcell_filter_harmony_exclude_BCR_gene.rdata')

b_scvelo_meta <- merge(bcell_filter@meta.data, Embeddings(bcell_filter,'umap'),by.x = 0, by.y=0,all =T)
write.csv(b_scvelo_meta, file = './scripts/scvelo/naive_B/all_B_seurat_meta.csv')

# naive_b from b_cell_slingshot.R
naive_b <- subset(bcell_filter, idents =c('B.IFN-response','B.naive','B.transition'))
naive_b_cells <- CellSelector(DimPlot(naive_b))
naive_b <- subset(naive_b, cells = naive_b_cells , invert =T)
naiveb_scvelo_meta <- merge(naive_b@meta.data, Embeddings(naive_b,'umap'),by.x = 0, by.y=0,all =T)
write.csv(naiveb_scvelo_meta, file = './scripts/scvelo/naive_B/navieB_seurat_meta.csv')
rm(naive_b,naive_b_cells,naiveb_scvelo_meta)
