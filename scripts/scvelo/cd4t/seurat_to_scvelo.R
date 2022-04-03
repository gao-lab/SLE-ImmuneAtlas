
load('./final/seurat/t_cell/04-') #name 'cd4_filter'

cd4_filter_meta <- merge(cd4_filter@meta.data, Embeddings(cd4_filter,'umap'),by.x = 0, by.y=0,all =T)
write.csv(cd4_filter_meta, file = './scripts/scvelo/cd4t/cd4t_seurat_meta.csv')
rm(cd4_filter_meta)
