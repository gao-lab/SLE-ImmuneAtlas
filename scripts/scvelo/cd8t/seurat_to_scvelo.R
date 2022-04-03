
load('./final/seurat/t_cell/04-') #name 'cd8_filter'

cd8_filter_meta <- merge(cd8_filter@meta.data, Embeddings(cd8_filter,'umap'),by.x = 0, by.y=0,all =T)
write.csv(cd8_filter_meta, file = './scripts/scvelo/cd8t/cd8t_seurat_meta.csv')
rm(cd8_filter_meta)
