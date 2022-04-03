
load('./final/seurat/mono_dc//03-mono_dc_anno_filter_harm.rdata') #name 'mono_dc_filter'

mono_dc_filter_meta <- merge(mono_dc_filter@meta.data, Embeddings(mono_dc_filter,'umap'),by.x = 0, by.y=0,all =T)
write.csv(mono_dc_filter_meta, file = './scripts/scvelo/mono//mono_dc_seurat_meta.csv')
rm(mono_dc_filter_meta)
