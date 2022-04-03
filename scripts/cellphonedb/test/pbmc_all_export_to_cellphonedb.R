# Note: not need this scripts anymore

# load the file 
library(Seurat)
library(data.table)
load('./output_file/seurat/all_pbmc/final_concat_pbmc_without_pSS.rdata')

pbmc_final_tmp <- CreateSeuratObject(pbmc_final@assays$RNA@counts, meta.data = pbmc_final@meta.data,
                                     min.cells = 30)
pbmc_final_tmp <- NormalizeData(pbmc_final_tmp, normalization.method = 'RC', scale.factor = 10000)
Idents(pbmc_final_tmp) <- 'subtype'
pbmc_final_downsample <-  subset(pbmc_final_tmp, downsample = 1000)
count_norm <- pbmc_final_downsample@assays$RNA@data
count_norm_dense <- sparse_to_dense(count_norm)


meta_data <- cbind(rownames(pbmc_final_downsample@meta.data), pbmc_final_downsample@meta.data[,'subtype', drop=F])

fwrite(count_norm_dense, './scripts/cellphonedb/cellphonedb_count.txt', sep='\t', quote=F,row.names = F, col.names = T,nThread = 32,buffMB = 256L)
write.table(meta_data, './scripts/cellphonedb/cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
