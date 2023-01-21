library(Seurat)
library(future)

pbmc_all$tmp <- 'none'
pbmc_all$tmp[pbmc_all$treatment %in% c('untreated')] <-'basic'
pbmc_all$tmp[pbmc_all$orig.ident %in% c('LL','WYY')] <-'compare'
pbmc_all$tmp %>% table()
plan("multiprocess", workers = 6)
options(future.globals.maxSize= 150 * 1000 * 1024 ^2 )
Idents(pbmc_all) <- 'tmp'
marker_tmp <- FindMarkers(pbmc_all, ident.1 = 'compare', ident.2 = 'basic')

marker_tmp.filter <- marker_tmp %>% filter(p_val_adj < 0.05 & abs(avg_log2FC)  > 0.3)

EnhancedVolcano(marker_tmp.filter,
                lab = rownames(marker_tmp.filter),
                x = 'avg_log2FC',
                y = 'p_val_adj', pCutoff = 10e-32,
                FCcutoff = 0.5) %>% print()
