

load('output_file/seurat/all_pbmc/modify_subtype_final_concat_pbmc_without_pSS.rdata')

marker_gene <- 'CXCR4'
# my_comparisons <- list( c("HC_HC", "SLE_pah_treated"), c("HC_HC", "SLE_pah_untreated"),
#                         c("SLE_pah_treated","SLE_pah_untreated"))
pbmc_final$compare_tmp <- paste0(pbmc_final$group,'-', pbmc_final$treatment)
Idents(pbmc_final) <- 'compare_tmp'
tmp <- subset(pbmc_final, idents = c('HC-HC','SLE_pah-untreated','SLE_pah-treated'), invert = F) 
# only B CD4 CD8 mono macrob NK DC
Idents(tmp) <- 'compare_meta'
tmp <- subset(tmp, idents = c('B cells','CD14.Mono','CD16.Mono','CD4 T cells',
                                     'CD8 T cells','Macrophage','NK cells','pDC','mDC'), invert = F) 

# tmp %>% VlnPlot( group.by = 'compare_meta', split.by = 'compare_tmp', features = marker_gene, pt.size = 0) 
    # stat_compare_means(aes(group = tmp$compare_group)) # + # Add pairwise comparisons p-value
    # stat_compare_means(label.y = 50)     # Add global p-value
    


tmp$compare_meta[which(tmp$compare_meta %in% c('CD14.Mono','CD16.Mono'))] <- 'Monocyte'
tmp$compare_meta[which(tmp$compare_meta %in% c('pDC','mDC'))] <- 'DC cells'
# arrange the order
tmp$compare_tmp <- factor(x =tmp$compare_tmp , levels = c('HC-HC','SLE_pah-untreated','SLE_pah-treated'))
tmp %>% VlnPlot( group.by = 'compare_meta', split.by = 'compare_tmp', 
                 features = marker_gene, pt.size = 0.05) + xlab('')

Idents(pbmc_final) <- 'compare_meta'
tmp2 <- subset(pbmc_final, idents = 'B cells')
VlnPlot(tmp2, features = 'CXCR4', group.by = 'subtype') + xlab('')

tmp2$compare_tmp <- paste0(tmp2$group,'-', tmp2$treatment)
Idents(tmp2) <- 'compare_tmp'
tmp2 <- subset(tmp2, idents = c('HC-HC','SLE_pah-untreated','SLE_pah-treated'), invert = F)
tmp2$compare_tmp <- factor(x =tmp2$compare_tmp , levels = c('HC-HC','SLE_pah-untreated','SLE_pah-treated'))
VlnPlot(tmp2, features = 'CXCR4', group.by = 'subtype', split.by = 'compare_tmp') + xlab('')
# rm(pbmc_final,tmp)

