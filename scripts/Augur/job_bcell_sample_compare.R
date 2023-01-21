library(Seurat)
library(Augur)
library(tidyverse)
library(ggpubr)
library(parallel)

load('/data/sle/scripts/Augur/bcell.rdata')

pairs <- names(bcell_filter$pair %>% table())[-3]
Idents(bcell_filter) <- 'group'
bcell_hc <- subset(bcell_filter, idents = 'HC')
Idents(bcell_filter) <- 'pair'
augur_list <- mclapply(pairs,FUN = function(pair){
    tmp <- subset(bcell_filter, idents = pair)
    Idents(tmp) <- 'treatment'
    before <- subset(tmp, idents = 'untreated')
    after <- subset(tmp ,idents = 'treated')
    before <- merge(before, bcell_hc)
    after <- merge(after, bcell_hc)
    mat_before <- before@assays$RNA@counts %>% as.matrix()
    mat_after <- after@assays$RNA@counts %>% as.matrix()
    before <- calculate_auc(input = mat_before,meta = before@meta.data,
                            label_col = 'treatment',
                            cell_type_col = 'subtype', n_threads=8)
    after <- calculate_auc(input = mat_after,meta = after@meta.data,
                           label_col = 'treatment',
                           cell_type_col = 'subtype', n_threads=8)
    augur_df <- data.frame(before$AUC) %>%  mutate(group = 'HC VS untreated') %>% 
        rbind(after$AUC %>%  mutate(group = 'HC VS treated'))
    augur_df$group <- factor(augur_df$group,levels = c('HC VS untreated','HC VS treated'))
    return(augur_df)
},mc.cores = 32)

save(augur_list,'./scripts/Augur/augur_list.bcell.rdata')

