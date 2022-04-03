library(tidyverse)

# ----------------------------- data pretreatment ------------------------------
# immcantation and Seurat use different barcode so I need do some transition 
SHM_info <- read.csv('./scripts/immcantaion/bcell_immacantation_SHM.csv')
SHM_info$barcode <- str_split_fixed(SHM_info$sequence_id, pattern = '_',n = 2)[,1] %>% 
    paste0('_',SHM_info$sample)

bcell_filter$barcode <- str_split_fixed(Cells(bcell_filter), pattern = '_',n = 2)[,1]  %>% 
    paste0('_',bcell_filter$orig.ident)
intersect(SHM_info$barcode , bcell_filter$barcode) %>% length()

# 20892/23322 = 90%: 90% of bcell detect heavy chain in immcantation results

# -------------------------- add SHM to Seurat meta ----------------------------
bcell_filter@meta.data %<>%  left_join(SHM_info, by = c('barcode' = 'barcode'))
bcell_filter$mu_freq_seq_r[is.na(bcell_filter$mu_freq_seq_r)] <- 0
FeaturePlot(bcell_filter, features = 'nCount_RNA')

# --------------------------- SHM ratio in cluster -----------------------------
# we only access the heavy chain SHM ratio 
SHM_info %>% group_by(barcode) %>% summarise(max = max(mu_freq_seq_r)) %>% dim()


SHM_cell <- bcell_filter@meta.data %>% left_join(SHM_info, by = c('barcode' = 'barcode')) %>%
    drop_na(mu_freq_seq_r)
ggboxplot(SHM_cell, "subtype", "mu_freq_seq_r", color = "disease",
          palette = c("#4169E1", "#FF8C00")) + ggtitle("Total mutations(replace)") +
    xlab("Isotype") + ylab("Mutation frequency") +
    theme_bw() + stat_compare_means(aes(group = disease), label = "p.signif") +
    facet_wrap(~c_call,scales = "free")

SHM_cell %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
ggboxplot( "subtype", "mu_freq_seq_r", color = "disease",
          palette = c("#4169E1", "#FF8C00")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() + stat_compare_means(aes(group = disease), label = "p.signif") +
    facet_wrap(~c_call,ncol = 4) + theme(axis.text.x=element_text(angle=30, hjust=1))
 
compare_group <-  list(c('HC','untreated'),c('treated','untreated'),c('HC','treated'))
SHM_cell$treatment <- factor(SHM_cell$treatment,levels = c('HC','untreated','treated')) 
SHM_cell %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHD','IGHM')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", color = "treatment",
               palette = c("#4169E1", "#FF8C00","#808000")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() + 
    facet_grid(subtype ~ c_call)  +
    stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                                    comparisons =compare_group,label.y =c(0.1,0.12,0.14)) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

SHM_cell %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHD','IGHM')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", color = "treatment",
               palette = c("#4169E1", "#FF8C00","#808000")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() +  facet_wrap(~c_call)  +
    stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                       comparisons =compare_group,label.y =c(0.13,0.14,0.15)) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

