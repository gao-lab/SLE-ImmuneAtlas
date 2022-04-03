library(tidyverse)

# ----------------------------- data pretreatment ------------------------------
# immcantation and Seurat use different barcode so I need do some transition 
SHM_info <- read.csv('./scripts/immcantaion/bcell_immacantation_SHM.csv')
SHM_info$barcode <- str_split_fixed(SHM_info$sequence_id, pattern = '_',n = 2)[,1] %>% 
    paste0('_',SHM_info$sample)

plasma_filter$barcode <- str_split_fixed(Cells(plasma_filter), pattern = '_',n = 2)[,1]  %>% 
    paste0('_',plasma_filter$orig.ident)
intersect(SHM_info$barcode , plasma_filter$barcode) %>% length()

# 3060/3154   = 96.7%: 96.7% of bcell detect heavy chain in immcantation results


# --------------------------- SHM ratio in cluster -----------------------------
# we only access the heavy chain SHM ratio 
compare_group <-  list(c('HC','untreated'),c('treated','untreated'),c('HC','treated'))
SHM_plasma <- plasma_filter@meta.data %>% left_join(SHM_info, by = c('barcode' = 'barcode')) %>%
    drop_na(mu_freq_seq_r)
SHM_plasma$treatment <- factor(SHM_plasma$treatment,levels = c('HC','untreated','treated')) 

# not use 
SHM_plasma %>% filter(!c_call %in% c('IGHG3','IGHG4','IGHD')) %>%
ggboxplot("c_call", "mu_freq_seq_r", color = "treatment",
          palette = c("#4169E1", "#FF8C00","808000")) + ggtitle("Total mutations(replace)") +
    xlab("Isotype") + ylab("Mutation frequency") +
    theme_bw() + stat_compare_means(aes(group = treatment), label = "p.signif", comparisons = compare_group) 
    # facet_wrap(~c_call,scales = "free")


SHM_plasma %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHG3','IGHD')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", color = "treatment",
               palette = c("#4169E1", "#FF8C00","#808000")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() + 
    facet_grid(subtype ~ c_call)  +
    stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                       comparisons =compare_group,label.y =c(0.1,0.12,0.14)) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

SHM_plasma %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHG3','IGHD')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", color = "treatment",
               palette = c("#4169E1", "#FF8C00","#808000")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") + theme_bw() +  
    facet_wrap(~c_call,ncol = 5)  +
    stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                       comparisons =compare_group,label.y =c(0.13,0.14,0.15)) +
    theme(axis.text.x=element_text(angle=30, hjust=1))


# -------------------------- add SHM to Seurat meta ----------------------------
# plasma_filter@meta.data %>%  left_join(SHM_info, by = c('barcode' = 'barcode')) %>% dim()
tmp <- plasma_filter
tmp_df <- merge(tmp@meta.data,SHM_plasma,by.x ='barcode' ,by.y = 'barcode', sort = F,all.x = T)
tmp$mu_freq_seq_r <- tmp_df$mu_freq_seq_r
tmp$mu_freq_seq_r[is.na(tmp$mu_freq_seq_r)] <- 0
FeaturePlot(tmp, features = 'mu_freq_seq_r', split.by = 'treatment',cols = c('grey','red'),pt.size = 0.1)
