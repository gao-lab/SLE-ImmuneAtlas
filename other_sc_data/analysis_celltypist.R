library(tidyverse)
library(tibble)
library(gridExtra)
library(ggplot2)
setwd('/data/sle')
celltypist_anno <- fread('./other_sc_data/extend_data_anno_result.csv',header = T) 
obs <- fread('./other_sc_data/predictions_obs.csv',header = T)

# obs$label <- 'unknown'
obs$label <- celltypist_anno$predicted_labels

# table(obs$sample)
# table(obs$group)
# table(obs$treatment)

obs$group <- factor(obs$group, levels = c('sle_treated','sle_flare','sle','hc','sle_child','hc_child','IFNbeta_stim'))


obs  %>% 
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    filter(!group %in% c('IFNbeta_stim','sle_treated') ) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette = get_color(7))+ 
    facet_wrap(~label,scales = "free_y",ncol = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    stat_compare_means(label = "p.signif" ,comparisons = list(c("sle", "sle_flare"),c('hc_child','sle_child'),c('sle','hc'),c('sle_flare','hc')),
                       method = 't.test')


stat_obs <- obs  %>% 
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') )
table(stat_obs$group)

stat_obs <- obs  %>%  select(sample, group,dataset) %>% group_by(group) %>% distinct() 
table(stat_obs$group ,stat_obs$dataset)


# try to merge some group 
obs$group_main <- obs$group 
obs$group_main[obs$group == 'hc_child'] <- 'hc'
obs$group_main[obs$group == 'sle_child'] <- 'sle'
obs$group_main[obs$group == 'sle_treated'] <- 'sle'

obs  %>% 
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28,55)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    ggpubr::ggboxplot(x='group_main',y='Ratio', fill = 'group_main',
                      palette = get_color(7))+ 
    facet_wrap(~label,scales = "free_y",ncol = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + 
    stat_compare_means(comparisons = list(c("sle", "hc"),c('sle_flare','hc'),c('sle_treated','hc')) )



#########################################################
# caculate ratio within the  main cell type 
#########################################################
unique(obs$label)
t_subtype <- unique(obs$label)[1:18]
mono_subtype <- unique(obs$label)[19:28]
b_subtype <- unique(obs$label)[29:35]
plasma_subtype <- unique(obs$label)[38:42]

# plot 
label_list <- list(t_subtype,mono_subtype,b_subtype,plasma_subtype)
p <- list()
for(i in 1:4){
    p[[i]] <- obs  %>% filter(label %in% label_list[[i]]) %>%
        group_by(sample,label) %>% summarise(sub_num = n()) %>% 
        mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
        left_join(obs[,c(5,28,55)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
        ggpubr::ggboxplot(x='group_main',y='Ratio', fill = 'group_main',
                          palette = get_color(7))+ 
        facet_wrap(~label,scales = "free_y",ncol = 7) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + 
        stat_compare_means(comparisons = list(c("sle", "hc"),c('sle_flare','hc'),c('sle_flare','sle')),label = "p.signif")
}
do.call(grid.arrange,p)
p


pbmc_all@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pbmc_all@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    ggpubr::ggboxplot(x='treatment',y='Ratio', fill = 'treatment',
                      palette =c("#B3D492","#DA9494","#9FB1D4"))+ 
    facet_wrap(~subtype,scales = "free",nrow = 4) + 
    stat_compare_means(comparisons = my_comparisons, hide.ns = F,label = "p.signif") 


################################################################################
#
# pDC cell number decrease in SLE
#
################################################################################
obs  %>% 
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28,55)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    filter(label == 'pDC') %>%
    ggpubr::ggboxplot(x='group_main',y='Ratio', fill = 'group_main',
                      palette = get_color(7))+ 
    facet_wrap(~label,scales = "free_y",ncol = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + 
    stat_compare_means(comparisons = list(c("sle", "hc"),c('sle_flare','hc') ))

