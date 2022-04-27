library(tidyverse)
library(tibble)
library(gridExtra)
library(ggplot2)
setwd('/data/sle')
celltypist_anno <- fread('./other_sc_data/extend_data_anno_result.csv',header = T) 
obs <- fread('./other_sc_data/predictions_obs.csv',header = T)

# obs$label <- 'unknown'
obs$label <- celltypist_anno$predicted_labels

table(obs$sample)
table(obs$group)
table(obs$treatment)

obs  %>% 
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    # filter(group != 'treated') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette = get_color(7))+ 
    facet_wrap(~label,scales = "free_y",ncol = 7) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )


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
    stat_compare_means(comparisons = list(c("sle", "hc"),c('sle_flare','hc'),c('sle_treated','hc')),  method = "t.test" )



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
        stat_compare_means(comparisons = list(c("sle", "hc"),c('sle_flare','hc'),c('sle_flare','sle')),  method = "t.test" )
}
do.call(grid.arrange,p)
p
