# test B cell conga ratio
library(gridExtra)

tmp <- read.csv('./scripts/conga/BCR/all/cong_bcell_adata_meta.csv')
tmp$conga <- 'no_significant'
tmp$conga[which(tmp$conga_scores < 0.001)] <- 'significant'

table(tmp$subtype)
hist(tmp$conga)
hist((tmp %>% filter(conga == 'no_significant'))$conga_scores,breaks = 100)
hist((tmp %>% filter(conga == 'significant'))$conga_scores,breaks = 100)

hist((tmp %>% filter(conga == 'no_significant'))$conga_fdr_value,breaks = 100)
hist((tmp %>% filter(conga == 'significant'))$conga_fdr_value,breaks = 100)

# only in B cell
tmp  %>% filter(!subtype %in% c('plasma','plasma.IgA','plasma.IgG','plasmablast')) %>%
    group_by(orig.ident,conga) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_bcr@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                      palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
    facet_wrap(~conga,scales = "free")+ stat_compare_means()

# only in plasma
tmp %>% filter(subtype %in% c('plasma','plasma.IgA','plasma.IgG','plasmablast')) %>%
    group_by(orig.ident,conga) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(plasma_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                      palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
    facet_wrap(~conga,scales = "free")+ stat_compare_means()

# in subtype
p <- list()
for (b_subtype in unique(tmp$subtype)) {
    print(b_subtype)
    p[[b_subtype]] <- tmp %>% filter(subtype == b_subtype) %>%
        group_by(orig.ident,conga) %>% 
        summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
        mutate(Ratio = sub_num/sample_num*100) %>%
        left_join(plasma_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
        filter(treatment != 'treated') %>%  filter(conga != 'significant') %>% 
        ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                          palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ ggtitle(b_subtype)+ 
        stat_compare_means(method = 't.test')
}
do.call(grid.arrange,p)

