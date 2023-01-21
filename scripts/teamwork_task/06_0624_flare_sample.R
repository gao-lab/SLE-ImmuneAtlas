# how do we distinguish flare sample and real treatment naive

# plasmablast
(table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>%
    filter(! treatment == 'HC') %>% 
    unique() %>% ggboxplot(x ='treatment' , y = 'Freq', color = 'treatment') + 
    stat_compare_means(paired = F, method = 't.test',label.x = 1.5,label = "p.signif")


t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>% 
    filter(subtype == 'T.prolife_T' & ! treatment =='HC' & ! orig.ident %in% c('LL','LL2') ) %>% 
    ggboxplot(x = 'treatment', y ='Ratio') + stat_compare_means(paired = F, method = 'wilcox.test',label.x = 1.5)

# all subtype
tmp <- pbmc_all@meta.data %>% group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>%
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pbmc_all@meta.data[,c(1,4,5,6)]) %>% filter(treatment == 'untreated') %>% unique() %>%
    mutate(flare = 'F', flare = if_else(grepl('WYY|LL',orig.ident), 'T', flare)) 

ggdotplot(tmp, x = 'flare', y = 'Ratio', color = 'flare', label = 'orig.ident') + facet_wrap(~subtype,scales = 'free')          
