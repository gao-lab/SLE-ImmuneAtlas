
#-----------all maintype change in the pbmc------------------

# By subtype
pbmc_all@meta.data  %>% 
    group_by(orig.ident,main_type) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07" ,"#00AFBB"))+ 
    facet_wrap(~main_type,scales = "free",ncol = 4) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )

# add  cd4 and cd8 label
pbmc_all$cd4_8 <- 'not'
pbmc_all$cd4_8[which( grepl('CD4',pbmc_all$subtype))] <- 'CD4'
pbmc_all$cd4_8[which( grepl('CD8',pbmc_all$subtype))] <- 'CD8'

pbmc_all@meta.data  %>% 
    group_by(orig.ident,cd4_8) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07" ,"#00AFBB"))+ 
    facet_wrap(~cd4_8,scales = "free",ncol = 4) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )
