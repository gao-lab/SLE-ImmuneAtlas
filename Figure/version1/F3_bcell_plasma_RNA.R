################################################################################
#
# Figure 3B: Marker of plasma
#
###################################################################3############
DotPlot2(plasma_filter, marker_list = c('CD38','IGHA1','IGHG1','MKI67'), group.by = 'subtype')


################################################################################
#
# Figure 3C: Plasmablast subtype ratio and after treatment
#
###################################################################3############
plasma_filter@meta.data  %>% group_by(orig.ident,subtype) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(plasma_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    filter(subtype %in% c("plasmablast")) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(subtype~.,scales = "free", ncol = 1)+ theme_cowplot() +
    stat_compare_means(label.x = 1.4,label = "p.signif") + 
    ylab('Ratio of Plamsa') + xlab('')

tmp1<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='treated') %>%
    filter(!orig.ident == 'XYY2')
tmp2<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='untreated')
data.frame(sample = tmp1$orig.ident, before=tmp2$Freq, after=tmp1$Freq ) %>%
    mutate(subtype = 'plasmablast') %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
              fill  = "condition", line.color = "gray", line.size = 0.4,
              palette =c('#DA9494','#9FB1D4')) +  
    stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.5, label = 'p.signif')+
    ylab('Ratio of plasma') + facet_wrap(subtype~.,scales = "free", ncol = 1)+ theme_cowplot()
rm(tmp1,tmp2)

################################################################################
#
# Extend Data Fig3a: IFN subtype ratio 
#
###################################################################3############
bcell_filter@meta.data  %>% group_by(orig.ident,subtype) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
                                             "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
    filter(subtype %in% c("B.IFN-response")) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(subtype~.,scales = "free", ncol = 1)+ theme_cowplot() +
    stat_compare_means(label.x = 1.5,label = "p.signif") + xlab('') +
    ylab('Ratio of B cell')

################################################################################
#
# Extend Data Fig33: DN subtype ratio 
#
###################################################################3############
bcell_filter@meta.data  %>% group_by(orig.ident,subtype) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
                                             "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
    filter(subtype %in% c("B.mem.CD27-")) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(subtype~.,scales = "free", ncol = 1)+ theme_cowplot() +
    stat_compare_means(label.x = 1.5,label = "p.format")+ xlab('') +
    ylab('Ratio of B cell')

