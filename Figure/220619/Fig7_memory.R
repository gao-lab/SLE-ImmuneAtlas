library(Seurat)

################################################################################
#
# memory B and CD8 cell proportion increase after treatment
#
################################################################################
# memory B
tmp1<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY')
tmp2<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY')
b_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                           before=tmp2$Ratio, after=tmp1$Ratio )
p1 <- b_pair_df %>% filter(subtype == 'B.mem') %>%
    ggpaired(b_pair_df, cond1 = 'before', cond2 = 'after',legend = 'right',point.size = 3,
             color = "condition", line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(label = "p.signif",paired = TRUE, method = 't.test',label.x = 1.5)+
    ylab('memory B cell / B cell (%)') + xlab('') + NoLegend() + ggtitle('memory B')

# memory CD8 t cell
tmp1<- cd8_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY')
tmp2<- cd8_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY')
cd8_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                        before=tmp2$Ratio, after=tmp1$Ratio )
p2 <- cd8_pair_df %>% filter(subtype == 'T.CD8.mem') %>%
    ggpaired(cd8_pair_df, cond1 = 'before', cond2 = 'after',legend = 'right',point.size = 3,
             color = "condition", line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(aes(label = paste0('p =', ..p.format..)),paired = TRUE, method = 't.test',label.x = 1.5)+
    ylab('CD8 T cell / CD8 cell (%)') + xlab('')+ ggtitle('memory CD8 T')

p1|p2

