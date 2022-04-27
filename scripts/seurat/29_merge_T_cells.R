library(Seurat)

t_cell_merge <- merge(cd4_filter, y=c(cd8_filter,other_T_filter,nk_filter,prolife_T_fliter))
back_run(do_harmony,out_name = 't_cell_merge', job_name = 't_cell_merge',
         seu_obj = t_cell_merge, harmony_slot = 'orig.ident',max.iter = 30,res =c(1.0,0.8),from_begin = T)


DimPlot(t_cell_merge, group.by = 'subtype',label = T)

save(t_cell_merge, file = './final/seurat/t_cell/04-t_cell_merge.rdata')


# Analysis
t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07" ,"#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free",ncol = 6) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )

tmp1<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')  %>% filter(!subtype =='unknown') 
tmp2<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')%>% filter(!subtype =='unknown')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
              fill  = "condition", line.color = "gray", line.size = 0.4,
              palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,label = 'p.format')+
    ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 6)
