library(Scillus)

setwd('/data/sle')
output_path <- './output_file/seurat/mono_dc/'
source('./scripts/function_R/utils.R')

#------------------------------ Re Analysis ------------------------------------
mono_dc_filter <- subset(mono_dc_harmony, idents = c(6,9,14,12,13,15,17), invert = T)
mono_dc_filter <- do_seurat(mono_dc_filter)
mono_dc_pic <- subset(mono_dc_harmony, idents = c(6,9,14,12,13,15,17))
save(mono_dc_pic,file=  paste0(output_path,'pic_mono_dc_harmony.rdata'))

#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
mono_dc_filter <- mono_dc_filter %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

mono_dc_filter <- mono_dc_filter %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()
save(mono_dc_filter, file = paste0(output_path, 'mono_dc_filter_harmony.rdata'))

DimPlot(mono_dc_filter,group.by = 'orig.ident') + DimPlot(mono_dc_filter, group.by = 'subtype', label = T) 
plot_scdata(mono_dc_filter, color_by = "treatment", pal_setup = 'Set2') + 
  plot_scdata(mono_dc_filter, color_by = "subtype", pal_setup = 'Set1')

# add meta
mono_dc_filter$disease <- 'SLE'
mono_dc_filter$disease[which(mono_dc_filter$group == 'HC')] <- 'HC'

#---------------------------- Vis the Marker Gene ------------------------------
DotPlot(mono_dc_filter, features = c('CD14','FCGR3A','CLEC4C',macro_sub_marker), group.by = 'subtype')+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('')

Idents(mono_dc_filter) <- 'subtype'
DotPlot(mono_dc_filter, features = c('IRF5','IRF9','IFNAR1','IFNAR2','TLR7','TLR9','IRF7'),
        split.by = 'treatment', cols = c("red", "blue", "green"))

#---------------------------- Vis the distribution------------------------------
tmp_df <- mono_dc_filter@meta.data %>% filter(group != 'pSS_pah')
ggplot(data = tmp_df, aes(x = tmp_df$subtype, fill =tmp_df$group))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tmp_df$group))) + 
  labs(fill="group") + coord_flip()


# load('./output_file/seurat/mono_dc/mono_dc_filter_harmony.rdata')
mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~subtype,scales = "free") 

mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'treatment',
                    palette =c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  facet_wrap(~subtype,scales = "free") 

#------------------------------- Ratio -----------------------------------------
# By Disease
mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
  filter(!treatment=='treated') %>%
  ggpubr::ggboxplot(x='disease',y='Ratio', fill = 'disease',
                    palette =c("#FC4E07", "#00AFBB"))+ 
  facet_wrap(~subtype,scales = "free") + stat_compare_means( method = 't.test',label.x = 1.2)

# paired test
tmp1<- mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')
tmp2<- mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')
mono_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                             before=tmp2$Ratio, after=tmp1$Ratio )
ggpaired(mono_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free" )


#------------------------------- Marker genes-----------------------------------
# marker_DC_pDC <- Find


#---------------------------------- Save ---------------------------------------
save(mono_dc_filter, file = './output_file/seurat/mono_dc/final_mono_dc.rdata')
