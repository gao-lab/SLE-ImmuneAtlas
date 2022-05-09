library(cowplot)
library(harmony)
library(RISC)
library(Seurat)
library(tidyverse)
library(ggpubr)

setwd('/data/sle')
output_path <- './output_file/seurat/mono_dc/'
source('./scripts/function_R/utils.R')

load('./final//seurat/mono_dc/01-mono_dc_raw.rdata')

#---------------------------- Harmony Batch Remove -----------------------------

back_run(do_harmony,out_name = 'mono_dc_harm', job_name = 'mono_dc_harm',
         seu_obj = mono_dc_harm, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.8,1.0), from_begin =T)
save(mono_dc_harm, file = 'final/seurat/mono_dc/02-mono_dc_harm.rdata')

#---------------------------- Remove Doublets -----------------------------
DimPlot(mono_dc_harm, label = T) + DimPlot(mono_dc_harm, group.by = 'group')
mono_dc_filter <- subset(mono_dc_harm, idents = c(6,11,15,18,20,9,16), invert = T)
DimPlot(mono_dc_filter, label = T) + DimPlot(mono_dc_filter, group.by = 'group')

mono_index <- CellSelector(DimPlot(mono_dc_filter))
mono_dc_filter <- subset(mono_dc_filter, cells = mono_index, invert = T)

marker_mono_c19 <- FindMarkers(mono_dc_harm, ident.1 = 19, only.pos = T, min.pct = 0.25)
marker_mono_c19   %<>%  filter(p_val_adj < 0.05) %>% arrange(-avg_log2FC )

DotPlot(mono_dc_filter, features = marker_list)+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
FeaturePlot(mono_dc_harm, features = feats[1:2])

DotPlot2(mono_dc_filter, marker_list = c('CD14','FCGR3A')) 
## HSPC MAST cell marker 
DotPlot2(mono_dc_filter, marker_list = c('TBSB2','TPSAB1','CD34','CD38')) 
## Macrophage cell marker 
DotPlot2(mono_dc_filter, marker_list = macro_sub_marker)
## Mono DC cell marker 
DotPlot2(mono_dc_filter, marker_list = macro_sub_marker)

# pub (ratio)
ggplot(data = mono_dc_filter@meta.data, aes(x = treatment, 
                                          fill =subtype))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(bcell_filter$subtype))) + xlab('')+
  labs(fill="") +  scale_fill_manual(values=get_color(10,set_len = 8)) + theme_bw() 

table(mono_dc_filter$seurat_clusters, mono_dc_filter$group) %>% 
  as.data.frame()  %>%group_by(Var2) %>% 
  mutate(Frequency = sum(Freq)) %>% mutate(Ratio = Freq/Frequency*100) %>%
ggplot() + geom_bar(mapping = aes(x = Var1 , y = Ratio, fill= Var2),
           stat = "identity",
           position = "dodge")
table(mono_dc_filter$seurat_clusters, mono_dc_filter$treatment) %>% 
  as.data.frame()  %>%group_by(Var2) %>% 
  mutate(Frequency = sum(Freq)) %>% mutate(Ratio = Freq/Frequency*100) %>%
  ggplot() + geom_bar(mapping = aes(x = Var1 , y = Ratio, fill= Var2),
                      stat = "identity",
                      position = "dodge")

#------------------------------- Finder markers -----------------------------------
plan('multiprocess', workers = 10)
options(future.globals.maxSize= 100 * 1000 * 1024 ^2 ) 
# marker_all_mono.dc.filter <- FindAllMarkers(object = mono_dc_filter, only.pos = T, logfc.threshold = 0.25)
back_run(func =FindAllMarkers,out_name = 'marker_all_mono.dc.filter',job_name = 'marker_all_mono.dc.filter',
         object = mono_dc_filter, only.pos = T, logfc.threshold = 0.25)
marker_all_mono.dc.filter %<>% filter(p_val_adj <0.05)
marker_all_mono.dc.filter %>% group_by(cluster) %>%top_n(avg_log2FC  , n = 10) %>% 
  arrange(desc(cluster),avg_log2FC)%>% View()

marker_mono.dc.filter_c034 <- FindMarkers(mono_dc_filter, ident.1 = c(0,3,4), only.pos = T,
                                          logfc.threshold = 0.25, min.diff.pct = 0.2)
marker_mono.dc.filter_cdc1.cdc2 <-  FindMarkers(mono_dc_filter,group.by = 'subtype', 
                                                ident.1 = 'cDC1',ident.2 = 'cDC2', only.pos = F,
                                                logfc.threshold = 0.25, min.diff.pct = 0.2)
marker_mono.dc.filter_cdc1.cdc2 %>%filter(p_val_adj <0.05) %>%slice_min(avg_log2FC  , n = 10) %>% 
  arrange(-avg_log2FC)%>% View()

plan('sequential')
#----------------------------- Anno Cell Sub type ------------------------------
# cdc1 marker CD1C and cdc2 marker CLEC9A
mono_dc_filter$subtype <- 'unknown'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(1,2,5))] <- 'Mono.CD14'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(0,10))] <- 'Mono.CD14.LGALS2+'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(3))] <- 'Mono.CD14.APOBEC3A+'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(4))] <- 'Mono.CD14.recruit'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(7,8,14))] <- 'Mono.CD16'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(17))] <- 'HPSC'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(12))] <- 'cDC1'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(19))] <- 'cDC2'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(13))] <- 'pDC'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c(14))] <- 'Macrophage'
mono_dc_filter$subtype[which(mono_dc_filter$seurat_clusters %in% c())] <- 'doubelt'
table(mono_dc_filter$subtype)

# pub
DimPlot(mono_dc_filter,group.by = 'subtype', label = T,cols = get_color(len = 10,set_len = 11,set = 'set1'),pt.size = 0.1) + NoAxes() 


#-------------------------- Compare at Sample Level ----------------------------
# Barplot prolife (NO prolife)
ggplot(mono_dc_filter@meta.data, aes(x=seurat_clusters, fill=Phase))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(names(table(mono_dc_filter$seurat_clusters))) + 
  labs(fill="Phase") + coord_flip()

# Cluster by group
mono_dc_filter@meta.data  %>% group_by(orig.ident,seurat_clusters) %>% 
  summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
  mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  filter(treatment != 'treated') %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free")+ stat_compare_means()

# pub
# By Group
mono_dc_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  filter(!treatment == 'treated') %>%
  ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                    palette =c("#DA9494", "#B4D493"))+ 
  facet_wrap(~subtype,scales = "free",nrow = 2) + stat_compare_means(label.x = 1.2)+
  theme_cowplot()


# By Treatment
mono_dc_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_harm@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'treatment',
                    palette =c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  facet_wrap(~subtype,scales = "free") 

# Box plot :(paired)cell subtype ratio between treatment
tmp1<- mono_dc_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')
tmp2<- mono_dc_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                             before=tmp2$Ratio, after=tmp1$Ratio ) %>%
ggpaired( cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 5 )

#------------------------ Save the Anno Objcet --------------------------------------
save(mono_dc_filter, file = './final/seurat/mono_dc/03-mono_dc_anno_filter_harm.rdata')
