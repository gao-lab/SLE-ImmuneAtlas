

library(cowplot)
library(harmony)
library(RISC)
library(future)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./output_file/seurat/t_cell/cd8_tcell_filter_harmony_01.rdata')
DimPlot(t.cd8.filter, label = T,cols = getPalette(17), pt.size = 0.3) + NoAxes() 
# DimPlot(t.cd8.filter, label = T) + DimPlot(t.cd8.filter,group.by = 'seurat_clusters')
VlnPlot(t.cd8.filter,features = c('CD4','CD8A','CD8B'), ncol = 3)

DotPlot2(t.cd8.filter, marker_list = marker_list)

t.double <- subset(t.cd8.filter, idents =c(7,11,14))
t.double_harmony <- do_harmony(t.double,from_begin = T,max.iter = 20,res = c(0.5,0.8))
back_run(func =do_harmony, out_name = 't.double_harmony',job_name = 't.double_harmony', 
         seu_obj = t.double, from_begin = T,max.iter = 30,res = c(0.5,0.8))

# add meta
t.cd8.filter$disease <- 'SLE'
t.cd8.filter$disease[which(t.cd8.filter$group == 'HC')] <- 'HC'

# filter cluster 15
select.cells <- CellSelector(plot = DimPlot(t.cd8.filter, label = T))
t.cd8.filter <-subset(t.cd8.filter,cells = select.cells, invert = T)

#--------------------------------- Vis Marker ----------------------------------
VlnPlot(t.cd8.filter, features = c('CD4','CD8A','CD8B','PDCD1'))
DotPlot2(t.cd8.filter, marker_list = t_cd8_all_marker)
DotPlot2(t.cd8.filter, marker_list = c('GZMB','GZMH','NKG7','CCL5','IFIT3',
                                       'CCR7','CD27','SELL','TCF7',
                                       'IL7R','CD28','GPR183',
                                       'TBX21','PRF1','CX3CR1','KLRG1','EOMES','FGFBP2',
                                       'CXCR4','GZMK',
                                       'XCL1','XCL2','MYADM','CD6',
                                       'KLRB1','RORC','SLC4A10',
                                       'CCL3','CCL4','IFNG'))

#------------------------------- Finder DEGs -----------------------------------
plan('multiprocess', workers = 16)
options(future.globals.maxSize= 200 * 1000 * 1024 ^2 ) 
# marker_all_cd8 <- FindAllMarkers(object = t.cd8.filter, only.pos = T, logfc.threshold = 0.25)
back_run(func =FindAllMarkers,out_name = 'marker_all_cd8',job_name = 'marker_all_cd8',
         object = t.cd8.filter, only.pos = T, logfc.threshold = 0.25)
marker_all_cd8.filter <-  marker_all_cd8 %>% filter(p_val_adj <0.05)
marker_all_cd8.filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10) %>% View()

#---------------------------------- Anno ---------------------------------------
# 5 6 7 9 
# cluster15 is red blood cell
t.cd8.filter$subtype <- 'unknown'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(10))] <- 'CD8_Tex2_Tcm'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(8))] <- 'CD8_Tex1'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(2,14,6,9))] <- 'CD8_Naive'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(0,1,16))] <- 'CD8_Temra'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(11))] <- 'CD8_Tcm'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(4))] <- 'CD8_Temra_active'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(12))] <- 'CD8_IFN-response_Teff'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(13))] <- 'CD8_IFN-response_Naive'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(5))] <- 'CD8_mem_GZMK'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(3))] <- 'CD8_Teff_GZMK'
# t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(9))] <- 'Tfh'
t.cd8.filter$subtype[which(t.cd8.filter$seurat_clusters %in% c(7))] <- 'CD8_low_quality'
DimPlot(t.cd8.filter, label = T, group.by = 'subtype',cols = getPalette(11), pt.size = 0.3) + NoAxes()

Idents(t.cd8.filter) <- 'subtype'


#---------------------------------- Ratio---------------------------------------
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# overall plasma ratio
t.cd8.filter@meta.data %>% filter(group != 'pSS_pah') %>% 
  ggplot( aes(x = disease , fill = seurat_clusters))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(t.cd8.filter$seurat_clusters))) + 
  scale_fill_manual(values =getPalette(17)) +
  labs(fill="group") + coord_flip() 

# By Group
t.cd8.filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t.cd8.filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='disease',y='Ratio', color = 'disease',
                    palette =c("orange", "green"))+ 
  facet_wrap(~subtype,scales = "free",ncol = 4) + 
  stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test",vjust = T )

# paired treatment
# paired test
tmp1<- t.cd8.filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t.cd8.filter@meta.data[,c(1,12,13,14,17,18)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')
tmp2<- t.cd8.filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t.cd8.filter@meta.data[,c(1,12,13,14,17,18)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')
pdc_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                          before=tmp2$Ratio, after=tmp1$Ratio )
ggpaired(pdc_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,tip.length = T)+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free" )
rm(tmp1,tmp2)

## test
t.cd8.filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t.cd8.filter@meta.data[,c(1,12,13,14,17,18)]  %>%  distinct()) %>%
  filter(treatment =='treated')%>% filter(!orig.ident == 'XYY2') %>%
  ggpaired( x = "supp", y = "len",
           color = "supp", line.color = "gray", line.size = 0.4,
           palette = "npg") + facet_wrap(~subtype,scales= "free" )

#---------------------------------- Save ---------------------------------------
save(t.cd8.filter,file = './output_file/seurat/t_cell/final_t.cd8.rdata')

