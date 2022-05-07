library(cowplot)
library(harmony)
library(Seurat)
library(future)

setwd('/data/sle')
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./final//seurat/t_cell/03-CD4_Tcell_raw_harm.rdata')
DimPlot(cd8_filter, label = T,cols = get_color(12), pt.size = 0.3) + NoAxes() 
# DimPlot(cd8_filter, label = T) + DimPlot(cd8_filter,group.by = 'seurat_clusters')
VlnPlot(cd8_filter,features = c('CD4','CD8A','CD8B'), ncol = 3,stack = T)

# filter cluster 15
cd8_filter <- subset(cd8_filter, idents =c(8,11), invert =T)
select.cells <- CellSelector(plot = DimPlot(cd8_filter, label = T))
cd8_filter <-subset(cd8_filter,cells = select.cells, invert = T)

#--------------------------------- Vis Marker ----------------------------------
VlnPlot(cd8_filter, features = c('CD4','CD8A','CD8B','PDCD1'))
DotPlot2(cd8_filter, marker_list = marker_list)
DotPlot2(cd8_filter, marker_list = t_cd8_all_marker)
DotPlot2(cd8_filter, marker_list = c('GZMB','GZMH','NKG7','CCL5','IFIT3',
                                       'CCR7','CD27','SELL','TCF7',
                                       'IL7R','CD28','GPR183',
                                       'TBX21','PRF1','CX3CR1','KLRG1','EOMES','FGFBP2',
                                       'CXCR4','GZMK',
                                       'XCL1','XCL2','MYADM','CD6',
                                       'KLRB1','RORC','SLC4A10',
                                       'CCL3','CCL4','IFNG'))

DotPlot2(cd8_filter, marker_list = c('IFIT3','ISG15','MX1','OAS1','IFI44L'))

#------------------------------- Finder markers -----------------------------------
plan('multiprocess', workers = 10)
options(future.globals.maxSize= 100 * 1000 * 1024 ^2 ) 
# marker_all_cd8 <- FindAllMarkers(object = cd8_filter, only.pos = T, logfc.threshold = 0.25)
back_run(func =FindAllMarkers,out_name = 'marker_all_cd8',job_name = 'marker_all_cd8',
         object = cd8_filter, only.pos = T, logfc.threshold = 0.25)
marker_all_cd8.filter <-  marker_all_cd8 %>% filter(p_val_adj <0.05)
marker_all_cd8.filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10) %>% View()

marker_cd8_c6 <- FindMarkers(cd8_filter, ident.1 = 6, only.pos = T,logfc.threshold = 0.25)
marker_cd8_c6 %<>% filter(p_val_adj <0.05) %>% arrange(-avg_log2FC )

marker_cd8_c8 <- FindMarkers(cd8_filter, ident.1 = 8, only.pos = T,logfc.threshold = 0.25)
marker_cd8_c8 %<>% filter(p_val_adj <0.05) %>% arrange(-avg_log2FC )

marker_cd8_c10 <- FindMarkers(cd8_filter, ident.1 = 10, only.pos = T,logfc.threshold = 0.25) %>% 
  filter(p_val_adj <0.05) %>% arrange(-avg_log2FC )


#---------------------------------- Anno ---------------------------------------
cd8_filter$subtype <- 'unknown'
cd8_filter$subtype[which(cd8_filter$seurat_clusters %in% c(9))] <- 'T.CD8.Tex'
cd8_filter$subtype[which(cd8_filter$seurat_clusters %in% c(4,5,6))] <- 'T.CD8.Naive'
cd8_filter$subtype[which(cd8_filter$seurat_clusters %in% c(10))] <- 'T.CD8.IFN-response'
cd8_filter$subtype[which(cd8_filter$seurat_clusters %in% c(2,7))] <- 'T.CD8.mem'
cd8_filter$subtype[which(cd8_filter$seurat_clusters %in% c(0,1,3))] <- 'T.CD8.Teff'

DimPlot(cd8_filter, label = T, group.by = 'subtype',cols = get_color(5), pt.size = 0.3) + NoAxes()
Idents(cd8_filter) <- 'subtype'

#---------------------------------- Ratio---------------------------------------
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# overall  ratio
cd8_filter@meta.data %>% filter(group != 'pSS_pah') %>% 
  ggplot( aes(x = disease , fill = seurat_clusters))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(cd8_filter$seurat_clusters))) + 
  scale_fill_manual(values =getPalette(17)) +
  labs(fill="group") + coord_flip() 

# By subtype  
cd8_filter@meta.data  %>% group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                    palette =c("#FC4E07","#00AFBB"),legend = 'right')+ 
  facet_wrap(~subtype,scales = "free", ncol = 5) + 
  stat_compare_means( method = "t.test" ,label.x = 1.2)

# paired test
tmp1<- cd8_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')
tmp2<- cd8_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                          before=tmp2$Ratio, after=tmp1$Ratio ) %>%
  ggpaired(cond1 = 'before', cond2 = 'after',fill = "condition", line.color = "gray", 
            line.size = 0.4,palette = "npg",legend = 'right') +  
  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,tip.length = T) +
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free", ncol = 5 )
rm(tmp1,tmp2)

## test
cd8_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd8_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(treatment =='treated')%>% filter(!orig.ident == 'XYY2') %>%
  ggpaired( x = "supp", y = "len",
           color = "supp", line.color = "gray", line.size = 0.4,
           palette = "npg") + facet_wrap(~subtype,scales= "free" )

#---------------------------------- Save ---------------------------------------
save(cd8_filter,file = './final/seurat/t_cell/04-CD8_Tcell_filter_anno.rdata')

