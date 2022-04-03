library(cowplot)
library(harmony)
# library(RISC)
library(Seurat)
library(tidyverse)
library(ggpubr)

setwd('/data/sle')
output_path <- './output_file/seurat/platelet/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/platelet/platelet_subset_01.rdata')

#------------------------------ Re Analysis ------------------------------------
platelet <- do_seurat(platelet)
DimPlot(platelet, label = T) + DimPlot(platelet, group.by = 'orig.ident')


#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
platelet_harmony <- platelet %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

platelet_harmony <- subset(platelet_harmony, idents = c(8,9 ), invert =10)

# save(platelet_harmony, file = paste0(output_path, 'platelet_harmony_14iter.rdata'))

DimPlot(platelet_harmony, label = T) + DimPlot(platelet_harmony, group.by = 'orig.ident')
plot_scdata(platelet_harmony, color_by = "treatment", pal_setup = 'Set2') + 
  plot_scdata(platelet_harmony, color_by = "seurat_clusters", pal_setup = 'Set1')


#------------------------------ Finder Marker ----------------------------------
back_run(func =FindAllMarkers,out_name = 'marker_all_platelet_filter',job_name = 'marker_all_platelet_filter',
         object = platelet_harmony, only.pos = T, logfc.threshold = 0.25)
marker_all_platelet_filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10) %>% View()

#----------------------------- Vis Marker Genes --------------------------------
VlnPlot(platelet_harmony, features = feats)

FeaturePlot(platelet_harmony, features = c('CD58','CD69','IFITM1'), pt.size = 0.7,
            cols = c('grey','red'), ncol = 3) &NoAxes()

DotPlot2(platelet_harmony, marker_list  =c('CD58','CD69','IFITM1','PRKRA'))
DotPlot(platelet_harmony, features =c('CD58','CD69','IFITM1','PRKRA') , 
        group.by = 'group')
DotPlot2(platelet_harmony, marker_list = c('CD58','CD69','IFITM1'), group.by = 'subtype')
VlnPlot(platelet_harmony, features =c('CD58','CD69','IFITM1'), stack = T , group.by = 'subtype')
marker_platelet_slepah <- FindMarkers(platelet_harmony, group.by = 'group', 
                                      ident.1 = 'SLE_pah')
marker_platelet_slepah <- marker_platelet_slepah%>% filter(p_val_adj<0.05)

# by cell type annotation result 
DotPlot(platelet_harmony, features = c('CD58','CD69','IFITM1','PRKRA'),
        group.by = 'subtype')

#---------------------------- Anno the Sub type --------------------------------
platelet_harmony$subtype <- 'Platelet'
# platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c(7))] <- 'Platelet.CD58+CD69+'
platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c(1,7,10))] <- 'Platelet.IFITM1+'
# platelet_harmony$subtype[which(platelet_harmony$seurat_clusters %in% c())] <- ''

table(platelet_harmony$subtype)
# pub
DimPlot(platelet_harmony, group.by = 'subtype', label = T, cols = get_color(2,set = 'Set2',set_len = 2)) + NoAxes()


#-------------------------- Compare at Sample Level ----------------------------
# plate to all pbmc

# By seurat_clusters
platelet_harmony@meta.data %>% 
  group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free") + stat_compare_means(method = 't.test',label.x = 1.4,label = 'p.format')

# By subtype 
platelet_harmony@meta.data %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                    palette =rev(c("#00AFBB", "#FC4E07")))+ 
  facet_wrap(~subtype ,scales = "free") + stat_compare_means(method = 't.test',label.x = 1.4,label = 'p.format')

# pub
# Box plot :(paired)cell subtype ratio between treatment
tmp1<- platelet_harmony@meta.data  %>% 
  group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')  %>% filter(!subtype =='unknown') 
tmp2<- platelet_harmony@meta.data  %>% 
  group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(platelet_harmony@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')%>% filter(!subtype =='unknown')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>%
  # mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
                                           # "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
  ggpaired( cond1 = 'before', cond2 = 'after',
            fill  = "condition", line.color = "gray", line.size = 0.4,
            palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,label = 'p.format')+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 4 )

# pub  gsea
gsea_reslut <-plot_gsea(platelet_harmony, group_by = 'subtype', focus = 'Platelet.IFITM1+',
          title = 'Platelet GSEA enrichment',category = 'H')


#-------------------------------- Save Files -----------------------------------
platelet_filter <- platelet_harmony
save(platelet_filter, file = './final/seurat/platelet/02-platelet_anno.rdata')
