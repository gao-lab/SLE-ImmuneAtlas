library(cowplot)
library(harmony)
library(RISC)
library(Seurat)
library(tidyverse)
library(ggpubr)

setwd('/data/sle')
output_path <- './output_file/seurat/mono_dc/'
source('./scripts/function_R/utils.R')

# load('./output_file/seurat/mono_dc_subset_01.rdata')

#------------------------------ Re Analysis ------------------------------------
# mono_dc <- do_seurat(mono_dc)
# DimPlot(mono_dc, label = T) + DimPlot(mono_dc, group.by = 'orig.ident')


#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
# mono_dc_harmony <- mono_dc %>% 
#   RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
# 
# mono_dc_harmony <- mono_dc_harmony %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.8) %>% 
#   identity()
back_run(do_harmony,out_name = 'mono_dc_harm', job_name = 'mono_dc_harm',
         seu_obj = mono_dc, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.8,1.0))


save(mono_dc_harmony, file = paste0(output_path, 'mono_dc_harmony_11iter.rdata'))

DimPlot(mono_dc_harmony, label = T) + DimPlot(mono_dc_harmony, group.by = 'group')

DotPlot(mono_dc_harmony, features = marker_list)+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FeaturePlot(mono_dc_harmony, features = feats[1:2])

## MAST cell marker 
DotPlot(mono_dc_harmony, features = c('TBSB2','TPSAB1'))+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Macrophage cell marker 
DotPlot(mono_dc_harmony, features = c('CD163','C1QC','SPP1'))+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('')

## Mono DC cell marker 
DotPlot(mono_dc_harmony, features = macro_sub_marker)+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('')

FeaturePlot(mono_dc_harmony, features = feats[1:2])


table(mono_dc_harmony$seurat_clusters, mono_dc_harmony$group) %>% 
  as.data.frame() %>% filter(Var2 !='pSS_pah') %>%group_by(Var2) %>% 
  mutate(Frequency = sum(Freq)) %>% mutate(Ratio = Freq/Frequency*100) %>%
ggplot() + geom_bar(mapping = aes(x = Var1 , y = Ratio, fill= Var2),
           stat = "identity",
           position = "dodge")
table(mono_dc_harmony$seurat_clusters, mono_dc_harmony$treatment) %>% 
  as.data.frame() %>% filter(Var2 !='pSS_pah') %>%group_by(Var2) %>% 
  mutate(Frequency = sum(Freq)) %>% mutate(Ratio = Freq/Frequency*100) %>%
  ggplot() + geom_bar(mapping = aes(x = Var1 , y = Ratio, fill= Var2),
                      stat = "identity",
                      position = "dodge")


#----------------------------- Anno Cell Sub type ------------------------------
mono_dc_harmony$subtype <- 'unknown'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(0,1,2))] <- 'CD14.Mono'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(4))] <- 'TR.Mono'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(5,8))] <- 'CD16.Mono'
# mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(8))] <- ''
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(3))] <- 'Pro-infla.Mono'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(7))] <- 'mDC'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(10))] <- 'pDC'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(16))] <- 'MAST'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(11))] <- 'Macrophage'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(6,9,14))] <- 'PIC_Mono_T'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(12))] <- 'PIC_Mono14_Platelet'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(13))] <- 'PIC_Mono16_Platelet'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(15))] <- 'PIC_mDC_T'
mono_dc_harmony$subtype[which(mono_dc_harmony$seurat_clusters %in% c(17))] <- 'PIC_pDC_T'
table(mono_dc_harmony$subtype)
DimPlot(mono_dc_harmony, group.by = 'subtype', label = T) + NoLegend()


#-------------------------- Compare at Sample Level ----------------------------
# By Group
mono_dc_harmony@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_harmony@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~subtype,scales = "free") 

# By Treatment
mono_dc_harmony@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(mono_dc_harmony@meta.data[,c(1,12,13,14)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'treatment',
                    palette =c("#00AFBB", "#E7B800", "#FC4E07"))+ 
  facet_wrap(~subtype,scales = "free") 

