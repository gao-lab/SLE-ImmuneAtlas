library(cowplot)
library(harmony)
library(tidyverse)
# library(RISC)
library(Seurat)
library(RColorBrewer)

setwd('/data/sle')
output_path <- './output_file/seurat/b_cell/'
source('./scripts/function_R/utils.R')

# load('./output_file/seurat/b_cell/bcell_subset_01.rdata')

#------------------------------ Re Analysis ------------------------------------
# bcell <- do_seurat(bcell)
# 
# DimPlot(bcell, label = T)


#------------------------------ Harmony Batch Remove ---------------------------
# bcell_harm <- bcell %>% 
#     RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
# 
# bcell_harm <- bcell_harm %>% 
#     RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#     FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#     FindClusters(resolution = 0.8) %>% 
#     identity()
back_run(do_harmony,out_name = 'bcell_harm', job_name = 'bcell_harm',
         seu_obj = bcell, harmony_slot = 'orig.ident',max.iter = 30, res = c(0.8, 1.0), from_begin = T)

DimPlot(bcell_harm, label = T)
VlnPlot(bcell_harm, features = marker_list, stack = T)

save(bcell_harm, file = './final/addtional/seurat/b_cell/02-b_cell_raw_harm.rdata')

#------------------------------- Exclude plasma --------------------------------
bcell_filter <- subset(bcell_harm, idents = c(10, 13, 15, 6, 8), invert = T)
# back_run(do_harmony,out_name = 'bcell_filter', job_name = 'bcell_filter', from_begin = T,
#          seu_obj = bcell_filter, harmony_slot = 'orig.ident',max.iter = 30,res =c(0.6,1.0,1.2,0.8))
bcell_filter <- do_harmony(seu_obj = bcell_filter, harmony_slot = 'orig.ident',
                           max.iter = 30, res =c(0.8), from_begin = T)

DimPlot(bcell_filter, label = T)
VlnPlot(bcell_filter, features = c(marker_list, feats, c('MKI67')), stack = T)
bcell_filter <- subset(bcell_filter, idents = c(10, 11, 12, 13), invert = T)

DimPlot(bcell_filter,group.by = 'subtype', label = T, cols = get_color(len = 7), pt.size = 0.1) + NoAxes() 

VlnPlot(bcell_filter, features = marker_list,stack = T) + NoLegend()
VlnPlot(bcell_filter, features = b_marker_list,stack = T) + NoLegend()
VlnPlot(bcell_filter, features = feats,stack = F) + NoLegend()
DotPlot2(bcell_filter, marker_list = marker_list)

#------------------------- Vis the Marker Gene list ----------------------------
DimPlot(bcell_filter, label = T) + DimPlot(bcell_filter, label = T,group.by = 'group') & NoLegend()

# marker list from biorxiv: Multi-regional characterization of renal cell 
# carcinoma and micro environment at single cell resolution

DotPlot2(bcell_filter, marker_list = c('TLR7','ISG15','IFI27','MME'))
DotPlot2(bcell_filter, marker_list = b_marker_list)
# DotPlot(bcell_filter, features = c('CD38','AICDA','TLR7'),group.by = 'RNA_snn_res.0.8')

plot_scdata(bcell_filter, color_by = "treatment", pal_setup = 'Set2') +
    plot_scdata(bcell_filter, color_by = "subtype", pal_setup = pal)

DimPlot(bcell_filter,label = T) + DimPlot(bcell_filter, group.by = 'subtype', label = T) 
DimPlot(bcell_filter, group.by = 'Phase')
DotPlot2(bcell_filter, marker_list = c('IL21R' ,'BACH2','BCL6','AICDA','IRF4' ,'PRDM1','PAX5','IFNAR1','IFNAR2','CD19'), group.by = 'seurat_clusters')

DotPlot2(bcell_filter, group.by = 'seurat_clusters', marker_list = c('PAX5','IL4R','TCL1A','IFIT1','IFIT3','ISG15','TLR7','IFI27','CD5','CD9','MME','IGHA2','IGHG1','AICDA'))
DotPlot2(bcell_filter,marker_list = c('IFIT1','IFIT3','ISG15','IFI27','TLR7'))
DotPlot2(subset(bcell_filter, idents = 5),marker_list = c('IFIT1','IFIT3','ISG15','IFI27','TLR7'), group.by = 'treatment')
DotPlot2(bcell_filter,marker_list = c('IFIT1','IFIT3','ISG15','IFI27','TLR7'), group.by = 'treatment')
DotPlot2(bcell_filter,marker_list = c('IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHA3','IGHE','IFIT1','IFIT3','ISG15'), group.by = 'treatment')
VlnPlot(bcell_filter, features =  c('IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHA3','IGHE','IFIT1','IFIT3','ISG15'), stack = T,group.by = 'treatment')

DotPlot2(bcell_filter, marker_list = b_marker_list)

#----------------------------- Finder Markers ----------------------------------
# back_run(FindMarkers, out_name = 'marker_b_sle.hc',job_name = 'marker_b_sle.hc',
#          bcell_filter,group.by = 'treatment', ident.1 = 'untreated', ident.2 = 'HC',
#          only.pos = T, min.pct = 0.3,logfc.threshold = 0.35)
# marker_b_sle.hc %<>% filter(p_val_adj <0.05) %>% arrange(-avg_log2FC)
# # marker_all_bcell_filter <- FindAllMarkers(object = bcell_filter, only.pos = T, logfc.threshold = 0.25)
# back_run(func =FindAllMarkers,out_name = 'marker_all_bcell_filter',job_name = 'marker_all_bcell_filter',
#          object = bcell_filter, only.pos = T, logfc.threshold = 0.25)
# marker_all_bcell_filter %<>% filter(p_val_adj < 0.05)
# marker_all_bcell_filter %>% group_by(cluster) %>% top_n(avg_log2FC, n = 10) %>% View()
# 
# marker_b_.c10 <- FindMarkers(bcell_filter, ident.1 = c(10),
#                              logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
# marker_b_.c10 %<>% filter(p_val_adj <0.05)  %>% arrange(-avg_log2FC)
# 
marker_b_.c12 <- FindMarkers(bcell_filter, ident.1 = c(12),
                              logfc.threshold = 0.25, min.pct = 0.2, only.pos = T)
marker_b_.c12 %<>% filter(p_val_adj <0.05) %>% arrange(-avg_log2FC)

#---------------------------- Anno the Sub type --------------------------------
bcell_filter$subtype <- 'unknown'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(2))] <- 'B.transition'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(0,1))] <- 'B.naive'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(3,9))] <- 'B.mem.CXCR3+'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(5))] <- 'B.mem.CD27-'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(7))] <- 'B.mem.IGHM+'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(6,8))] <- 'B.mem'
bcell_filter$subtype[which(bcell_filter$RNA_snn_res.0.8 %in% c(4))] <- 'B.IFN-response'

table(bcell_filter$subtype)
DimPlot(bcell_filter, group.by = 'subtype', label = T) + NoLegend()

bcell_filter <- subset(bcell_filter, idents = 9, invert = T)


# pub plot 
bcell_filter$subtype <- factor(bcell_filter$subtype, 
                               levels = rev(c("B.transition","B.naive","B.IFN-response",
                                              "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-")))
DotPlot2(bcell_filter, marker_list = c('MME','IGHD','CXCR4','CCR7','TCL1A','IGHM','IFITM1','STAT1','CD40','CD24','CD22','CD69','VPREB3',
                                       'CD27','IGHG1','CXCR3','IGHA1','IGHA2','MS4A1','FGR','FCRL5')
         , group.by = 'subtype') +   theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))
# pub plot 
bcell_filter$subtype <- factor(bcell_filter$subtype, 
                               levels = rev(c("B.transition","B.naive","B.IFN-response",
                                              "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-")))
DotPlot2(bcell_filter, marker_list = c('STAT2','STAT1','IRF9','IFIT1','IFI44L','IFI6','ISG15',
                                       'EPSTI1','OAS2','IRF1','MX1','IRF7','MX2','IFIT3','IFIT2','ISG20')
         , group.by = 'subtype') +   theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))

#------------------------ Ratio of B cell in Sample ----------------------------
# Bar plot :B cell clusters distribution in subtypes
ggplot(data = bcell_filter@meta.data, aes(x = bcell_filter$treatment, 
                                          fill =bcell_filter$subtype))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(bcell_filter$subtype))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values=get_color(7,set_len = 8)) + theme_bw() 

# Cluster ratio by group
bcell_filter@meta.data  %>% group_by(orig.ident,seurat_clusters) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                      palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
    facet_wrap(~seurat_clusters,scales = "free")+ stat_compare_means(method = 't.test')

# pub
# subtype  ratio  by group
bcell_filter@meta.data  %>% group_by(orig.ident,subtype) %>% 
    summarise(sub_num = n()) %>% mutate(sample_num = sum(sub_num)) %>% 
    mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(treatment != 'treated') %>%
    mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
                                                   "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07", "#3CB371"))+ 
    facet_wrap(subtype~.,scales = "free", ncol = 4)+ stat_compare_means(method = 't.test',label.x = 1.4,label = 'p.format')

# Bar plot :cell type sample distribution 
bcell_filter@meta.data %>% ggplot( aes(x = subtype , fill = orig.ident))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels=  unique(bcell_filter$orig.ident)) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

# Box plot :cell subtype ratio between groups
bcell_filter@meta.data %>% group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free")  + stat_compare_means(label.x = 1.2,method = 't.test')

# Box plot :(unpaired) cell subtype ratio between treatment
bcell_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>%
    ggpubr::ggboxplot(x='treatment',y='Ratio', color = 'treatment',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free") + stat_compare_means( method = 't.test',label.x = 1.2)

# pub
# Box plot :(paired)cell subtype ratio between treatment
tmp1<- bcell_filter@meta.data  %>% 
    group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')  %>% filter(!subtype =='unknown') 
tmp2<- bcell_filter@meta.data  %>% 
    group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')%>% filter(!subtype =='unknown')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>%
    mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
                                             "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
              fill  = "condition", line.color = "gray", line.size = 0.4,
              palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,label = 'p.format')+
    ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 4 )

# Heatmap: sample clustering based on B cell sub type 
bcell_sample_ratio_df <- table(bcell_filter$subtype, bcell_filter$orig.ident) %>% 
    as.data.frame() %>% filter(Var2 !='pSS_pah') %>%group_by(Var2) %>% 
    mutate(Frequency = sum(Freq)) %>% mutate(Ratio = Freq/Frequency*100) %>%
    select(Var1,Var2,Ratio) %>% dcast(Var1~Var2) %>% column_to_rownames('Var1')

heatmap(as.matrix(bcell_sample_ratio_df))
# heatmap(as.matrix(bcell_sample_ratio_df),scale = c('col'))


#------------------------ Save the Anno Object --------------------------------------
save(bcell_filter, file = './final/addtional/seurat/b_cell/03-b_cell_anno_filter_harm.rdata')

