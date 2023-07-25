library(cowplot)
library(harmony)

setwd('/data/sle')
source('./scripts/function_R/utils.R')
output_path <- './final/seurat/t_cell/'

#------------------------------ Re Analysis ------------------------------------
# tcell <- do_seurat(seu_obj = tcell,res = c())
tcell <- back_run(do_seurat,out_name = 'tcell', job_name = 'tcell',
                  seu_obj = tcell,res = c(1.0,0.8))
# save(tcell, file = 'final/seurat/t_cell/02-t_cell_raw_cluster.rdata')

#------------------------------ Harmony Batch Remove ---------------------------

# tcell <- tcell %>% 
#   RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
# 
# tcell <- tcell %>% 
#   RunUMAP(reduction = "harmony", dims = 1:20) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
#   FindClusters(resolution = 0.8) %>% 
#   identity()
# do_harmony(seu_obj = tcell,harmony_slot = 'orig.ident',max.iter = 30,res =c(0.8,1.0) )
back_run(do_harmony,out_name = 'tcell_harm', job_name = 'tcell_harm',
         seu_obj = tcell_harm, harmony_slot = 'orig.ident',max.iter = 30,res =c(1.0,0.8),from_begin = T)
save(tcell_harm, file = 'final/seurat/t_cell/02-t_cell_raw_harm.rdata')

DimPlot(tcell_harm, reduction = "umap", label = TRUE, pt.size = .1)
DimPlot(tcell_harm, reduction = "umap", label = TRUE, group.by = 'Phase')
DimPlot(tcell_harm, reduction = "umap", label = TRUE, group.by = 'orig.ident')
VlnPlot(tcell_harm, features = feats, ncol = 2)

#----------------------------- Finder doublets ---------------------------------
# can not run this in the docker and I do it on jupyter with R 4.0.3
if(F){
library(sceasy)
use_condaenv('sceasy')
tcell_tmp <- tcell_harm
tcell_tmp@assays$RNA@data <- tcell_tmp@assays$RNA@counts
convertFormat(tcell_tmp, from="seurat", to="anndata",
                      outFile='./final/seurat/t_cell/02-t_cell_raw_harm.sceasy.h5ad')
rm(tcell_tmp)
}
# library(feather)

tcell_doublet_index <- read.csv('./scripts///scrublet/tcell_doublet_result.csv')
tcell_harm$scrublet_doublet <- tcell_doublet_index$doublet

# Barplot : doublet ratio of clusters
DimPlot(tcell_harm, group.by = 'scrublet_doublet',raster=FALSE) + DimPlot(tcell_harm,label = T) & NoLegend()
DimPlot(tcell_harm, group.by = 'RNA_snn_res.0.8',raster=T, label = T, split.by = 'RNA_snn_res.0.8')
ggplot(data = tcell_harm@meta.data, aes(x = RNA_snn_res.0.8, fill =scrublet_doublet))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'doublet proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tcell_harm$RNA_snn_res.0.8))) + coord_flip()

VlnPlot(tcell_harm, features = c('CD4','CD8A','CD8B',t_cell_marker),stack = T)
VlnPlot(tcell_harm, features = c('CD4','CD8A','CD8B',marker_list),stack = T)

#-----------------------------  Finder Markers ---------------------------------
library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize= 100 * 1000 * 1024 ^2 )
marker_tcell_c19 <- FindMarkers(tcell_harm, ident.1 = 19, min.pct = 0.25,min.diff.pct = 0.1)
marker_tcell_c19 <- marker_tcell_c19 %>% filter(p_val_adj < 0.05) %>% arrange(-avg_log2FC )

marker_tcell_c15 <- FindMarkers(tcell_harm, ident.1 = 15, min.pct = 0.25,min.diff.pct = 0.1, only.pos = T)
marker_tcell_c15 <- marker_tcell_c15 %>% filter(p_val_adj < 0.05) %>% arrange(-avg_log2FC )

plan("sequential")

#------------------------------- Vis Marker Genes-------------------------------
# Barplot prolife
ggplot(tcell_harm@meta.data, aes(x=seurat_clusters, fill=Phase))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(names(table(tcell_harm$seurat_clusters))) + 
  labs(fill="Phase") + coord_flip()


DotPlot(tcell_harm, features = c('CD4','CD8A','CD8B')) +
  VlnPlot(tcell_harm, features = c('CD4','CD8A','CD8B'), stack = T)

FeaturePlot(tcell_harm, features = c('CD4','CD8A','CD8B'), ncol = 3) + 
  DimPlot(tcell_harm, reduction = "umap", label = TRUE, pt.size = .1)

VlnPlot(tcell_harm, features = feats,stack = F, ncol = 2)
# invariant NKT and MAIT and yd T
DotPlot(tcell_harm, features = c(t_invarient_NKT_marker, t_MAIT_marker, t_yd_marker)) +
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  ylab('') + xlab('')

# check TRA and TRB
DotPlot(tcell_harm, features = grep(rownames(tcell_harm),pattern = "^TRA[VDJ]", value = T, invert = F, perl = T)) +
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('') + xlab('')
DotPlot(tcell_harm, features = grep(rownames(tcell_harm),pattern = "^TRB[VDJ]", value = T, invert = F, perl = T)) +
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  scale_size(breaks = c(0, 25, 50, 75)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('') + xlab('')

#---------------------------- Anno the Sub type --------------------------------
# cluster  ratio is correlated with cell number
tcell_harm@meta.data %>% group_by(orig.ident) %>% filter(seurat_clusters == 11) %>% 
  summarise(nkt=n()) %>% 
  left_join(y = all_cell@meta.data%>%group_by(orig.ident)%>% summarise(all=n()))%>%
  mutate(ratio = nkt/all*100) %>% ggplot(aes(x=all,y=ratio)) +  geom_point()+ stat_smooth(method = 'lm')

# Boxplot : T cell subtype radio in HC and SLE by cluster 
tcell_harm@meta.data  %>% group_by(orig.ident,seurat_clusters) %>% 
  summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(tcell_harm@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  # filter(treatment != 'treated') %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free") + stat_compare_means(method = 't.test')


tcell_harm$subtype <- 'unknown'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(0,3,6,17))] <- 'T.CD8.cyto'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(4,14))] <- 'T.CD8.naive'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(1,2,8))] <- 'T.CD4.naive'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(5,16))] <- 'NK'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(9,12))] <- 'T.CD4.cyto'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(13,18))] <- 'T.prolife'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(7))] <- 'T.yd'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(10))] <- 'MAIT'
# tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c())] <- 'low_quality'
tcell_harm$subtype[which(tcell_harm$seurat_clusters %in% c(11,15,19))] <- 'doublet'
table(tcell_harm$subtype)
DimPlot(tcell_harm, group.by = 'subtype', label = T) + NoLegend()

# beautiful 
DimPlot(tcell_harm, group.by = 'subtype', label = T,
        pt.size = 0.6, cols = brewer.pal(10, "Set3"),raster=FALSE) + NoLegend() + NoAxes()


#-------------------------- Celltype distribution ------------------------------
# Barplot : T cell clusters distribution in groups
ggplot(data = tcell_harm@meta.data, aes(x = subtype, fill =group))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tcell_harm$group))) + 
  labs(fill="group") + coord_flip()

# Barplot : T cell clusters distribution in sample
ggplot(data = tcell_harm@meta.data, aes(x =orig.ident, fill =subtype))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tcell_harm$subtype))) + 
  labs(fill="group") + coord_flip()

# Boxplot : T cell subtype radio in HC and SLE by group 
tcell_harm@meta.data  %>% group_by(orig.ident,subtype) %>% 
  summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(tcell_harm@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  filter(treatment != 'treated') %>%
  ggpubr::ggboxplot(x='group',y='Ratio', color = 'group',
                    palette =c("#FC4E07", "#E7B800", "#00AFBB"))+ 
  facet_wrap(~subtype,scales = "free") + stat_compare_means(method = 't.test')

# Boxplot: (paired)T prolife/T ratio 
tmp1<- tcell_harm@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(tcell_harm@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')
tmp2<- tcell_harm@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(tcell_harm@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')
Tprolife_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                             before=tmp2$Ratio, after=tmp1$Ratio )
ggpaired(Tprolife_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free" )

#-------------------------------- save the file --------------------------------
#NOTE not remove the doublet 
save(tcell_harm, file = 'final/seurat/t_cell/03_t_cell_anno_harm.rdata')

#------------------------------ RISC Batch Remove (abandon) ------------------------------
# see tmp job
t_cell_list <- SplitObject(tcell_harm, split.by = 'orig.ident')
save(t_cell_list, file = './tmp/RISC_input_t_cell_list.rdata')

# NOTE: please see job 
# risc_list <- list()
# for(i in c(1:length(t_cell_list))){
#   data = t_cell_list[[i]]@assays$RNA@counts
#   coldat <- t_cell_list[[i]]@meta.data
#   rawdata <- data.frame(Symbol = rownames(data), RNA = "Gene Expression", row.names = rownames(data))
#   risc_list[[i]] <- readscdata(count = data, cell = coldat, gene = rawdata)
# }
# 
# process0 <- function(obj0){
#   # Filter cells and genes
#   obj0 = scFilter(obj0, min.UMI = 1000, max.UMI = 40000, min.gene = 200, min.cell = 3)
#   # Normalize the raw counts
#   obj0 = scNormalize(obj0)
#   # Find highly variable genes
#   obj0 = scDisperse(obj0)
#   # print(length(obj0@vargene))
#   return(obj0)
# }
# 
# options(repr.plot.width=6, repr.plot.height=5)
# for(i in c(1:length(risc_list))){
#   risc_list[[i]] <- process0(risc_list[[i]])
# }
# 
# var0 = risc_list[[1]]@rowdata$Symbol
# for(i in c(1:length(risc_list)) ){
#   var0 = intersect(var0,risc_list[[i]]@rowdata$Symbol )
# }
# 
# InPlot(risc_list,var.gene = var0, Std.cut = 0.99, ncore = 32, minPC = 16, nPC = 40)
# 
# data0 = scMultiIntegrate(
#   objects = risc_list, eigens = 18, add.Id = NULL, var.gene = var0, 
#   method = "RPCI", align = 'OLS', npc = 50, adjust = TRUE, 
#   ncore = 32, do.fast = "AUTO"
# )
# 
# data0 = scUMAP(data0, npc = 18, use = 'PLS')
# 
# data0 = scCluster(data0, slot = "cell.pls", neighbor = 3, method = "louvain", npc = 20)
# print(table(data0@coldata$Cluster))


