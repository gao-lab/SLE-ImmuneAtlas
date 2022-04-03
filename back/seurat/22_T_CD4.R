
library(cowplot)
library(harmony)
library(RISC)
library(future)

setwd('/data/sle')
output_path <- './output_file/seurat/t_cell/'
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./output_file/seurat/t_cell/CD4_tcell_filter_harmony_01.rdata')
DimPlot(t.cd4.filter, label = T) + DimPlot(t.cd4.filter,group.by = 'subtype')
FeaturePlot(t.cd4.filter,features = c('CD4','CD8A','CD8B'), ncol = 3)

t.cd4.filter <- subset(t.cd4.filter, idents = 11, invert = T)

# add meta
t.cd4.filter$disease <- 'SLE'
t.cd4.filter$disease[which(t.cd4.filter$group == 'HC')] <- 'HC'


DimPlot(t.cd4.filter, label = T, group.by = 'subtype',cols = getPalette(9), pt.size = 0.3) + NoAxes() 
  # title('CD4+ T cell')
plot_scdata(t.cd4.filter, color_by = "treatment", pal_setup = 'Set2') + 
  plot_scdata(t.cd4.filter, color_by = "seurat_clusters", pal_setup = 'Set1')

#---------------------------------- projectTIL anno ----------------------------
library(ProjecTILs)
library(scRepertoire)
ref <- load.reference.map('./source/projectTIL/ref_TILAtlas_mouse_v1.rds')


down_sample <- T
size <- 1000
ncores= 4
for(i in c(0:1)){
  tmp <- subset(t.cd4.filter, idents = i)
  if(down_sample){tmp <- tmp[, sample(colnames(tmp), size =size, replace=F)]}
  out_name =  paste0('project.cd4.cluster_',i)
  print(out_name)
  back_run(func = make.projection,out_name = out_name,job_name = out_name,
           tmp = tmp, ref=ref, filter.cells = F,future.maxSize=1024*20, ncores = ncores)
}

# sample_data <- t.cd4.filter[, sample(colnames(t.cd4.filter), size =1000, replace=F)]
sample_data <- subset(t.cd4.filter, idents =4)
# querydata <- ProjecTILs::query_example_seurat
query.projected <- make.projection(sample_data, ref=ref,filter.cells = F,ncores = 32, future.maxSize = 1024*200)
plot.projection(ref, query.projected)
query.projected <- cellstate.predict(ref=ref, query=query.projected)
plot.statepred.composition(ref, query.projected,metric = "Percent")
plot.states.radar(ref, query=query.projected, min.cells=30)

#--------------------------------- Vis Marker ----------------------------------
VlnPlot(t.cd4.filter, features = c('CD4','CD8A','CD8B','PDCD1'))
DotPlot2(t.cd4.filter, marker_list = t_sub_marker)
DotPlot2(t.cd4.filter, marker_list = t_cd4_all_marker)
DotPlot2(t.cd4.filter, marker_list = c('IFIT3','ISG15','MX1','OAS1','IFI44L'))
DotPlot(t.cd4.filter, features = marker_list)

DotPlot2(t.cd4.filter, marker_list = IL_genes)
DotPlot2(t.cd4.filter, marker_list = c('LRRN3','CCR7','NPM1','SELL','PASK','IL7R',
                                       'GZMH','GZMA','GNLY','IFIT2','ISG15','IFI27',
                                       'FOXP3'))
DotPlot2(t.cd4.filter, marker_list = c('CXCR3','CCL5','GZMK'))

#------------------------------- Find Markers ----------------------------------
plan('multiprocess', workers = 24)
options(future.globals.maxSize= 250 * 1000 * 1024 ^2 ) 
marker_all_cd4 <- FindAllMarkers(t.cd4.filter, only.pos = T, logfc.threshold = 0.25)
marker_all_cd4.filter <-  marker_all_cd4 %>% filter(p_val_adj <0.05)
marker_all_cd4.filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10) %>% View()

VlnPlot(t.cd4.filter,features = c('PDCD1','ICOS'), split.by = 'disease',stack = T)
#---------------------------------- Anno ---------------------------------------
t.cd4.filter$subtype <- 'unknown'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(8))] <- 'T0_naive_ifn'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(6,10))] <- 'Treg'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(3))] <- 'Th1_Temra'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(7))] <- 'Th2'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(9))] <- 'Th17'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(2))] <- 'Tcm'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(4))] <- 'Tem'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(5))] <- 'Tfh'
# t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(9))] <- 'Tfh'
t.cd4.filter$subtype[which(t.cd4.filter$seurat_clusters %in% c(0,1))] <- 'Th0_naive'
DimPlot(t.cd4.filter, label = T, group.by = 'subtype')

Idents(t.cd4.filter) <- 'subtype'

#---------------------------------- Ratio---------------------------------------
colourCount = length(11)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# overall plasma ratio
t.cd4.filter@meta.data %>% filter(group != 'pSS_pah') %>% 
  ggplot( aes(x = disease , fill = seurat_clusters))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(t.cd4.filter$seurat_clusters))) + 
  scale_fill_manual(values =getPalette(11)) +
  labs(fill="group") + coord_flip() 

# By Group
t.cd4.filter@meta.data %>% filter(group !='pSS_pah') %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t.cd4.filter@meta.data[,c(1,12,13,14,19)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='disease',y='Ratio', color = 'disease',
                    palette =c("orange", "green"))+ 
  facet_wrap(~subtype,scales = "free") + 
  stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )


#---------------------------------- Save ---------------------------------------
save(t.cd4.filter,file = './output_file/seurat/t_cell/CD4_tcell_filter_harmony_anno.rdata')
save(t.cd4.filter,file = './output_file/seurat/t_cell/final_t.cd4.rdata')
