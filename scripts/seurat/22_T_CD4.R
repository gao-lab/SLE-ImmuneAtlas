
library(cowplot)
library(harmony)
library(RColorBrewer)
library(future)

setwd('/data/sle')
source('./scripts/function_R/utils.R')

# Already Harmony 
load('./final/seurat/t_cell/03-CD4_Tcell_raw_harm.rdata')
DimPlot(cd4_filter, label = T) + DimPlot(cd4_filter,group.by = 'subtype')
FeaturePlot(cd4_filter,features = c('CD4','CD8A','CD8B'), ncol = 3)

cd4_filter <- subset(cd4_filter, idents = 11, invert = T)
cd4_filter_index <-CellSelector(DimPlot(cd4_filter))
cd4_filter <- subset(cd4_filter, invert=  T, cells = cd4_filter_index)


DimPlot(cd4_filter, label = T,cols = rev(get_color(11)), group.by = 'subtype', pt.size = 0.3) + NoAxes() 
  # title('CD4+ T cell')
# plot_scdata(cd4_filter, color_by = "treatment", pal_setup = 'Set2') + 
#   plot_scdata(cd4_filter, color_by = "seurat_clusters", pal_setup = 'Set1')

FeaturePlot(cd4_filter, features = c('IL17A','IFNG','IL17F','IL17B','IL1B'))

#--------------------------------- Vis Marker ----------------------------------
VlnPlot(cd4_filter, features = c('CD4','CD8A','CD8B','PDCD1'))
DotPlot2(cd4_filter, marker_list = t_sub_marker)
DotPlot2(cd4_filter, marker_list = t_cd4_all_marker)
DotPlot2(cd4_filter, marker_list = c('IFIT3','ISG15','MX1','OAS1','IFI44L'))
DotPlot2(cd4_filter, marker_list = marker_list)

DotPlot2(cd4_filter, marker_list = IL_genes)
DotPlot2(cd4_filter, marker_list = c('LRRN3','CCR7','NPM1','SELL','PASK','IL7R',
                                       'GZMH','GZMA','GNLY','IFIT2','ISG15','IFI27',
                                       'FOXP3'))
DotPlot2(cd4_filter, marker_list = c('CXCR3','CCL5','GZMK'))

#------------------------------- Find Markers ----------------------------------
plan('multiprocess', workers = 16)
options(future.globals.maxSize= 250 * 1000 * 1024 ^2 ) 
marker_all_cd4 <- FindAllMarkers(cd4_filter, only.pos = T, logfc.threshold = 0.25)
marker_all_cd4.filter <-  marker_all_cd4 %>% filter(p_val_adj <0.05)
marker_all_cd4.filter %>% group_by(cluster) %>% top_n(avg_log2FC  , n = 10) %>% View()

marker_cd4_c1 <- FindMarkers(cd4_filter, ident.1 = 1, only.pos = T,logfc.threshold = 0.25)

marker_cd4_c10 <- FindMarkers(cd4_filter, ident.1 = 10, only.pos = T,logfc.threshold = 0.25) %>% 
  filter(p_val_adj <0.05) %>% arrange(-avg_log2FC )

VlnPlot(cd4_filter,features = c('PDCD1','ICOS'), split.by = 'group',stack = T)
#---------------------------------- Anno ---------------------------------------
cd4_filter$subtype <- 'unknown'

# we find that we can not explain the Th1 ratio up with Th17 ratio down in SLE
# So I decided to anno them as Th1-liked cell totally
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(10))] <- 'T.CD4.IFN-response'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(7,9))] <- 'T.CD4.Treg'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(3))] <- 'T.CD4.Th1.cxcr3'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(5))] <- 'T.CD4.Th1'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(6))] <- 'T.CD4.Th2'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(4))] <- 'T.CD4.Th17'
# # cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(2))] <- 'Tcm'
# # cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(4))] <- 'Tem'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(8))] <- 'T.CD4.Tfh'
# # cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(9))] <- 'Tfh'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(0,1,2))] <- 'T.CD4.naive'
# DimPlot(cd4_filter, label = T, group.by = 'subtype')


cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(10))] <- 'T.CD4.IFN-response'
cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(7,9))] <- 'T.CD4.Treg'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(3))] <- 'T.CD4.Th1.cxcr3'
cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(3,4,5))] <- 'T.CD4.Th1'
cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(6))] <- 'T.CD4.Th2'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(4))] <- 'T.CD4.Th17'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(2))] <- 'Tcm'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(4))] <- 'Tem'
cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(8))] <- 'T.CD4.Tfh'
# cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(9))] <- 'Tfh'
cd4_filter$subtype[which(cd4_filter$seurat_clusters %in% c(0,1,2))] <- 'T.CD4.naive'
DimPlot(cd4_filter, label = T, group.by = 'subtype')

Idents(cd4_filter) <- 'subtype'

#---------------------------------- Ratio---------------------------------------
colourCount = length(11)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

# overall  ratio
cd4_filter@meta.data %>% filter(group != 'pSS_pah') %>% 
  ggplot( aes(x = group , fill = seurat_clusters))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(cd4_filter$seurat_clusters))) + 
  scale_fill_manual(values =getPalette(11)) +
  labs(fill="group") + coord_flip() 

# By subtype
cd4_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  filter(group != 'treated') %>% 
  ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                    palette =c("#FC4E07" ,"#00AFBB"))+ 
  facet_wrap(~subtype,scales = "free",ncol = 8) + 
  stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )

# By cluster
cd4_filter@meta.data  %>% 
  group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
  ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                    palette =c("#FC4E07","#00AFBB"))+ 
  facet_wrap(~seurat_clusters,scales = "free") + 
  stat_compare_means(  method = "t.test" )

# treatment paired
tmp1<- cd4_filter@meta.data  %>% 
  group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')  %>% filter(!subtype =='unknown') 
tmp2<- cd4_filter@meta.data  %>% 
  group_by(orig.ident,subtype,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')%>% filter(!subtype =='unknown')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>%
  # mutate(across(subtype,factor, levels = c("B.transition","B.naive","B.IFN-response",
  #                                          "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
  ggpaired( cond1 = 'before', cond2 = 'after',
            fill  = "condition", line.color = "gray", line.size = 0.4,
            palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,label = 'p.format')+
  ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 8 )

# treatment paired (cluster)
tmp1<- cd4_filter@meta.data  %>% 
  group_by(orig.ident,seurat_clusters,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
  filter(!orig.ident == 'XYY2')  %>% filter(!seurat_clusters =='unknown') 
tmp2<- cd4_filter@meta.data  %>% 
  group_by(orig.ident,seurat_clusters,.drop = FALSE) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
  filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
  filter(!orig.ident == 'XYY2')%>% filter(!seurat_clusters =='unknown')
data.frame(sample = tmp1$orig.ident, seurat_clusters= tmp1$seurat_clusters, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>%
  # mutate(across(seurat_clusters,factor, levels = c("B.transition","B.naive","B.IFN-response",
  #                                          "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-"))) %>%
  ggpaired( cond1 = 'before', cond2 = 'after',
            fill  = "condition", line.color = "gray", line.size = 0.4,
            palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4,label = 'p.format')+
  ylab('Prolife Plasma ratio') + facet_wrap(~seurat_clusters,scales= "free",ncol = 4 )

#---------------------------------- Save ---------------------------------------
save(cd4_filter,file = './final/seurat/t_cell/04-CD4_Tcell_filter_anno.rdata')



#----------------------------- projectTIL anno (abandon)---------------------------
library(ProjecTILs)
ref <- load.reference.map('./source/projectTIL/ref_TILAtlas_mouse_v1.rds')


down_sample <- T
size <- 1000
ncores= 4
for(i in c(0:1)){
  tmp <- subset(cd4_filter, idents = i)
  if(down_sample){tmp <- tmp[, sample(colnames(tmp), size =size, replace=F)]}
  out_name =  paste0('project.cd4.cluster_',i)
  print(out_name)
  back_run(func = make.projection,out_name = out_name,job_name = out_name,
           tmp = tmp, ref=ref, filter.cells = F,future.maxSize=1024*20, ncores = ncores)
}

# sample_data <- cd4_filter[, sample(colnames(cd4_filter), size =1000, replace=F)]
sample_data <- subset(cd4_filter, idents =4)
# querydata <- ProjecTILs::query_example_seurat
query.projected <- make.projection(sample_data, ref=ref,filter.cells = F,ncores = 32, future.maxSize = 1024*200)
plot.projection(ref, query.projected)
query.projected <- cellstate.predict(ref=ref, query=query.projected)
plot.statepred.composition(ref, query.projected,metric = "Percent")
plot.states.radar(ref, query=query.projected, min.cells=30)



