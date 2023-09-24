library(Seurat)
library(stringr)
library(future)
library(celldex)
library(stringr)
library(SingleR)
library(BiocParallel)
library(pheatmap)
library(tidyverse)
library(feather) 
library(ggpubr)
setwd('/data/sle')
source('./scripts/function_R/utils.R')

output_path <- './addtional/seurat/pbmc/'

#---------------------------- Read in Files-------------------------------------
samples_list <- list.files('./data/10x_rna/',pattern = "h5$")
add_samples_list <- list.files('./data/addtional/10x_rna/',pattern = "h5$")

meta <- read.csv('vdj/immcatation/all_sample_meta.csv')
length(samples_list)
samples_list

seu_list <- c()
for (i in samples_list) {
  print(i)
  tmp_seu_obj <- CreateSeuratObject(Read10X_h5(paste0('./data/10x_rna/',i)),
                                    min.cells = 3, min.features = 200, 
                                    project = str_split(i,"_")[[1]][1])
  seu_list <- append(seu_list, tmp_seu_obj)
  rm(tmp_seu_obj)
}
for (i in add_samples_list) {
    print(i)
    tmp_seu_obj <- CreateSeuratObject(Read10X_h5(paste0('./data/addtional/10x_rna/',i)),
                                      min.cells = 3, min.features = 200, 
                                      project = str_split(i,"_")[[1]][1])
    seu_list <- append(seu_list, tmp_seu_obj)
    rm(tmp_seu_obj)
}

all_cell <- merge(seu_list[[1]], c(seu_list[[2]],seu_list[[3]],seu_list[[4]],
                                   seu_list[[5]],seu_list[[6]],seu_list[[7]],
                                   seu_list[[8]],seu_list[[9]],seu_list[[10]],
                                   seu_list[[11]],seu_list[[12]],seu_list[[13]],
                                   seu_list[[14]],seu_list[[15]],seu_list[[16]],
                                   seu_list[[17]],seu_list[[18]],seu_list[[19]],
                                   seu_list[[20]],seu_list[[21]],seu_list[[22]],
                                   seu_list[[23]],seu_list[[24]],seu_list[[25]],
                                   seu_list[[26]]
                                   ))
rm(seu_list)


#------------------------------ Add Meta ---------------------------------------
# add group meta
all_cell$group <- 'SLE' 
all_cell$group[which(all_cell$orig.ident %in% c('ZMY1','QJY','ZS','ZH'))] <- 'HC'
# all_cell$group[which(all_cell$orig.ident %in% 
#                        c('XH','LL','XYY','XYY2','LL2'))] <- 'SLE_pah'

# add treatment(5 treated patients with 6 samples: XYY&XYY2 --> XH)
all_cell$treatment <- 'untreated'
all_cell$treatment[which(all_cell$orig.ident %in% c('ZMY1','QJY','ZS','ZH'))] <- 'HC'
all_cell$treatment[which(all_cell$orig.ident %in% 
                           c('XYY','XYY2','LL2','ZPP2','WYF2','HXR2'))] <- 'treated'

# add paired meta
all_cell$pair <- 'unpaired' 
all_cell$pair[which(all_cell$orig.ident %in% c('XH','XYY','XYY2'))] <- 'XH_pair'
all_cell$pair[which(all_cell$orig.ident %in% c('LL','LL2'))] <- 'LL_pair'
all_cell$pair[which(all_cell$orig.ident %in% c('WH','WH2'))] <- 'WH_pair'
all_cell$pair[which(all_cell$orig.ident %in% c('ZPP','ZPP2'))] <- 'ZPP_pair'
all_cell$pair[which(all_cell$orig.ident %in% c('WYF','WYF2'))] <- 'WYF_pair'
all_cell$pair[which(all_cell$orig.ident %in% c('HXR','HXR2'))] <- 'HXR_pair'

# system.time(sceasy::convertFormat(all_cell, from="seurat", to="anndata",
#                                   outFile='./final/seurat/pbmc/all_cell_seu_without_filter.h5ad'))


#---------------------------- Quality Control ----------------------------------
all_cell <- PercentageFeatureSet(all_cell, "^MT-", col.name = "percent_mito")
all_cell <- PercentageFeatureSet(all_cell, "^RP[SL]", col.name = "percent_ribo")
# VlnPlot(all_cell, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend()

# add HBB gene filter
HBB_ratio <- all_cell@assays$RNA@counts['HBB',]/all_cell$nCount_RNA
filter_HBB_barcode <- Cells(all_cell)[HBB_ratio > 0.6]
all_cell <- subset(all_cell, cells = filter_HBB_barcode, invert = T )

all_cell <- subset(all_cell, subset = nFeature_RNA > 200 
                   & nFeature_RNA < 6000 
                   & percent_mito < 20
                   & nCount_RNA > 500
                   & nCount_RNA < 50000)
table(all_cell$orig.ident)
# VlnPlot(all_cell, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 4) + NoLegend()

# system.time(SaveH5Seurat(all_cell, filename = 'all_cell_seu_filter.h5Seurat'))


#---------------------------- Seurat Pipeline ----------------------------------
all_cell <- do_seurat(all_cell, res = c(0.8, 1.0))
# back_run(do_seurat, job_name = 'pbmc', out_name = 'all_cell',
#          seu_obj = all_cell, res = c(0.8, 1.0))
all_cell <- do_harmony(all_cell)

DimPlot(all_cell, label = T) + NoLegend()
DimPlot(all_cell, group.by = 'Phase')
DimPlot(all_cell, group.by = c('orig.ident', 'seurat_clusters'), label = T)
FeaturePlot(all_cell, features = feats, ncol = 2)

# Barplot prolife
ggplot(all_cell@meta.data, aes(x=all_cell$seurat_clusters, fill=all_cell$Phase))+
  geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
  scale_fill_discrete(names(table(all_cell$seurat_clusters))) + 
  labs(fill="Phase") + 
  coord_flip()

# vis all marker genes
VlnPlot(all_cell,features = c(marker_list,feats), stack = T)

#---------------------------- marker gene  ----------------------------------
marker_all_c16.c12 <- FindMarkers(all_cell, ident.1 = 16, ident.2 = 12, min.pct = 0.2, only.pos = T)

#---------------------------- marker gene  ----------------------------------
all_cell$main_type <- 'unknown'
# typical cell type
all_cell$main_type[which(all_cell$seurat_clusters %in% c(0,3,11,14,19))] <- 'T.cyto'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(2,4,7))] <- 'T.naive'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(15))] <- 'T.prolife'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(1,22))] <- 'Mono.CD14'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(12))] <- 'Mono.CD16'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(6,9))] <- 'Bcell'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(17,20))] <- 'Plasma'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(5))] <- 'NK'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(24))] <- 'pDC'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(21))] <- 'mDC'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(10))] <- 'Platelet'
# PIC/doublet
all_cell$main_type[which(all_cell$seurat_clusters %in% c(8, 18, 26))] <- 'Low_quality'
all_cell$main_type[which(all_cell$seurat_clusters %in% c(13, 16, 20, 23, 25, 27))] <- 'doublets'
table(all_cell$main_type)

DimPlot(all_cell,label = T, group.by = 'main_type') + NoLegend()
VlnPlot(all_cell, features = feats, ncol = 2)

#------------------------------- Save Files ------------------------------------
# less than 1 min
save(all_cell, file = './final/addtional/seurat/pbmc/02-pbmc_filter_cluster_main_type.rdata')
write.csv(all_cell@meta.data, file ='./final/addtional/seurat/pbmc/02-pbmc_filter_cluster_main_type.meta.csv')

# subset and save cell types 
Idents(all_cell) <- 'main_type'
tcell <- subset(all_cell,idents = c('T.cyto','T.naive','T.prolife','NK'))
mono_dc <- subset(all_cell, idents = c('Mono.CD14','Mono.CD16','pDC','mDC'))
bcell <- subset(all_cell,idents = c('Bcell','Plasma'))

# save
save(tcell, file = './final/addtional/seurat/t_cell/01-t_cell_raw.rdata')
save(mono_dc, file = './final/addtional/seurat/mono_dc/01-mono_dc_raw.rdata')
save(bcell, file = './final/addtional/seurat/b_cell/01-b_cell_raw.rdata')


#----------------------------- Basic Statistic ----------------------------------
tmp_df <- all_cell@meta.data %>% filter(group != 'pSS_pah') %>% 
  filter(group != 'pSS_pah') %>% filter(!str_detect(main_type, '_'))
ggplot(data = tmp_df, aes(x = main_type, fill = group))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tmp_df$group))) + 
  scale_fill_manual(values =c('green','red','orange')) +
  labs(fill="group") + coord_flip()

ggplot(data = tmp_df, aes(x = main_type, fill = treatment))+
  geom_bar(stat = 'count',position = 'fill')+
  labs(y = 'proportions', x = "") + 
  scale_fill_discrete(labels= names(table(tmp_df$group))) + 
  scale_fill_manual(values =c('green','orange','red')) +
  labs(fill="group") + coord_flip()

all_cell@meta.data %>% filter(PIC == 'no') %>% filter(main_type != 'Low_quality')%>% 
  group_by(orig.ident, main_type,PIC) %>% summarise(num = n()) %>% 
  left_join(table(all_cell$orig.ident) %>% as.data.frame(),by = c('orig.ident' = 'Var1')) %>%
  mutate(ratio = num/Freq * 100) %>% left_join(meta, c('orig.ident'='name')) %>%
  ggboxplot('main_type','ratio',fill = 'group') + xlab('')+ ylab('percentage') +
  rotate_x_text(angle = 45, hjust = 0.9, vjust = NULL) +
  stat_compare_means(aes(group = group), label = "p.signif",method = 't.test') 





