library(Seurat)
library(tidyverse)

#---------------------- all maintype change in the pbmc ------------------------
# By subtype
pbmc_all@meta.data  %>% 
    group_by(orig.ident,main_type) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07" ,"#00AFBB"))+ 
    facet_wrap(~main_type,scales = "free",ncol = 4) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )

# add  cd4 and cd8 label
pbmc_all$cd4_8 <- 'not'
pbmc_all$cd4_8[which( grepl('CD4',pbmc_all$subtype))] <- 'CD4'
pbmc_all$cd4_8data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg==[which( grepl('CD8',pbmc_all$subtype))] <- 'CD8'

pbmc_all@meta.data  %>% 
    group_by(orig.ident,cd4_8) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#FC4E07" ,"#00AFBB"))+ 
    facet_wrap(~cd4_8,scales = "free",ncol = 4) + 
    stat_compare_means(comparisons = list(c("SLE", "HC")),  method = "t.test" )


#------------------------ Mono.CD14.recruit function ---------------------------
# cluster == 4 means subset the Mono.CD14.recruit marker only
marker_mono_recruit <- filter(marker_all_mono.dc.filter, cluster ==4)
marker_mono_recruit %<>% arrange(-avg_log2FC )
keggo_plot(marker_mono_recruit$gene)
plot_reactome(marker_mono_recruit)

# P <- plot_gsea(mono_dc_filter,group_by = 'subtype', focus = 'Mono.CD14.recruit', 
#                title = 'GSEA',category = 'C5')

# test gsva(slow: ~20min in 2500 cells)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
# library(org.Mm.eg.db)
library(dplyr)
library(msigdbr)

mat <- subset(mono_dc_filter, idents = 'Mono.CD14.recruit')
mat <- as.matrix(mat@assays$RNA@counts)
meta <- mono_dc_filter@meta.data[,c("orig.ident", "group", "treatment")]
species = 'Homo sapiens'; category = 'C5' ; subcategory =NULL

gene_set <- msigdbr(species = species, category = category, 
                    subcategory = subcategory) %>% split(x = .$gene_symbol, f = .$gs_name)
gsva_res <- gsva(expr = mat, 
               gset.idx.list = gene_set,
               kcdf="Poisson",
               parallel.sz = 24)
back_run(gsva,out_name = 'gsva_res',job_name = 'gsva_res',
         expr = mat, gset.idx.list = gene_set,kcdf="Poisson",parallel.sz = 24)










