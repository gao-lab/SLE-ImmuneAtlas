library(Augur)
library(ggplot2)
library(ggpubr)

setwd('/data/sle')
source('./scripts/function_R/utils.R')
load('./output_file/seurat/all_pbmc/downsample_final_concat_pbmc_without_pSS.rdata')
load('./output_file/seurat/all_pbmc/final_concat_pbmc_without_pSS.rdata')


pbmc_final$subtype[which(pbmc_final$subtype ==  c('T0_naive_ifn'))] <- 'CD4_T0_naive_ifn'
pbmc_final$subtype[which(pbmc_final$subtype == c('Treg'))] <- 'CD4_Treg'
pbmc_final$subtype[which(pbmc_final$subtype == c('Th1_Temra'))] <- 'CD4_Th1_Temra'
pbmc_final$subtype[which(pbmc_final$subtype == c('Th2'))] <- 'CD4_Th2'
pbmc_final$subtype[which(pbmc_final$subtype == c('Th17'))] <- 'CD4_Th17'
pbmc_final$subtype[which(pbmc_final$subtype == c('Tcm'))] <- 'CD4_Tcm'
pbmc_final$subtype[which(pbmc_final$subtype == c('Tem'))] <- 'CD4_Tem'
pbmc_final$subtype[which(pbmc_final$subtype == c('Tfh'))] <- 'CD4_Tfh'
# pbmc_final$subtype[which(pbmc_final$subtype == c('Tfh'))] <- 'Tfh'
pbmc_final$subtype[which(pbmc_final$subtype == c('Th0_naive'))] <- 'CD4_Th0_naive'


# add the main meta to make data comparable with IFN-stimulated PBMC data
pbmc_final$subtype %>% table()
pbmc_final$compare_meta <- pbmc_final$subtype 
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = '^B')] <- 'B cells'
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = '^CD8')] <- 'CD8 T cells'
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = '^CD4')] <- 'CD4 T cells'
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = 'NK')] <- 'NK cells'
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = 'Platelet')] <- 'Megakaryocytes'
pbmc_final$compare_meta[str_detect(pbmc_final$subtype,pattern = 'Platelet')] <- 'Megakaryocytes'

pbmc_final$compare_meta %>% table()

save(pbmc_final,file = './output_file/seurat/all_pbmc/modify_subtype_final_concat_pbmc_without_pSS.rdata')


augur = calculate_auc(input = pbmc_final_downsample@assays$RNA@counts %>% as.matrix(),
                      meta = pbmc_final_downsample@meta.data,
                      label_col = 'disease',
                      cell_type_col = 'subtype', n_threads=16)
# back_run(func = calculate_auc,out_name = 'augur',job_name ='augur' ,
#          input = pbmc_final_downsample@assays$RNA@counts %>% as.matrix(),
#          meta = pbmc_final_downsample@meta.data,
#          label_col = 'disease',
#          cell_type_col = 'subtype', n_threads=16)

save(augur,file = './output_file/Augur/augur_downsample_pbmc.rdata')
plot_lollipop(augur)

plot_umap(augur, pbmc_final_downsample)
plot_umap(augur, pbmc_final_downsample, mode = 'rank')
plot_scatterplot()


# ----------------------------- Run in HPC -------------------------------------

augur_sub = calculate_auc(input = expr_mat,
                      meta = pbmc_final@meta.data,
                      label_col = 'disease',
                      cell_type_col = 'subtype', n_threads=24)
save(augur_sub, file = './augur_sub.rdata')

augur_main = calculate_auc(input = expr_mat,
                           meta = pbmc_final@meta.data,
                           label_col = 'disease',
                           cell_type_col = 'compare_meta', n_threads=24)
save(augur_main, file = './augur_main.rdata')

pdf(file = './Augur_pics.pdf')
print(plot_lollipop(augur_sub))
print(plot_lollipop(augur_main))
print(plot_umap(augur_sub, pbmc_final))
print(plot_umap(augur_sub, pbmc_final, mode = 'rank'))
print(plot_umap(augur_main, pbmc_final))
print(plot_umap(augur_main, pbmc_final, mode = 'rank'))
dev.off()


