library(Augur)
library(Suerat)
library(tidyverse)
source('./utils.R')
load('../data/modify_subtype_final_concat_pbmc_without_pSS.rdata')
# pbmc_final

expr_mat  <-  sparse_to_dense(pbmc_final@assays$RNA@counts,formart = 'matrix',threads = 20)

print('-----------------start augur----------------------')
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

pdf(file = './Augur_picsM.pdf')
print(plot_lollipop(augur_sub))
print(plot_lollipop(augur_main))
print(plot_umap(augur_sub, pbmc_final))
print(plot_umap(augur_sub, pbmc_final, mode = 'rank'))
print(plot_umap(augur_main, pbmc_final))
print(plot_umap(augur_main, pbmc_final, mode = 'rank'))
dev.off()





