library(Seurat)
library(tidyverse)
library(iTALK)
library(data.table)


setwd('/data/sle')
source('./scripts/function_R/utils.R')
output_path <- 'output_file/seurat/all_pbmc/'
load('./output_file/seurat/all_pbmc/final_concat_pbmc_without_pSS.rdata')

PBMC_list <- SplitObject(pbmc_final,split.by = 'disease')
sle_pbmc <- PBMC_list$SLE
hc_pbmc <- PBMC_list$HC
rm(PBMC_list)

Idents(sle_pbmc) <- 'subtype'
Idents(hc_pbmc) <- 'subtype'
sle_pbmc_downsample <- subset(sle_pbmc, downsample =1000)
hc_pbmc_downsample <- subset(hc_pbmc, downsample =1000)
rm(sle_pbmc,hc_pbmc)

# # sle_iTALK <- seu_to_iTALK(seu_obj = sle_pbmc)
# # hc_iTALK <- seu_to_iTALK(hc_pbmc)
# back_run(func = seu_to_iTALK, out_name = 'sle_data_iTALK', job_name = 'sle_iTALK',
#          seu_obj = sle_pbmc_downsample)
# back_run(func = seu_to_iTALK, out_name = 'hc_data_iTALK', job_name = 'hc_iTALK',
#          seu_obj = hc_pbmc_downsample)
downsample_count_mat <- sparse_to_dense(dgc_matrix = sle_pbmc_downsample@assays$RNA@counts,formart = 'DF')
sle_iTalk_data <- transpose(downsample_count_mat)
colnames(sle_iTalk_data) <- rownames(downsample_count_mat)
rownames(sle_iTalk_data) <- colnames(downsample_count_mat)
sle_iTalk_data <- cbind(data.frame(cell_type = sle_pbmc_downsample$subtype), sle_iTalk_data)
sle_iTalk_data <- cbind(data.frame(compare_group = sle_pbmc_downsample$disease), sle_iTalk_data)

downsample_count_mat <- sparse_to_dense(dgc_matrix = hc_pbmc_downsample@assays$RNA@counts,formart = 'DF')
hc_iTalk_data <- transpose(downsample_count_mat)
colnames(hc_iTalk_data) <- rownames(downsample_count_mat)
rownames(hc_iTalk_data) <- colnames(downsample_count_mat)
hc_iTalk_data <- cbind(data.frame(cell_type = hc_pbmc_downsample$subtype), hc_iTalk_data)
hc_iTalk_data <- cbind(data.frame(compare_group = hc_pbmc_downsample$disease), hc_iTalk_data)

# highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
back_run(func = rawParse,out_name = 'hc_iTALK_hvg', job_name = 'hc_iTALK_hvg',
         input = hc_iTalk_data)
back_run(func = rawParse,out_name = 'sle_iTalk_hvg', job_name = 'sle_iTalk_hvg',
         input = sle_iTalk_data)


