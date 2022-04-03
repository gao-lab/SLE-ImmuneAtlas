
# Popalign is bulid on python so need do some preparation to import the data
library(DropletUtils)
setwd('/data/sle')

write10xCounts(path ='./data/10x_mtx_format/' ,x = pbmc_final@assays$RNA@counts,
               version = '3')
back_run(func =write10xCounts , out_name = 'tmp',job_name = 'wtite_10x',
         path ='./data/10x_mtx_format/' ,x = pbmc_final@assays$RNA@counts,
         version = '3')
meta <- pbmc_final@meta.data
meta$cell_barcode <- rownames(meta)
meta$sample_id <- meta$orig.ident

write.csv(meta, file = './data/10x_mtx_format/meta.csv')

# -------------------------- Set Python Env ------------------------------------
library(reticulate)
use_miniconda(condaenv = '/root/miniconda3/envs/Popalign/' ,required = T)
