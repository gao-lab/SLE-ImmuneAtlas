library(sceasy)
library(Seurat)
# use_condaenv('sceasy')

load(snakemake@input[["seu_obj"]])
print(ls()[1])

tmp <- get(ls()[1])

tmp@assays$RNA@data <- tmp@assays$RNA@counts
convertFormat(tmp, from="seurat", to="anndata",
                      outFile=snakemake@output[["h5ad"]])