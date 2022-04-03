
library(scater)
library(Seurat)

load('./output_file/seurat/all_pbmc/modify_subtype_final_concat_pbmc_without_pSS.rdata')
pbmc.sce <- as.SingleCellExperiment(pbmc_final)
pbmc.sce$cell <- group_by()
for(n in assayNames(pbmc.sce)) {
    tmp <- runPCA(pbmc.sce, exprs_values = n, ncomponents = 20)
    
    print(
        plotPCA(
            tmp,
            colour_by = "treatment",
            size_by = "nFeature_RNA",
            shape_by = "disease"
        ) +
            ggtitle(n)
    )
}

