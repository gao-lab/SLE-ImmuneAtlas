library(Seurat)
library(tidyverse)
library(harmony)

setwd('/data/sle')
source('./scripts/function_R/utils.R')

#---------------------------------- Load ---------------------------------------
load('./other_sc_data/data/17_NBT_8_SLE_patients_Demuxlet/nbt_17_meta.rdata')
load('./other_sc_data/data/19_PNAS_neutrophil_MultiOmics/pnas_19_meta.rdata')
load('./other_sc_data/data/21_EBioMedicine_IFN/ebi_21_meta.rdata')
load('./other_sc_data/data/19_NC_1M_eQTL/nc_19_sle_flare.rdata')
load('./other_sc_data/data/19_NC_1M_eQTL/nc_19_sle.rdata')
load('./other_sc_data/data/20_NatImm_children/nat_imm_anno_all_cell.Rdata')
natim_20 <- all_cell
rm(all_cell)

natim_20 <- CreateSeuratObject(natim_20@assays$RNA@counts,meta.data = natim_20@meta.data, min.cells = 10)
nc_sle <- CreateSeuratObject(sle@assays$RNA@counts, meta.data = sle@meta.data, min.cells = 10)
nc_sle.flare <- CreateSeuratObject(sle_flare@assays$RNA@counts, meta.data = sle_flare@meta.data, min.cells = 10)
rm(sle,sle_flare)
gc()

#------------------------------  batch by sample ------------------------------
ebi_21$dataset <- 'ebi_21'
natim_20$dataset <- 'natim_20'
nbt_17$dataset <- 'nbt_17'
nc_sle$dataset <- 'nc_sle'
nc_sle.flare$dataset <- 'nc_sle.flare'
pnas_19$dataset <- 'pnas_19'
pbmc_all$dataset <- 'pbmc_all'
collect_cell <- merge(natim_20,y = c(ebi_21,nbt_17,nc_sle,nc_sle.flare,pnas_19,pbmc_all))

save(collect_cell, file = 'other_sc_data/merge/pbmc_collect_meta_raw.rdata')
load( file = 'other_sc_data/merge/pbmc_collect_meta_raw.rdata')

# back_run(func = do_harmony, out_name = 'collect_cell', job_name = 'collect_cell',
        # seu_obj = collect_cell,harmony_slot = 'orig.ident',max.iter = 30,from_begin = T)
do_harmony(seu_obj = collect_cell,harmony_slot = 'orig.ident',max.iter = 30,
           from_begin = T,res = c(1,1.2,0.8))

save(collect_cell, file = 'other_sc_data/merge/pbmc_collect_meta_harm.rdata')

DimPlot(collect_cell)

#------------------------------ Transfer label ---------------------------------
load('./final/seurat/pbmc/03-pbmc_all_concat.rdata')

ebi_21.anchors <- FindTransferAnchors(reference = pbmc_all, query = ebi_21,
                                        dims = 1:30, reference.reduction = "pca")
ebi_21.predictions <- TransferData(anchorset = ebi_21.anchors, refdata = pbmc_all$subtype,
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

trans_labbel <- function(reference, query){
  anchors <- FindTransferAnchors(reference = reference, query = ebi_21,
                                 dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = reference$subtype,
                              dims = 1:30)
  return(predictions)
}


#------------------------------ Transfer UMAP ---------------------------------








