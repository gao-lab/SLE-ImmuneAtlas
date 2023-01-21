# highlight the cell with TCR and BCR

load('./final//scRepertoire/BCR/combined_bcr.rdata') # combined_bcr
load('./final//scRepertoire/TCR/combined_tcr.rdata') # combined_tcr


# -------------------------- for B cell ----------------------------------------
# need to unify the barcode in seurat and scPepertoire
bcr_sample <- list.files('./data/10x_bcr/',pattern = 'csv$')
bcr_sample_name <- str_split_fixed(bcr_sample,pattern = '_',2)[,1]

b_vdj_barcode <- pbmc_all@meta.data %>% rownames_to_column('barcode_sc')%>% 
  mutate(new_barcode = str_split_fixed(barcode_sc ,'_',2)[,1])%>% 
  mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 
pbmc_all$b_vdj_barcode <- b_vdj_barcode$scRepertoire
pbmc_all$vdj <- 'lack'

for (sample in unique(paste0(pbmc_all$orig.ident,'_',pbmc_all$group))) {
  tmp <- combined_bcr[[sample]][['barcode']]
  pbmc_all$vdj[which(pbmc_all$b_vdj_barcode %in% tmp)] <- 'BCR'
}

pbmc_all$vdj %>% table()

# -------------------------- for T cell ----------------------------------------
# tcr and bcr share same barcode 

for (sample in unique(paste0(pbmc_all$orig.ident,'_',pbmc_all$group))) {
  tmp <- combined_tcr[[sample]][['barcode']]
  pbmc_all$vdj[which(pbmc_all$b_vdj_barcode %in% tmp)] <- 'TCR'
}

pbmc_all$vdj %>% table()

# -------------------------------- Plot ----------------------------------------
pbmc_all$vdj <- factor(pbmc_all$vdj, levels =c('BCR','TCR','lack'))
DimPlot(pbmc_all, group.by = 'vdj', cols = c('red','blue','grey'),raster=FALSE) + NoAxes()


