################################################################################
# In this scripts we concat all of the cell subsets together to draw the whole 
# landscape of the PBMC, we will use below subtypes
# - T cell
#    - CD 4
#    - CD 8
#    - Other T cells
# - B cell
# - Plasma
# - Monocyte and DC
# - Platelet
# We do not include any PICs here
################################################################################

# -------------------------------- Initialize ----------------------------------
setwd('/data/sle')
source('./scripts/function_R/utils.R')
output_path <- 'output_file/seurat/all_pbmc/'

# read in the final version of the subtype annotation file 
load('./output_file/seurat/t_cell/final_t.cd4.rdata')  # t.cd4.filter
load('./output_file/seurat/t_cell/final_t.cd8.rdata')  # t.cd8.filter
load('./output_file/seurat/t_cell/tcell_harmony_anno.rdata') # tcell
# doublet exclusion
tcell_doublet_index <- read.csv('./scripts/scrublet/tcell_doublet_result.csv')
tcell$scrublet_doublet <- (data.frame(barcode = Cells(tcell)) %>% left_join(tcell_doublet_index))[,2]
Idents(tcell) <- 'scrublet_doublet'
tcell <- subset(tcell, idents = 'False')
Idents(tcell) <- 'subtype'
other_t <- subset(tcell, idents = c('MAIT','NK','NKT','T.prolife','T.yd'))

load('./output_file/seurat/mono_dc/final_mono_dc.rdata')
load('./output_file/seurat/b_cell/final_bcell.rdata')
load('./output_file/seurat/platelet/final_platelet.rdata')

# --------------------------- Filter and Harmony -------------------------------
pbmc_final <- merge(other_t,y = c(t.cd4.filter,t.cd8.filter,bcell_filter,
                                  mono_dc_filter,platelet_harmony),
                    merge.data = F,merge.dr = F)
# exclude the pSS_pah sample
Idents(pbmc_final) <- 'orig.ident'
pbmc_final <- subset(pbmc_final, idents = c('WH1','WH2'), invert = T)

# pbmc_final <- do_harmony(seu_obj = pbmc_final,harmony_slot = 'orig.ident',max.iter = 30,from_begin = T)
back_run(func = do_harmony, out_name = 'pbmc_final', job_name = 'pbmc_final',
         seu_obj = pbmc_final,harmony_slot = 'orig.ident',max.iter = 30,from_begin = T)

# -------------------------------- Analysis ------------------------------------
DimPlot(pbmc_final, group.by = 'subtype', label = T, raster = F,
        cols = get_color(length(table(pbmc_final$subtype))))

# -------------------------------- Modify meta ---------------------------------

pbmc_final$disease <- 'SLE'
pbmc_final$disease[which(pbmc_final$group == 'HC')] <- 'HC'

#---------------------------------- Save ---------------------------------------
save(pbmc_final, file = './output_file/seurat/all_pbmc/final_concat_pbmc_without_pSS.rdata')


