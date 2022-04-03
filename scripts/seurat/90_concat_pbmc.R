################################################################################
# In this scripts we concat all of the cell subsets together to draw the whole 
# landscape of the PBMC, we will use below subtypes
# - T cell
#    - CD 4
#    - CD 8
#    - Other T cells
#    - NK
# - B cell
# - Plasma
# - Monocyte and DC
# - Platelet
# We do not include any PICs here
################################################################################

# -------------------------------- Initialize ----------------------------------
setwd('/data/sle')
source('./scripts/function_R/utils.R')


# read in the final version of the subtype annotation file 
load('./final/seurat/t_cell/04-t_cell_merge.rdata')  # t_cell_merge
load('./final/seurat/t_cell/04-CD4_Tcell_filter_anno.rdata')  # cd4_filter
load('./final/seurat/t_cell/04-CD8_Tcell_filter_anno.rdata')  # cd8_filter
load('./final/seurat/t_cell/04-Other_Tcell_filter_anno.rdata') # other_T_filter
load('./final/seurat/t_cell/04-NK_Tcell_filter_anno.rdata') # nk_filter

load('./final/seurat/mono_dc/03-mono_dc_anno_filter_harm.rdata') # mono_dc_filter
load('./final/seurat/b_cell/03-b_cell_anno_filter_harm.rdata')  # bcell_filter
load('./final/seurat/platelet/02-platelet_anno.rdata') # platelet_filter
load('./final/seurat/plasma/03-plasma_anno_filter_No_harm.rdata') # plasma_filter

# --------------------------- Filter and Harmony -------------------------------
pbmc_all <- merge(t_cell_merge,y = c(mono_dc_filter,bcell_filter,
                                     platelet_filter,plasma_filter),
                    merge.data = F,merge.dr = F)
# exclude the pSS_pah sample
Idents(pbmc_all) <- 'orig.ident'

# pbmc_all <- do_harmony(seu_obj = pbmc_all,harmony_slot = 'orig.ident',max.iter = 30,from_begin = T)
back_run(func = do_harmony, out_name = 'pbmc_all', job_name = 'pbmc_all',
         seu_obj = pbmc_all,harmony_slot = 'orig.ident',max.iter = 30,from_begin = T)

# -------------------------------- Analysis ------------------------------------
DimPlot(pbmc_all, group.by = 'subtype', label = T, raster = F,
        cols = get_color(length(table(pbmc_all$subtype))),repel = TRUE)
DimPlot(pbmc_all, group.by = 'main_type', label = T, raster = F,
        cols = get_color(length(table(pbmc_all$subtype))))
# -------------------------------- Modify meta ---------------------------------



#---------------------------------- Save ---------------------------------------
save(pbmc_all, file = './final/seurat/pbmc/03-pbmc_all_concat.rdata')


