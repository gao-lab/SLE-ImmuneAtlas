
setwd('/data/sle')
output_path <- './output_file/seurat/b_cell//'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/plasma/plasma_fliter_remove_ss_recluster_without_harmony.rdata')
all_bcr <- readRDS('./output_file/scRepertotre/BCR/all_BCR_combined.rds')

# We need to handle the barcode difference between scRepertoire and seurat
head(all_bcr$GW_SLE$barcode)
# SAMPLE + GROUP + BARCODE
Cells(bcell_filter) %>% head()
scRepertoire_barcode <- bcell_filter@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 

bcr_bcell <- bcell_filter
bcr_bcell <- RenameCells(bcr_bcell, new.names = scRepertoire_barcode$scRepertoire)
bcr_bcell <- combineExpression(all_bcr, bcr_bcell, group.by='sample',cloneCall="gene",
                                cloneTypes = c(Single =1, Small=10,Medium =20, Large=30, Hyperexpanded=100),
                                proportion = F)

# ------------------------------- Analysis -------------------------------------
slot(bcr_bcell, "meta.data")$cloneType <- factor(slot(bcr_bcell, "meta.data")$cloneType,
                                                  levels = c("Hyperexpanded (30 < X <= 100)",
                                                             "Large (20 < X <= 30)",
                                                             "Medium (10 < X <= 20)",
                                                             "Small (1 < X <= 10)",
                                                             "Single (0 < X <= 1)", NA))
DimPlot(bcr_bcell, group.by = "cloneType",split.by = 'treatment',pt.size = 0.4) +
    # plot_scdata(bcr_bcell, color_by = "cloneType", split_by = "treatment") +
    scale_color_manual(values = c('#B22222','#FF8C00','#FFDAB9','#FDF5E6','#F5F5F5'), na.value="#F5F5F5") + NoAxes()

# clonetype freq
## HC and sle
ggboxplot(bcr_bcell@meta.data, x = 'disease',y = 'Frequency',fill = 'disease') + 
    ylim(c(0,60)) + xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("SLE", "HC")),  method = "wilcox.test",
                       label.y = 58,bracket.size = 0.5)

# before and after treatment 
bcr_bcell@meta.data %>% filter(!treatment == 'HC') %>% group_by(orig.ident) %>%
    # mutate(cell_num = n()) %>% mutate(clone_num = cell_num * Frequency ) 
    ggboxplot( x = 'treatment',y = 'Frequency',fill = 'treatment') + 
    ylim(c(0,60)) + 
    xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("untreated", "treated")),method = 'wilcox.test'
                       ,label.y = 58,bracket.size = 0.5)
## paried sample 
bcr_bcell@meta.data  %>% filter(!treatment == 'HC') %>% group_by(orig.ident) %>%
    filter(!orig.ident %in% c('XYY2','WH','WH2')) %>%filter(!pair == 'unpaired')  %>% 
    ggboxplot( x = 'pair',y = 'Frequency',color = 'treatment',palette = "jco") + 
    ylim(c(0,30)) + 
    xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(aes(group = treatment),hide.ns = T)
# stat_compare_means(comparisons =  list(c("untreated", "treated")), 
# method = "t.test" ,label.y = 8,bracket.size = 0.5) +

#---------------------------------- Save ---------------------------------------
save(bcr_bcell, file = './output_file/scRepertotre/BCR/bcr_combine_seurat.rdata')


