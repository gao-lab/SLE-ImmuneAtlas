
setwd('/data/sle')
output_path <- './output_file/seurat/plasma/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/plasma/plasma_fliter_remove_ss_recluster_without_harmony.rdata')
all_bcr <- readRDS('./output_file/scRepertotre/BCR/all_BCR_combined.rds')

# We need to handle the barcode difference between scRepertoire and seurat
head(all_bcr$GW_SLE$barcode)
# SAMPLE + GROUP + BARCODE
Cells(plasma_fliter) %>% head()
scRepertoire_barcode <- plasma_fliter@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 
    
bcr_plasma <- plasma_fliter
bcr_plasma <- RenameCells(bcr_plasma, new.names = scRepertoire_barcode$scRepertoire)
bcr_plasma <- combineExpression(all_bcr, bcr_plasma, group.by='sample',cloneCall="gene",
                                cloneTypes = c(Single =1, Small=5, Large=10, Hyperexpanded=15),
                                proportion = F)

colorblind_vector <- colorRampPalette((c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
# colorblind_vector <- colorRampPalette(rev(c('red','green','blue','yellow')))

slot(bcr_plasma, "meta.data")$cloneType <- factor(slot(bcr_plasma, "meta.data")$cloneType,
                                              levels = c("Hyperexpanded (10 < X <= 20)",
                                                         "Large (5 < X <= 10)",
                                                         "Small (1 < X <= 5)",
                                                         "Single (0 < X <= 1)", NA))
DimPlot(bcr_plasma, group.by = "cloneType",split.by = 'treatment') +
# plot_scdata(bcr_plasma, color_by = "cloneType", split_by = "treatment") +
    scale_color_manual(values = c('#FF8C00','#FFDEAD','#FDF5E6'), na.value="#F5F5F5") + NoAxes()

clonalOverlay(bcr_plasma, reduction = "umap", 
              freq.cutpoint = 3, bins = 4, facet = "treatment") + 
    guides(color = FALSE)

# SLE and HC
ggboxplot(bcr_plasma@meta.data, x = 'disease',y = 'Frequency',fill = 'disease') + 
    ylim(c(0,20)) + xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("SLE", "HC")),  method = "wilcox.test",
                       label.y = 18,bracket.size = 0.5)
# remove Freq =1 (abondon)
bcr_plasma@meta.data %>% filter(Frequency > 1) %>%
ggboxplot(, x = 'disease',y = 'Frequency',fill = 'disease') + 
    ylim(c(0,20)) + xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("SLE", "HC")),  method = "wilcox.test",
                       label.y = 18,bracket.size = 0.5)


# before and after treatment 
bcr_plasma@meta.data %>% filter(!treatment == 'HC') %>% group_by(orig.ident) %>%
    # mutate(cell_num = n()) %>% mutate(clone_num = cell_num * Frequency ) 
ggboxplot( x = 'treatment',y = 'Frequency',fill = 'treatment') + 
    ylim(c(0,20)) + 
    xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("untreated", "treated")),method = 'wilcox.test'
                        ,label.y = 18,bracket.size = 0.5)

# remove Freq =1 (abondon)
bcr_plasma@meta.data %>% filter(!treatment == 'HC') %>% group_by(orig.ident) %>%
    filter(Frequency > 1) %>%
    # mutate(cell_num = n()) %>% mutate(clone_num = cell_num * Frequency ) 
    ggboxplot( x = 'treatment',y = 'Frequency',fill = 'treatment') + 
    ylim(c(0,20)) + 
    xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(comparisons =  list(c("untreated", "treated")), 
                       method = "t.test" ,label.y = 18,bracket.size = 0.5)

bcr_plasma@meta.data  %>% filter(!treatment == 'HC') %>% group_by(orig.ident) %>%
    filter(! orig.ident == 'XYY2') %>%filter(!pair == 'unpaired')  %>% 
    ggboxplot( x = 'pair',y = 'Frequency',color = 'treatment',palette = "jco") + 
    ylim(c(0,10)) + 
    xlab('') + ylab('Plasma clone type distribution') +
    stat_compare_means(aes(group = treatment),hide.ns = T)
    # stat_compare_means(comparisons =  list(c("untreated", "treated")), 
                       # method = "t.test" ,label.y = 8,bracket.size = 0.5) +
    
