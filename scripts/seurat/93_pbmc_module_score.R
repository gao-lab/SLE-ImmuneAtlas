
library(presto)
library(edgeR)

# ------------------------------- Define genes --------------------------------
allegs = get('GO:0034340', org.Hs.egGO2ALLEGS)
ifn_genes = unlist(mget(allegs,org.Hs.egSYMBOL)) %>% unique()
# IFN_genes
ifn_genes.2 <- intersect(genes, row.names(pbmc_all))
pbmc_all <- AddModuleScore(pbmc_all,features = list(ifn_genes.2),name = 'IFN_score')

# ------------------------------Main type --------------------------------------
IFN_df <- pbmc_all@meta.data %>% select(group,main_type,treatment,IFN_score1,subtype)
# res <- wilcox.test(IFN_score1 ~ group, data = IFN_df,
#                    exact = FALSE)
IFN_count <- IFN_df %>% select(IFN_score1)
IFN_count$IFN_score1 <- (IFN_count$IFN_score1 - mean(IFN_count$IFN_score1) )/sd(IFN_count$IFN_score1)
IFN_seu <- CreateSeuratObject( as.data.frame(t(as.matrix(IFN_count))), meta.data = IFN_df)
Idents(IFN_seu) <- 'treatment'
FindMarkers(IFN_seu, ident.1 = 'untreated', ident.2 = 'HC',logfc.threshold = 0,min.cells.feature = 1)
FindMarkers(IFN_seu, ident.1 = 'treated', ident.2 = 'HC',logfc.threshold = 0,min.cells.feature = 1)

IFN_logFC_df <- data.frame(matrix(NA,nrow = 0 ,ncol = 5 ))
colnames(IFN_logFC_df) <- c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')
for(cell_type in unique(IFN_seu$main_type)){
    Idents(IFN_seu) <- 'main_type'
    tmp <- subset(IFN_seu, idents = cell_type)
    Idents(tmp) <- 'treatment'
    tmp1 <- FindMarkers(tmp, ident.1 = 'untreated', ident.2 = 'HC',logfc.threshold = 0,
                        min.cells.feature = 0,min.pct = 0,min.diff.pct = 0)
    tmp2 <- FindMarkers(tmp, ident.1 = 'treated', ident.2 = 'HC',logfc.threshold = 0,
                        min.cells.feature = 0, min.pct = 0,min.diff.pct = 0)
    rownames(tmp1) <- paste0(cell_type,'_untreated')
    rownames(tmp2) <- paste0(cell_type,'_treated')
    IFN_logFC_df <- rbind(IFN_logFC_df,tmp1)
    IFN_logFC_df <- rbind(IFN_logFC_df,tmp2)
}

IFN_logFC_df.new <- IFN_logFC_df %>% rownames_to_column( var = "cell_type") %>%
    separate(cell_type,c('subtype','group'),sep = '_') %>% mutate(FoldChange = exp(avg_log2FC))
# pub
ggplot(IFN_logFC_df.new, aes(x=subtype,y=group)) + geom_tile(aes(fill= avg_log2FC)) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 1) +
    theme(axis.line=element_blank(),
          # axis.text.x=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())

# ------------------------------Sub type --------------------------------------
IFN_df.2 <- pbmc_all@meta.data %>% select(group,subtype,treatment,IFN_score1,subtype)
# res <- wilcox.test(IFN_score1 ~ group, data = IFN_df.2,
#                    exact = FALSE)
IFN_count.2 <- IFN_df.2 %>% select(IFN_score1)
IFN_count.2$IFN_score1 <- (IFN_count.2$IFN_score1 - mean(IFN_count.2$IFN_score1) )/sd(IFN_count.2$IFN_score1)
IFN_seu.2 <- CreateSeuratObject( as.data.frame(t(as.matrix(IFN_count.2))), meta.data = IFN_df.2)
Idents(IFN_seu.2) <- 'treatment'
# FindMarkers(IFN_seu.2, ident.1 = 'untreated', ident.2 = 'HC',logfc.threshold = 0,min.cells.feature = 1)
# FindMarkers(IFN_seu.2, ident.1 = 'treated', ident.2 = 'HC',logfc.threshold = 0,min.cells.feature = 1)

IFN_logFC_df.2 <- data.frame(matrix(NA,nrow = 0 ,ncol = 5 ))
colnames(IFN_logFC_df.2) <- c('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')
for(cell_type in unique(IFN_seu.2$subtype)){
    Idents(IFN_seu.2) <- 'subtype'
    tmp <- subset(IFN_seu.2, idents = cell_type)
    Idents(tmp) <- 'treatment'
    tmp1 <- FindMarkers(tmp, ident.1 = 'untreated', ident.2 = 'HC',logfc.threshold = 0,
                        min.cells.feature = 0,min.pct = 0,min.diff.pct = 0)
    tmp2 <- FindMarkers(tmp, ident.1 = 'treated', ident.2 = 'HC',logfc.threshold = 0,
                        min.cells.feature = 0, min.pct = 0,min.diff.pct = 0)
    rownames(tmp1) <- paste0(cell_type,'_untreated')
    rownames(tmp2) <- paste0(cell_type,'_treated')
    IFN_logFC_df.2 <- rbind(IFN_logFC_df.2,tmp1)
    IFN_logFC_df.2 <- rbind(IFN_logFC_df.2,tmp2)
}

IFN_logFC_df.2.new <- IFN_logFC_df.2 %>% rownames_to_column( var = "cell_type") %>%
    separate(cell_type,c('subtype','group'),sep = '_') %>% mutate(FoldChange = exp(avg_log2FC))

ggplot(IFN_logFC_df.2.new, aes(x=subtype,y=group)) + geom_tile(aes(fill= avg_log2FC)) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 1) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    theme(axis.line=element_blank(),
          # axis.text.x=element_blank(),
          # axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          # legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
