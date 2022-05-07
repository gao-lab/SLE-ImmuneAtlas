library(Seurat)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(org.Hs.eg.db)
library(topGO)

# different between SLE Treg and health Treg
Idents(cd4_filter) <- 'subtype'
treg <- subset(cd4_filter, idents = c('T.CD4.Treg'))
Idents(treg) <- 'treatment'
marker_Treg_sle <- FindMarkers(treg, ident.1 = 'untreated',ident.2 = 'HC', only.pos = T)
marker_Treg_sle %<>% filter(p_val_adj <0.05)  %>% arrange(-avg_log2FC) 
# rm(treg)

Idents(cd4_filter) <- 'subtype'
plan("multiprocess", workers = 10)
back_run(func =FindMarkers,out_name = 'marker_t_treg',job_name = 'marker_t_treg',
         object = cd4_filter, ident.1 = 'T.CD4.Treg', only.pos = T, logfc.threshold = 0.25)

plot_gsea(treg %>% subset(idents = 'treated', invert= T),group_by = 'treatment',focus = 'untreated',category = 'C2')

# clusterProfiler
gene =  bitr(rownames(marker_Treg_sle %>% head(150)), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(gene = gene$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 1)
dotplot(kegg,title="Enrichment KEGG_dot")

ego_ALL <- enrichGO(gene = gene$ENTREZID, OrgDb = org.Hs.eg.db, 
                    ont = "ALL",pAdjustMethod = "BH", readable = TRUE) 
dotplot(ego_ALL,showCategory=10,split='ONTOLOGY') +  facet_grid(ONTOLOGY~.,scale="free")

# key markers
DotPlot(treg, features = c('IL2RA','FOXP3','CTLA4'))
DotPlot2(treg, marker_list  = c('IL2RA','FOXP3','CTLA4'))
VlnPlot(treg, features = c('IL2RA','FOXP3','CTLA4'),stack = T)
DoHeatmap2(treg, marker_list = c('IL2RA','FOXP3','CTLA4'))


################################################################################
#
# Extend Data5: Treg score 
#
################################################################################
treg_score_genes <- c('CTLA4','IL2RA','LAG3','FOXP3')
DotPlot(treg, features = treg_score_genes, group.by = 'treatment')

treg_features <- list(treg_score_genes)
treg <- AddModuleScore(object = treg,features = treg_features,name = 'treg_score',ctrl = 10)
head(treg$treg_score1)

VlnPlot(treg, features = 'treg_score1')

treg$treatment <- as.vector(treg$treatment)
VlnPlot(treg, features = 'treg_score1', cols = get_color(3,set_len = 3)) + ggtitle('Treg Score')
