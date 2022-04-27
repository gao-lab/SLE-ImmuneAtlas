
DotPlot(pdc_filter, features = TLR_genes)

VlnPlot(pdc_filter, features = c('HLA-C'))


# clusterProfiler
gene =  bitr(rownames(marker_pdc_sle %>% head(50)), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(gene = gene$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 1)
dotplot(kegg,title="Enrichment KEGG_dot")

ego_ALL <- enrichGO(gene = gene$ENTREZID, OrgDb = org.Hs.eg.db, 
                    ont = "ALL",pAdjustMethod = "BH", readable = TRUE) 
dotplot(ego_ALL,showCategory=10,split='ONTOLOGY') +  facet_grid(ONTOLOGY~.,scale="free")
