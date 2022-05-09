################################################################################
#
# Figure 5A
#
################################################################################

# UMAP of T cell main type
DimPlot(cd4_filter, group.by = 'subtype',raster=FALSE, cols = get_color(8,set = 'Paired',set_len = 8)) + NoAxes()


################################################################################
#
# Extened Data Fig5: IFN-response CD4 ratio 
#
################################################################################
cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>% filter(subtype == 'T.CD4.Treg') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(~subtype,scales = "free",ncol = 8) + 
    stat_compare_means( label = 'p.signif',label.x = 1.5 ) + 
    xlab('') + ylab('Ratio of CD4 T cell') + theme_cowplot()


################################################################################
#
# Figure 5B
#
################################################################################
DotPlot2(cd4_filter, marker_list =  head(rownames(marker_cd4_c10)), group.by = 'subtype')



################################################################################
#
# Figure 5C : rario change of treg
#
################################################################################
# in core datastet
cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>% filter(subtype %in% c('T.CD4.Treg')) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#E3292C" ,"#4CAB45"))+ 
    facet_wrap(~subtype,scales = "free",ncol = 3) + xlab('') + ylab('Ratio of CD4+ T cells') +
    stat_compare_means(comparisons = list(c("SLE", "HC")),label = "p.signif" ) +  theme(legend.position="right")
# we need consider compare with whole T cell or CD4 T cell only

# in extend dataset
obs %>% group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    filter(label %in% c('T.CD4.Treg')) %>%
    filter(!group %in% c('IFNbeta_stim','sle_treated') )  %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group', palette = 'npg')+ 
    facet_wrap(~label,scales = "free",ncol = 3) + xlab('') + ylab('Ratio of CD4+ T cells') +
    stat_compare_means(label = "p.signif" ,comparisons = list(c("sle", "sle_flare"),c('hc_child','sle_child'),c('sle','hc'),c('sle_flare','hc')))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1),legend.position="right")


# compare with treatment 
tmp1<- cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>% filter(subtype %in% c('T.CD4.Treg')) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',fill  = "condition", 
              line.color = "gray", line.size = 0.4, palette = 'npg') +
    stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Ratio of CD4+ T cells') + facet_wrap(~subtype,scales= "free",ncol = 3 ) +
    theme(legend.position="right")


################################################################################
#
# Figure 5e: kegg of Treg
#
################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
gene <-  bitr(rownames(marker_Treg_sle %>% head(50)), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',pvalueCutoff = 0.1)
dotplot(kegg,title="Enrichment KEGG_dot")
ego_ALL <- enrichGO(gene = gene$ENTREZID, OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH", readable = TRUE)
dotplot(ego_ALL,showCategory=10,split='ONTOLOGY') +  facet_grid(ONTOLOGY~.,scale="free")


