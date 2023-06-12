library(Seurat)
library(tidyverse)
library(ggpubr)
library(car)
library(biomaRt)
library(pheatmap)
library(presto)
library(msigdbr)
library(fgsea)
library(magrittr)
library(EnhancedVolcano)
library(scRepertoire)
library(viridis)
library(patchwork)

source("./scripts/function_R/utils.R")


##################################################################
#                           Fig1
##################################################################
load("./final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata")

# Fig. 1b
f1b <- DimPlot(pbmc_all, group.by = "main_type",label = T, cols = get_color(12, set = "Paired", set_len = 12),
               raster = F, label.size =5, raster.dpi = c(2048, 2048)) + NoAxes() +
  theme(text = element_text(size = 15, colour = "black")) + labs(fill="") + ggtitle("")
ggsave("./Figure/f1b.svg", f1b, dpi = 600, width = 8, height = 6)

tmp <- DimPlot(pbmc_all, group.by = "subtype",label = T, cols =  get_color(41 ,set = "Paired",set_len = 12),
               raster = T, repel = T) + NoAxes()
ggsave("./Figure/tmp.svg", tmp, dpi = 600, width = 13, height = 6)


# Fig. 1c
f1c <- ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$sample_name, 
                                      fill =pbmc_all$main_type))+
  geom_bar(stat = "count", position = "fill") + labs(y = "proportions", x = "") + 
  scale_fill_discrete(labels= names(table(pbmc_all$main_type))) + xlab("")+ labs(fill="") + 
  scale_fill_manual(values=get_color(12 ,set = "Paired",set_len = 12)) +
  # scale_color_npg() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
         text = element_text(size = 15,colour = "black",face="bold"),
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("./Figure/f1c.svg", f1c, dpi = 600, width = 8, height = 6)


# Fig. 1d
Idents(pbmc_all) <- "subtype"
marker_pbmc <- c("CD3D","CCR7","IL7R","GZMA","NKG7","MKI67","KLRC1","CD79A","MZB1","LILRA4",
                 "ITGAX","CD86","CD14","FCGR3A","CD34","PPBP")
# table(pbmc_all$main_type)
# pbmc_all$main_type <- factor(pbmc_all$main_type, 
#                              levels = rev(c("T.naive","T.cyto","T.prolife","NK",
#                                             "Bcell","Plasma","pDC","cDC","Mono.CD14","Mono.CD16","HPSC",
#                                             "Platelet")))
# 
# f1d <- VlnPlot(pbmc_all, features = marker_pbmc, stack = T, group.by = "main_type", 
#              cols = get_color(16, set = "Paired", set_len = 12), combine = T) + 
#       NoLegend() + NoGrid() + ylab("") + xlab("")  
#   #theme(axis.text.y=element_text(angle=300, hjust=0.5, face = "bold", size = 12), axis.text.x=element_text(angle=270, hjust=0.5))
# 

pbmc_all@assays$RNA@data -> A
A[marker_pbmc, ] -> A
as.data.table(A, keep.rownames=TRUE) -> A
setnames(A, "rn", "Feat") -> A
part.expr.dt <- melt(A, id.vars="Feat", variable.name="Cell", value.name="Expr")

part.ident.dt <- data.table(Cell=names(pbmc_all$main_type), Ident=pbmc_all$main_type)
to.plot.dt <- merge(x=part.expr.dt, y=part.ident.dt, by="Cell", all=TRUE)

to.plot.dt$Ident <- factor(to.plot.dt$Ident, 
                           levels = c("T.naive","T.cyto","T.prolife","NK",
                                      "Bcell","Plasma","pDC","cDC","Mono.CD14",
                                      "Mono.CD16","HPSC", "Platelet"))
to.plot.dt$Feat <- factor(to.plot.dt$Feat, 
                          levels = c("CD3D","CCR7","IL7R","GZMA","NKG7","MKI67","KLRC1","CD79A","MZB1","LILRA4",
                                    "ITGAX","CD86","CD14","FCGR3A","CD34","PPBP"))

f1d <- ggplot(to.plot.dt, aes(x=Feat, y=Expr, fill=Feat)) +
    geom_violin(scale="width") +
    scale_fill_manual(values = get_color(16, set = "Paired", set_len = 12)) +
    facet_grid(Ident ~ ., switch="y", scales="free_y") + 
    theme_cowplot(font_size = 12) +
    theme(
        axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(),
        axis.text.x=element_text(angle=90, hjust=0.5),
        legend.position = "none", 
        panel.spacing = unit(0, "lines"), 
        panel.background = element_rect(fill = NA, color = "black"), 
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0)) +
    labs(y="", x="")
ggsave("./Figure/f1d.svg", f1d, dpi = 600, width = 8, height = 6)

celltypist_anno <- fread('./other_sc_data/extend_data_anno_result.csv',header = T) 
obs <- fread('./other_sc_data/predictions_obs.csv',header = T)

# obs$label <- 'unknown'
obs$label <- celltypist_anno$predicted_labels
obs$group <- factor(obs$group, levels = c('sle_treated','sle_flare','sle','hc','sle_child','hc_child','IFNbeta_stim'))

f1e <- obs %>%
    group_by(sample,label) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    filter(!group %in% c('IFNbeta_stim','sle_treated') ) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette = c("#e51a1c", "#ff7f00", "#3e8e93", "#e1c62f", "#7e6e85"))+ 
    facet_wrap(~label,scales = "free_y", ncol = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    stat_compare_means(label = "p.signif" ,comparisons = list(c("sle", "sle_flare"),c('hc_child','sle_child'),c('sle','hc'),c('sle_flare','hc')),
                       method = 't.test')
ggsave("./Figure/f1e.svg", f1e, dpi = 600, width = 20, height =12)

##################################################################
#                             Fig2
##################################################################
# Fig. 2a
# please see vdj/immcatation/analysis/SHM_sample.ipynb

# Fig. 2c
# please see vdj/immcatation/analysis/TCR_gene_usage.ipynb


##################################################################
#                             Fig3
##################################################################
# Fig. 3a
SLEDAI <- c(6,2,4,1,11,4,7,4,3,2,2,4,7,8,5,8,2,2,0,0,0,0)
names(SLEDAI) <- c("XH","XYY","XYY2","GW","WYF","WYF2","HXR","HXR2","ZPP","ZPP2",
                   "LGY","HXX","GZR","MXY","SQ","LL","LL2","WYY","QJY","ZMY1","ZH","ZS")
SLEDAI <- as.data.frame(SLEDAI)
meta <- unique(pbmc_all@meta.data[,c("treatment","orig.ident","pair","group")]) %>% remove_rownames() %>% column_to_rownames("orig.ident")
meta <-left_join(rownames_to_column(meta),rownames_to_column(SLEDAI),
                 by = c("rowname" = "rowname")) %>% column_to_rownames("rowname")
# paired 
paired_sledai <- data.frame(colnames= c("HXR","LL","WYF","XYY","ZPP") ,
                            before =c(7,8,11,6,3), after = c(4,2,4,2,2))
f3a <- paired_sledai %>% ggpaired(cond1 = "before", cond2 ="after", color = "condition",
                           palette = "npg", line.color = "gray", line.size = 0.5,
                           font.label = list(size = 14, color = "black")) +
    ylab("SLEDAI") + xlab("") +
    scale_x_discrete(breaks=c("before","after"), labels=c("Untreated", "Treated")) +
    stat_compare_means(label = "p.signif", label.x = 1.5, label.y = 10, paired = T, method = "t.test") +
    scale_fill_discrete(name = "Treatment", labels = c("Untreated", "Treated"))
f3a <- ggpar(f3a, legend.title = "Treatment", legend = "right") +
    scale_fill_discrete(name = "Treatment", labels = c("Untreated", "Treated"))
ggsave("./Figure/f3a.svg", f3a, dpi = 600, width = 6, height = 5)

# Fig. 3c
load("./final/seurat/t_cell/05-tcell_merge_final.rdata")
load("./final/seurat/plasma/03-plasma_anno_filter_No_harm.rdata")

t_prolife_ratio <-t_cell_merge@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(t_cell_merge@meta.data[,c(1,4)]  %>%  distinct()) %>%
  filter(subtype =="T.prolife_T") %>% left_join(meta %>% rownames_to_column(), by = c("orig.ident" = "rowname"))

plasmablast_ratio <-plasma_filter@meta.data  %>% 
  group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
  mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
  left_join(plasma_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
  filter(subtype =="plasmablast") %>% left_join(meta %>% rownames_to_column(), by = c("orig.ident" = "rowname"))

sledai_colrelation_df <- data.frame(row.names = rownames(t_prolife_ratio), SLEDAI = t_prolife_ratio$SLEDAI,
                                    `T prolifing` = t_prolife_ratio$Ratio, Plasmablast = plasmablast_ratio$Ratio,
                                    treatment =t_prolife_ratio$treatment)
# manually export as pdf
scatterplotMatrix(sledai_colrelation_df[,1:3],groups =sledai_colrelation_df$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main="Correlation with SLEDAI", by.groups = F,
                  pch = c(16,16,16),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend =F , cex = 1, cex.labels = 1.5, cex.axis = 1.5,  text.width =0.01)
f3c_part1 <- ggscatter(sledai_colrelation_df, x = "Plasmablast", y = "SLEDAI", color = 'treatment',
                       palette =  c("#88a16f", "#9DB0D3", "#DA9494"),
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE) + xlab('Plasmablast / Plasma ratio (%)') +
    stat_cor(method = "pearson", label.x = 10, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01)

f3c_part2 <- ggscatter(sledai_colrelation_df, x = "T.prolifing", y = "SLEDAI", color = 'treatment',
                       palette =  c("#88a16f", "#9DB0D3", "#DA9494"),
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + xlab('T.prolifing / T cell ratio (%)') +
    stat_cor(method = "pearson", label.x = 2, label.y = 12.5, p.accuracy = 0.001, r.accuracy = 0.01)
f3c <- f3c_part1 | f3c_part2
ggsave("./Figure/f3c.svg", f3c, dpi = 600, width = 8, height = 4)


# Fig. 3d
Idents(plasma_filter) <- "treatment"
plasma_filter_treat <- subset(plasma_filter, idents = "treated")
plasma_filter_untreat <- subset(plasma_filter, idents = "untreated")

marker_plasma_treat <- FindMarkers(plasma_filter, ident.1 = "treated", ident.2 = "untreated")
marker_plasma_treat %>% arrange(avg_log2FC) %>% head(20)
marker_plasma_treat %>% arrange(avg_log2FC) %>% tail(20)
marker_plasma_treat_filter <- marker_plasma_treat[!startsWith(rownames(marker_plasma_treat),"MT"),]
marker_plasma_treat_filter <- marker_plasma_treat_filter[!startsWith(rownames(marker_plasma_treat_filter),"RP"),]
f3d <- EnhancedVolcano(marker_plasma_treat_filter,
                lab = rownames(marker_plasma_treat_filter),
                x = "avg_log2FC", pointSize = 4, xlim = c(-3,3),
                y = "p_val_adj", FCcutoff = 0.7, legendPosition = "right",
                title = "Plasmablast", subtitle ="differential expression",caption="") + NoLegend()
ggsave("./Figure/f3d.svg", f3d, dpi = 600, width = 6, height = 8)


##################################################################
#                             Fig4
##################################################################
# Fig. 4a
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0034340
gene.data <- getBM(attributes=c("hgnc_symbol", "ensembl_transcript_id", "go_id"),
                   filters = "go", values = "GO:0034340", mart = ensembl)

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
GO_IFN_genes <- gene.data$hgnc_symbol %>% unique() %>% intersect(row.names(pbmc_all))
GO_IFN_genes = c("IFI27","IFI35","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2",
                  "IFITM3","IFNAR1","IFNAR2","IFNB1","IKBKE","IP6K2","IRAK1","IRF1",
                 "IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8",
                  "IRF9","ISG15","ISG20","JAK1","MUL1","MX1","MX2",
                  "OAS1","OAS2","OAS3","OASL",
                  "STAT1","STAT2","TREX1","USP18")

# IFN_genes
pbmc_all <- AddModuleScore(pbmc_all, features = list(GO_IFN_genes),name = "IFN_score")

pbmc_all@meta.data %>% as_tibble() %>%
  dplyr::select("treatment","subtype","orig.ident","IFN_score1") %>% 
  group_by(treatment,subtype) %>% summarise(mean_IFN = mean(IFN_score1)) %>%
  spread(treatment,mean_IFN) %>% mutate(treated_HC = log2(treated/HC), 
                                        untreated_HC = log2(untreated/HC) ,
                                        untreated_treated = log2(untreated/treated)) %>%
  column_to_rownames("subtype") %>% dplyr::select(4,5,6) %>%
  pheatmap(scale ="none",cluster_rows =F, cluster_cols = F,color =  c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
           breaks = bk)

pbmc_all@meta.data %>% 
  dplyr::select("treatment","subtype","orig.ident","IFN_score1") %>% 
  group_by(treatment,subtype) %>% summarise(mean_IFN= mean(IFN_score1)) %>%
  spread(treatment,mean_IFN) %>% mutate(treated_HC = 10*(treated-HC), 
                                        untreated_HC = 10*(untreated-HC) ,
                                        untreated_treated = 10*(untreated-treated)) %>%
  column_to_rownames("subtype") %>% 
  # select(1,2,3) %>%
  dplyr::select(5,4) %>%
  pheatmap(scale ="none",cluster_rows =F, cluster_cols = F, color = c(colorRampPalette(colors = c("blue","white"))(25),
                                                                      colorRampPalette(colors = c("white","red"))(25)),
           breaks = c(seq(-5,5,by=0.2)))


# Fig. 4b
load("./final/seurat/t_cell/04-CD4_Tcell_filter_anno.rdata")
load("./final/seurat/t_cell/04-CD8_Tcell_filter_anno.rdata")
load("./final/seurat/b_cell/04-final-b_cell_filter.rdata")
Idents(cd4_filter) <- "subtype"
Idents(cd8_filter) <- "subtype"
Idents(bcell_filter) <- "subtype"

p11 <- DimPlot(bcell_filter, cells.highlight = Cells(subset(bcell_filter, idents ="B.IFN-response"))) + NoLegend() + NoAxes() + ggtitle("B cell")
p12 <- DimPlot(cd4_filter,  cells.highlight =  Cells(subset(cd4_filter, idents ="T.CD4.IFN-response")))+ NoLegend()+ NoAxes()+ ggtitle("CD4 T cell")
p13 <- DimPlot(cd8_filter,  cells.highlight =  Cells(subset(cd8_filter, idents ="T.CD8.IFN-response" )))+ NoLegend()+ NoAxes()+ ggtitle("CD8 T cell")

p21 <- FeaturePlot(bcell_filter, features = c("IFI6"), min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p22 <- FeaturePlot(cd4_filter, features = c("IFI6"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p23 <- FeaturePlot(cd8_filter, features = c("IFI6"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoAxes()

p31 <- FeaturePlot(bcell_filter, features = c("MX1"), min.cutoff = 2, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p32 <- FeaturePlot(cd4_filter, features = c("MX1"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p33 <- FeaturePlot(cd8_filter, features = c("MX1"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoAxes()

p41 <- FeaturePlot(bcell_filter, features = c("ISG15"), min.cutoff = 2, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p42 <- FeaturePlot(cd4_filter, features = c("ISG15"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoLegend()+ NoAxes()
p43 <- FeaturePlot(cd8_filter, features = c("ISG15"),min.cutoff = 1, cols = c("grey","red"),pt.size = 0.6)+ NoAxes()

f4b <- (p11 | p12 | p13 ) /(p21 |p22 |p23) /(p31 |p32 |p33) /(p41 |p42 |p43)
ggsave("./Figure/f4b.svg", f4b, dpi = 600, width = 14, height = 13)

# Fig. 4c
tmp <- table(pbmc_all$orig.ident, pbmc_all$subtype) %>% data.frame()
tmp <- reshape2::dcast(tmp,formula = Var1~ Var2, value.var = 'Freq' ) %>% 
    column_to_rownames('Var1') %>% mutate(all = rowSums(.)) 
# for(i in c(1:dim(tmp)[1])){
#     for(j in c(1:(dim(tmp)[2]-1))){
#         # print(tmp[i,j])
#         tmp[i,j] <- tmp[i,j]/tmp[i,dim(tmp)[2]]
#     }
# }
# tmp <- select(tmp, -all)
rowSums(tmp)
tmp2 <-left_join(rownames_to_column(tmp),rownames_to_column(SLEDAI),
                 by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
# cor(tmp2) %>% View()

b_sum <- table(bcell_filter$orig.ident) %>% as.data.frame()
cd4_sum <- table(cd4_filter$orig.ident) %>% as.data.frame()
cd8_sum <- table(cd8_filter$orig.ident) %>% as.data.frame()
b_sum$Var1 == cd4_sum$Var1
cd4_sum$Var1 == cd8_sum$Var1
tmp <- tmp[match(b_sum$Var1 %>% as.vector(), rownames(tmp)),]
rownames(tmp) == cd4_sum$Var1

tmp$`B.IFN_response` <- tmp$`B.IFN-response`/b_sum$Freq
tmp$`T.CD4.IFN_response` <- tmp$`T.CD4.IFN-response`/cd4_sum$Freq
tmp$`T.CD8.IFN_response` <- tmp$`T.CD8.IFN-response`/cd8_sum$Freq

# pheatmap(tmp, scale = 'row')
# pheatmap(tmp, scale = 'column')
# pheatmap(tmp, scale = 'none')

ifn_ratio_df <- tmp[,colnames(tmp) %>% str_detect('IFN')] * 100


# ifn_ratio_df <- read.csv(file = "./tmp/ifn_ratio_df.csv", row.names = 1)
ifn_ratio_df <-left_join(rownames_to_column(ifn_ratio_df), rownames_to_column(meta),
                         by = c("rowname" = "rowname")) %>% column_to_rownames("rowname")
# scatterplotMatrix(ifn_ratio_df[1:3], groups = ifn_ratio_df$treatment ,
#                   smooth = list(spread = T, lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
#                   spread = FALSE, main="IFN-response and SLEDAI", by.groups = F,
#                   pch = c(15,16,17), col = c("#88a16f", "#9DB0D3", "#DA9494"), diagonal=list(method ="boxplot"), 
#                   legend = F)
# fit<-lm(sledai~`B.IFN-response` + `T.CD8.IFN-response` + `T.CD4.IFN-response` ,data=ifn_ratio_df)
# summary(fit)
# plot(fit)

f4c_part1 <- ggscatter(ifn_ratio_df, x = "B.IFN_response", y = "SLEDAI", color = 'treatment',
                       palette =  c("#88a16f", "#9DB0D3", "#DA9494"),
                       add = "reg.line",  # Add regressin line
                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                       conf.int = TRUE) + xlab('B.IFN-response / B cell ratio') +
    stat_cor(method = "pearson", label.x = 1, label.y = 11, p.accuracy = 0.001, r.accuracy = 0.01)

f4c_part2 <- ggscatter(ifn_ratio_df, x = "T.CD4.IFN_response", y = "SLEDAI", color = 'treatment',
                       palette =  c("#88a16f", "#9DB0D3", "#DA9494"),
                       add = "reg.line",  # Add regressin line
                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                       conf.int = TRUE) + xlab('T.CD4.IFN-response / CD4 T ratio (%)') +
    stat_cor(method = "pearson", label.x = 2, label.y = 11, p.accuracy = 0.001, r.accuracy = 0.01)

f4c_part3 <- ggscatter(ifn_ratio_df, x = "T.CD8.IFN_response", y = "SLEDAI", color = 'treatment',
                       palette =  c("#88a16f", "#9DB0D3", "#DA9494"),
                       add = "reg.line",  # Add regressin line
                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                       conf.int = TRUE) + xlab('T.CD8.IFN-response / CD8 T ratio (%)') +
    stat_cor(method = "pearson", label.x = 0.3, label.y = 11, p.accuracy = 0.001, r.accuracy = 0.01)
f4c <- f4c_part1 | f4c_part2 | f4c_part3
ggsave("./Figure/f4c.svg", f4c, dpi = 600, width = 12, height = 4)

f4c_part4 <- ggboxplot(ifn_ratio_df, x = 'treatment', y = 'B.IFN_response', fill = 'treatment',
                        palette = c('#da9494','#9db0d3','#8aa370')) + xlab('') + ylab('B.IFN-response / B cell ratio (%)') + NoLegend() +
                stat_compare_means(comparisons = list(c('HC','untreated'),c('untreated','treated')), method = 'wilcox', label = "p.signif") 
f4c_part5 <- ggboxplot(ifn_ratio_df, x = 'treatment', y = 'T.CD4.IFN_response', fill = 'treatment',
          palette =  c('#da9494','#9db0d3','#8aa370')) + xlab('') + ylab('T.CD4.IFN-response / CD4 T ratio (%)') + NoLegend() +
    stat_compare_means(comparisons = list(c('HC','untreated'),c('untreated','treated')), method = 'wilcox', label = "p.signif")
f4c_part6 <- ggboxplot(ifn_ratio_df, x = 'treatment', y = 'T.CD8.IFN_response', fill = 'treatment',
          palette =  c('#da9494','#9db0d3','#8aa370')) + xlab('') + ylab('T.CD8.IFN-response / CD8 T ratio (%)') + NoLegend() +
    stat_compare_means(comparisons = list(c('HC','untreated'),c('untreated','treated')), method = 'wilcox', label = "p.signif")
f4c_down <- f4c_part4 | f4c_part5 | f4c_part6
ggsave("./Figure/f4c_down.svg", f4c_down, dpi = 600, width = 12, height = 4)


# Fig. 4d 
## Define functions 
makeGeneList <- function(filename){
  gl <- read.table(filename)
  # gl <- readr::read_csv(filename)
  y <- grepl("^RPS|^RPL|^MRPL|^MRPS|^MT-", gl$feature)
  gl <- gl[!y, ]
  gl <- gl  %>% filter(!group == "HC") %>% # not use HC
    dplyr::select(feature, logFC) %>% arrange(desc(logFC)) %>% deframe()
  return(gl)
}

plotGSEA_Hallmark <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL, keep_significant_only = FALSE) {
  require(ggplot2)
  if(!is.null(cols)){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols = gg_color_hue(dplyr::n_distinct(gsea$group, na.rm = TRUE))
  } else {
    cols = c("#696969","#F5DEB3","#800000")
  }   
  
  gsea$NES[which(is.na(gsea$NES))] <- 0
  gsea$pval[which(is.na(gsea$pval))] <- 1
  gsea$padj[which(is.na(gsea$padj))] <- 1
  gsea$ranking[which(is.na(gsea$ranking))] <- 0
  gsea <- gsea[order(gsea$ranking),]      
  # gsea_spl <- split(gsea, gsea$group)
  # if(!is.null(group_ref)){
  #     gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
  #     gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
  # } else {
  #     gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
  # }
  # names(gsea_spl) <- NULL
  # 
  # gsea <- do.call(rbind, gsea_spl)
  
  if (keep_significant_only){
    gseax <- split(gsea, gsea$pathway)
    for (i in 1:length(gseax)){
      if (all(gseax[[i]]$pval >= 0.05)|all(gseax[[i]]$padj >=0.25)){
        gseax[[i]] <- NA        
      }
    }
    gseax <- gseax[!is.na(gseax)]
    gsea <- do.call(rbind, gseax)
    cols = c("#F5DEB3","#800000")
  }
  # gsea <- gsea[order(gsea$ranking), ]
  gsea %<>% arrange(desc(group),-NES)
  gsea$pathway <- gsub("HALLMARK_|", "", gsea$pathway)
  gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.25)] <- "NotSig"
  gsea$group[which(gsea$pval < 0.05 & gsea$padj >= 0.25)] <- "NotSig"
  gsea$group[which(gsea$pval >= 0.05)] <- "NotSig"
  gsea$group <- factor(gsea$group, levels = c("NotSig", "untreated", "treated"))
  
  x_lim_min <- abs(ceiling(min(-log10(gsea$padj))))
  x_lim_max <- abs(ceiling(max(-log10(gsea$padj))))
  
  if(x_lim_min > x_lim_max){
    xval1 <- x_lim_min * -1
    xval2 <- x_lim_min
  } else {
    xval1 <- x_lim_max * -1
    xval2 <- x_lim_max
  }
  
  g <- ggplot(gsea, aes(x = -log10(padj)*sign(NES), y = reorder(pathway, ranking), fill  = group, size = abs(NES))) + 
    scale_fill_manual(values = cols) +
    geom_point(shape = 21, colour = "black",alpha = 0.7,position = position_jitter(width = 0,height = 0)) +
    labs(x = expression(paste("Signed", " -log" ["10"], " adjusted P value")), y = "Hallmarks") +
    theme_bw() +
    geom_vline(xintercept = 0) + geom_vline(xintercept = -log10(0.25),linetype="dashed",colour = "gray") +
    geom_vline(xintercept = -log10(0.25)*-1,linetype="dashed",colour = "gray") + xlim(xval1, xval2) +
    scale_size_area(max_size = 7, limits = c(0,3)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(),
          text = element_text(face = "bold", size = 16, color = "black"))+
    labs(fill="Group") 
  
  g$data <- g$data[order(g$data$group, na.last = TRUE), ]
  return(g)
}

comparisons = list()
res = list()
h <- as.list(kelvinny::parse_gmt("/data/sle/source/gsea/h.all.v7.2.symbols.gmt"))
for (c in celltypes){
  for (g in groups){
    comparisons[[c]][[g]] = makeGeneList(paste0("/data/sle/final/gsea/bcell/",celltype,"_",group,"gsea.csv"))
    res[[c]][[g]] <- fgsea(pathways = h, stats = comparisons[[c]][[g]], nperm = 10000, minSize = 0, maxSize =1000) %>%
      as.data.frame()
  }
}

result = res
# names(result) <- gsub(".csv", "",files)
for(i in 1:length(comparisons)){
  result[[i]] <- lapply(result[[i]], function(x){
    x$ranking <- -log10(x$pval)*sign(x$NES) 
    x <- x[order(x$pathway), ]
    return(x)
  })
}

result <- lapply(result ,function(x){
  x[["untreated"]]$group = "untreated"
  x[["treated"]]$group = "treated"
  # x[["Moderate"]]$group = "Moderate"
  # x[["Mild"]]$group = "Mild"
  # x[["Asymptomatic"]]$group = "Asymptomatic"
  return(x)
})
result2 <- lapply(result, function(x) {
  y <- do.call(rbind, x)
  return(y)
})

f4d <- plotGSEA_Hallmark(result2$`B.IFN-response`,keep_significant_only = T, group_ref = "untreated")
ggsave("./Figure/f4d.svg", f4d, dpi = 600, width =10, height = 4)

#Fig. 4e
f4e <- DotPlot2(bcell_filter, marker_list = c("STAT2","STAT1","IRF9","IFIT1","IFI44L",
                                       "IFI6","ISG15","EPSTI1","OAS2","IRF1","MX1",
                                       "IRF7","MX2","IFIT3","IFIT2","ISG20")) + 
  theme(axis.text=element_text(colour="black", size = 14))
ggsave("./Figure/f4e.svg", f4e, dpi = 600, width =10, height = 4)


##################################################################
#                             Fig5
##################################################################
# Fig. 5a and b
load("./final/scRepertoire/BCR/combined_bcr.rdata")
scRepertoire_barcode <- plasma_filter@meta.data %>% rownames_to_column("barcode_new")%>% 
  mutate(new_barcode = str_split_fixed(barcode_new ,"_",2)[,1])%>% 
  mutate(scRepertoire=paste0(orig.ident,"_",group,"_",new_barcode) )%>% dplyr::select(scRepertoire)
plasma_filter_renamed <- RenameCells(plasma_filter, new.names = scRepertoire_barcode$scRepertoire)
bcr_df <- rbindlist(combined_bcr)
plasma_filter_renamed <- combineExpression(combined_bcr, plasma_filter_renamed, addLabel =T,
                                   cloneCall="gene+nt", group.by = "sample", proportion = FALSE, 
                                   cloneTypes=c(Single=1, Small=10, Medium=10, Large=1000))
plasma_filter_renamed$cloneSize <- "Mutiple"
plasma_filter_renamed$cloneSize[plasma_filter_renamed$cloneType=="Single (0 < X <= 1)"] <- "Single"
f5a_part1 <-DimPlot(plasma_filter_renamed, group.by = "cloneSize",cols = rev(c("#92C5DE","#D6604D")), 
        split.by = "treatment",pt.size = 1.5) + NoAxes()
ggsave("./Figure/f5a_part1.svg", f5a_part1, dpi = 600, width =18, height = 6)


f5b_part1 <- ggplot(data = plasma_filter_renamed@meta.data, aes(x = treatment, fill = cloneSize))+
  geom_bar(stat = "count",position = "fill") + labs(y = "proportions", x = "") + 
  scale_fill_discrete(labels= names(table(plasma_filter_renamed$cloneSize))) + xlab("")+
  labs(fill="") +  scale_fill_manual(values= rev(c("#92C5DE","#D6604D"))) + theme_classic() +
  labs(fill="Clone size") +
  theme(text = element_text(size = 17), axis.text.x=element_text(colour="black"),
        legend.title = element_text(size = 15)) + ggtitle("Plasma")
ggsave("./Figure/f5b_part1.svg", f5b_part1, dpi = 600, width =6, height = 6)

bcell_filter_renamed <- combineExpression(combined_bcr, bcell_filter, 
                                  cloneCall="gene+nt", group.by = "sample", proportion = FALSE,
                                  cloneTypes=c(Single=1, Small=10, Medium=10, Large=1000))
bcell_filter_renamed$cloneSize <- "Mutiple"
bcell_filter_renamed$cloneSize[bcell_filter_renamed$cloneType=="Single (0 < X <= 1)"] <- "Single"
f5a_part2 <-DimPlot(bcell_filter_renamed, group.by = "cloneSize",cols = rev(c("#92C5DE","#D6604D")), 
                    split.by = "treatment",pt.size = 1.5) + NoAxes()
ggsave("./Figure/f5a_part2.svg", f5a_part2, dpi = 600, width =18, height = 6)

f5b_part2 <- ggplot(data = bcell_filter_renamed@meta.data, aes(x = treatment, fill = cloneSize))+
  geom_bar(stat = "count",position = "fill") + labs(y = "proportions", x = "") + 
  scale_fill_discrete(labels= names(table(plasma_filter_renamed$cloneSize))) + xlab("")+
  labs(fill="") + scale_fill_manual(values= rev(c("#92C5DE","#D6604D"))) + theme_classic() +
  labs(fill="Clone size") +
  theme(text = element_text(size = 17), axis.text.x=element_text(colour="black"),
        legend.title = element_text(size = 15)) + ggtitle("B cell")
ggsave("./Figure/f5b_part2.svg", f5b_part2, dpi = 600, width =6, height = 6)

# Fig. 5c
all_b <- merge(x= bcell_filter_renamed, y = plasma_filter_renamed)
all_b_bcr <- combineExpression(combined_bcr, all_b, 
                               cloneCall="gene+nt", group.by = "sample", proportion = FALSE, 
                               cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000))
bcr_df_filter <- bcr_df %>% filter(barcode %in% Cells(all_b_bcr))

Idents(all_b_bcr) <- "treatment"
hc_bcr <- subset(all_b_bcr,idents = "HC")
Idents(hc_bcr) <- "subtype"
f5c_part1 <- clonalOverlap2(hc_bcr, cloneCall="gene+nt", method="jaccard",title = "Health Control",limit = c(0,0.045)) +
  theme(axis.text=element_text(colour="black", size = 15), text = element_text(colour="black", size = 15))

untreat_bcr <- subset(all_b_bcr,idents = "untreated")
Idents(untreat_bcr) <- "subtype"
f5c_part2 <- clonalOverlap2(untreat_bcr, cloneCall="gene+nt", method="jaccard",title = "SLE Untreated",limit = c(0,0.045))

treat_bcr <- subset(all_b_bcr,idents = "treated")
Idents(treat_bcr) <- "subtype"
f5c_part3 <- clonalOverlap2(treat_bcr, cloneCall="gene+nt", method="jaccard",title = "SLE Treated",limit = c(0,0.045))
rm(hc_bcr, untreat_bcr, treat_bcr)

ggsave("./Figure/f5c_part1.svg", f5c_part1, dpi = 600, width =12, height = 6)
ggsave("./Figure/f5c_part2.svg", f5c_part2, dpi = 600, width =12, height = 6)
ggsave("./Figure/f5c_part3.svg", f5c_part3, dpi = 600, width =12, height = 6)

# Fig. 5d
SHM_info <- read.csv("./scripts/immcantaion/bcell_immacantation_SHM.csv")
contig <- str_split_fixed(SHM_info$sequence_id, pattern = "_",n = 2)[,1] 
SHM_info$barcode <-paste0(SHM_info$sample,"_",SHM_info$disease,"_", contig)

bcell_filter$barcode <- Cells(bcell_filter)
intersect(SHM_info$barcode , Cells(bcell_filter)) %>% length()

bcell_filter@meta.data %<>%  left_join(SHM_info, by = c("barcode" = "barcode"))
bcell_filter$mu_freq_seq_r[is.na(bcell_filter$mu_freq_seq_r)] <- 0

SHM_bcell <- bcell_filter@meta.data
SHM_bcell$treatment <- factor(SHM_bcell$treatment,levels = c("HC","untreated","treated"))
SHM_bcell <- SHM_bcell %>% drop_na(sample)

SHM_plasma <- plasma_filter_bcr@meta.data %>% left_join(SHM_info, by = c("barcode" = "barcode"))
SHM_plasma$treatment <- factor(plasma_filter_bcr$treatment, levels = c("HC","untreated","treated"))
SHM_plasma <- SHM_plasma %>% drop_na(sample)

shared_cols <- intersect(colnames(SHM_bcell),colnames(SHM_plasma))
SHM_all <- rbind(SHM_bcell[,shared_cols], SHM_plasma[,shared_cols])
SHM_all$subtype[grepl(SHM_all$subtype,pattern = "plasma")] <- "plasma"

compare_group <-  list(c("HC","untreated"),c("treated","untreated"),c("HC","treated"))
f5d <- SHM_all %>% filter(!subtype %in% c("B.naive","B.IFN-response","B.transition")) %>%
  filter(!c_call %in% c("IGHG4","IGHD")) %>% 
  ggboxplot( "treatment", "mu_freq_seq_r", fill = "treatment",
             palette = c("#B4D493", "#DA9494","#9FB1D4")) + ggtitle("Total mutations(replace)") +
  xlab("B cell subtype") + ylab("Mutation frequency") +
  theme_bw() + 
  facet_grid(subtype ~ c_call)  +
  stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                     comparisons =compare_group,label.y =c(0.1,0.12,0.14)) +
      xlab("") + theme_cowplot()+ theme(axis.text.x=element_text(angle=30, hjust=1), text = element_text(size = 17))
ggsave("./Figure/f5d.svg", f5d, dpi = 600, width =15, height = 12)

# Fig. 5e
# Please see vdj/immcatation/analysis/SHM_sample.ipynb


##################################################################
#                             Fig6
##################################################################
# Fig. 6a
tmp <- load(file = "final/scRepertoire/TCR/cd8_tcr.rdata")
combined_tcr.df.filter <-  filter(combined_tcr.df, barcode %in% Cells(cd8_tcr))
combined_tcr.filter <- split(combined_tcr.df.filter, f = combined_tcr.df.filter$sample)
combined_tcr.filter <- addVariable(combined_tcr.filter, name = "group", 
                                   variables = c("before", "before","before","after", "before", "before", "before", "after", "before", "HC", "before",
                                                 "before", "after","before", "before", "after", "after", "HC", "HC","before", "after", "HC"))
f6a_part1 <- DimPlot(cd8_tcr, group.by = "cloneType", split.by = "treatment",pt.size = 0.5)  + 
  scale_color_brewer(palette = "RdBu", direction = 1) + NoAxes()
ggsave("./Figure/f6a_part1.svg", f6a_part1, dpi = 600, width =22, height = 8)


meta <- occupiedscRepertoire(cd8_tcr, x.axis = "subtype", facet.by = "treatment",label = F,proportion = T, exportTable = T) 
f6a_part2 <- meta %>%ggplot(aes(x =subtype, y = value, fill = cloneType)) + geom_bar(stat = "identity", position = "fill") + 
  facet_grid(. ~ meta[, "treatment"]) + ylab("Proportion of Cells") + xlab("") +
  theme_bw() + labs(fill = "Clone Size") + 
  scale_fill_brewer(palette = "RdBu") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),panel.background = element_blank(),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 16))
ggsave("./Figure/f6a_part2.svg", f6a_part2, dpi = 600, width =14, height = 7)


# Fig. 6a
p1 <- scatterClonotype(combined_tcr.filter, cloneCall = "gene+nt",
                       x.axis = "ZPP", 
                       y.axis = "ZPP2",
                       dot.size = "total",
                       graph = "count",exportTable = F)
p2 <- scatterClonotype(combined_tcr.filter, cloneCall = "gene+nt",
                       x.axis = "WYF", 
                       y.axis = "WYF2",
                       dot.size = "total",
                       graph = "count",exportTable = F)
p3 <- scatterClonotype(combined_tcr.filter, cloneCall = "gene+nt",
                       x.axis = "HXR", 
                       y.axis = "HXR2",
                       dot.size = "total",
                       graph = "count",exportTable = F)

p4 <- scatterClonotype(combined_tcr.filter, cloneCall = "gene+nt",
                       x.axis = "XH", 
                       y.axis = "XYY",
                       dot.size = "total",
                       graph = "count",exportTable = F)

e5c <- scatterClonotype(combined_tcr.filter, cloneCall = "gene+nt",
                       x.axis = "LL", 
                       y.axis = "LL2",
                       dot.size = "total",
                       graph = "count",exportTable = F) +
  scale_color_manual(labels = c("Dual expand", "Untreat expand","Untreat single",
                                "Treat expand","Treat single"), values = c(get_color(5))) +
  xlab("Untreated") + ylab("Treated") +theme_cowplot() + ggtitle('Pair5') +
  theme(text = element_text(size = 18))

for (i in c(1:4)){
  name <- paste0("p",i)
  assign(name, get(name) +  
           scale_color_manual(labels = c("Dual expand", "Untreat expand","Untreat single",
                                         "Treat expand","Treat single"), values = c(get_color(5))) +
           xlab("Untreated") + ylab("Treated") +theme_cowplot() + ggtitle(paste0("Pair ",i)+
            theme(text = element_text(size = 18)))
  )
}

f6b <- (p1+p2) / (p3+p4)
ggsave("./Figure/f6b.svg", f6b, dpi = 600, width = 14, height = 9)
ggsave("./Figure/e5c.svg", e5c, dpi = 600, width = 7, height = 5)


##################################################################
#                             Fig7
##################################################################
# Fig. 7b
mip1_beta <- read.csv('/data/sle/data/blood_cytokine/MCP-1.csv')
g1 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL3', fill = 'group', palette = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test', label = "p.signif") + xlab('') +
  ylab('CCL3 in blood') + NoLegend()

g2 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL4', fill = 'group', palette = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test', label = "p.signif") + xlab('') +
  ylab('CCL4 in blood') + NoLegend()
g3 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL5', fill = 'group', palette = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test', label = "p.signif") + xlab('') +
  ylab('CCL5 in blood') + NoLegend()
g1 <- ggpar(g1, font.y = 16, font.xtickslab = 16)
g2 <- ggpar(g2, font.y = 16, font.xtickslab = 16)
g3 <- ggpar(g3, legend = 'right', legend.title = 'Group', font.y = 16, font.xtickslab = 16, font.legend = 16)
f7b <- g1 + g2 + g3
ggsave("./Figure/f7b.svg", f7b, dpi = 600, width = 14, height = 5)

# Fig. 7d
Idents(pbmc_all) <- 'subtype'
f7d_part1 <- DotPlot(pbmc_all, features = c('ESR1','AR','OR5AR1','PGR'), group.by = 'subtype')+ 
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  # scale_size(breaks = c(0, 25, 50, 75)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "grey", size=1),
                    axis.text = element_text(colour = 'black',size=12)) + xlab('') + 
  theme(axis.text.x = element_text(angle = 90))

f7d_part2 <- DotPlot(pbmc_all %>% subset(idents = 'pDC'), features = c('ESR1','AR'), group.by = 'treatment')+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "grey", size=1),
                    axis.text = element_text(colour = 'black',size=14)) + xlab('') + 
  theme(axis.text.x = element_text(angle = 90))

f7d_part3 <- DotPlot(pbmc_all %>% subset(idents = c('plasmablast','plasma','plasma.IgG','plasma.IgA')), 
        features = c('ESR1','AR'), group.by = 'treatment')+
  scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
  theme_bw() +theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "grey", size=1),
                    axis.text = element_text(colour = 'black',size=14)) + xlab('') + 
  theme(axis.text.x = element_text(angle = 90))
ggsave("./Figure/f7d_part1.svg", f7d_part1, dpi = 600, width = 5, height = 11)
ggsave("./Figure/f7d_part2.svg", f7d_part2, dpi = 600, width = 4, height = 4)
ggsave("./Figure/f7d_part3.svg", f7d_part3, dpi = 600, width = 4, height = 4)


##################################################################
#                             Supp1.5
##################################################################
# Supp. 2a
s2a <- DimPlot(bcell_filter, group.by = 'subtype', label = T, 
        cols = get_color(len = 7, set = 'Paired', set_len = 10), pt.size = 0.1) + 
  theme(plot.title = element_text(size = 5)) + 
  NoAxes() + ggtitle('')
ggsave("./Figure/s2a.svg", s2a, dpi = 600, width = 7, height = 5)


bcell_filter$subtype <- factor(bcell_filter$subtype, 
                               levels = rev(c("B.transition","B.naive","B.IFN-response",
                                              "B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-")))
# Supp. 2b
s1_5b <- DotPlot2(bcell_filter,
                marker_list = c('MME','IGHD','CXCR4','CCR7','TCL1A','IGHM','IFITM1',
                                'STAT1','CD40','CD24','CD22','CD69','VPREB3','CD27',
                                'IGHG1','CXCR3','IGHA1','IGHA2','MS4A1','FGR','FCRL5','TBX21','ITGAX'),
                group.by = 'subtype') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(colour="black", size = 12),
        panel.background = element_rect(colour = "grey", size=1),
        axis.text.x = element_text(angle = 90, color = 'black')) + xlab('')
ggsave("./Figure/s1_5b.svg", s1_5b, dpi = 600, width = 9, height =4)


# Supp. 2c
# load('final/seurat/t_cell/04-t_cell_merge.rdata')
# Idents(t_cell_merge) <- 'subtype'
# t_cell_merge <- subset(t_cell_merge, idents ="doublet", invert =TRUE)
# t_cell_merge <- subset(t_cell_merge, idents ="MAIT", invert =TRUE)
# t_cell_merge@meta.data %<>% mutate(across('subtype', str_replace, 'CD8.Tex', 'T.CD8.Tex'))
# s15c <- DimPlot(t_cell_merge, group.by = 'subtype', label = T, 
#         cols = get_color(len = 23, set = 'Paired', set_len = 10), pt.size = 0.1,raster=FALSE) + 
#     theme(plot.title = element_text(size = 5)) + 
#     NoAxes() + ggtitle('')
# ggsave("./Figure/s1_5c.png", s15c, dpi = 500, width = 15, height = 8)
# 
# s2d <- DotPlot2(t_cell_merge,
#                 marker_list = c('CD4','CD8A','KLRD1', 'NCAM1','NKG7',
#                                 'SELL', 'LEF1', 'ISG15','CCF7', 'RORC','FOXP3',
#                                 'CCL5', 'EOMES', 'TRAV1-2','TOP2A','TRVD1','TRVD2'
#                                 
#                                 
#                                 
#                                 # 'CD4',"CXCR3","TBX21","IFNG","RUNX2","GATA3" ,"IL13","IL5","IL17RB","LTB4R",
#                                 # "CCR8" ,"RORC","CCR6","RUNX1","PALLD","CXCR5","BCL6" , # cd4 
#                                 # "TOX","TOX2","BTLA","ICOS","CD40LG" ,"ICA1","FOXP3","CCL5",
#                                 # 'CD8','GZMK','GZMB','KLRD1','CD69','NKG7',
#                                 ),
#                 group.by = 'subtype') + theme_bw() +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.text=element_text(colour="black", size = 12),
#           panel.background = element_rect(colour = "grey", size=1),
#           axis.text.x = element_text(angle = 90)) + xlab('')

load('final/seurat/t_cell/04-CD4_Tcell_filter_anno.rdata')
s15c <- DimPlot(cd4_filter, group.by = 'subtype', label = T, 
                cols = get_color(len = 8, set = 'Paired', set_len = 10), pt.size = 0.1,raster=FALSE) + 
    theme(plot.title = element_text(size = 5)) + 
    NoAxes() + ggtitle('')
ggsave("./Figure/s1_5c.png", s15c, dpi = 500, width = 10, height = 6)

s15d <- DotPlot2(cd4_filter,
                marker_list = c('IFIT2','ISG15', 'TCF7', 'LEF1', 'PDCD1','TBX21',
                                'CXCR3', 'RORC', 'GATA3', 'GLNY', 'FOXP3'),
                group.by = 'subtype') + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(colour="black", size = 12),
          panel.background = element_rect(colour = "black", size=1),
          axis.text.x = element_text(angle = 90, color = 'black')) + xlab('')
ggsave("./Figure/s1_5d.svg", s15d, dpi = 600, width = 7, height = 3.7)

load('final/seurat/t_cell/04-CD8_Tcell_filter_anno.rdata')
s15e <- DimPlot(cd8_filter, group.by = 'subtype', label = T, 
                cols = get_color(len = 8, set = 'Paired', set_len = 10), pt.size = 0.1,raster=FALSE) + 
    theme(plot.title = element_text(size = 5)) + 
    NoAxes() + ggtitle('')
ggsave("./Figure/s1_5e.png", s15e, dpi = 500, width = 10, height = 6)

s15f <- DotPlot2(cd8_filter,
                 marker_list = c('IFIT2','ISG15','EOMES', 'TCF7', 'LEF1', 'KLRD1', 'CTLA4','HAVCR2'),
                 group.by = 'subtype') + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(colour="black", size = 12),
          panel.background = element_rect(colour = "black", size=1),
          axis.text.x = element_text(angle = 90, color = 'black')) + xlab('')
ggsave("./Figure/s1_5f.svg", s15f, dpi = 600, width = 6, height = 3.6)


load('final/seurat/mono_dc/03-mono_dc_anno_filter_harm.rdata')
s15g <- DimPlot(mono_dc_filter, group.by = 'subtype', label = T, 
                cols = get_color(len = 10, set = 'Paired', set_len = 10), pt.size = 0.1,raster=FALSE) + 
    theme(plot.title = element_text(size = 5)) + 
    NoAxes() + ggtitle('')
ggsave("./Figure/s1_5g.png", s15g, dpi = 500, width = 10, height = 6)

s15h <- DotPlot2(mono_dc_filter,
                 marker_list = c( "CD1C",  "FCER1A","CLEC9A","CST3", 'CD34', 'AVP' ,"C1QC","CSF1R","CD14",'APOBEC3A','LGALS2','HERC5',"FCGR3A",
                                  "ITGAX","CD86","IRF7","LILRA4","CLEC4C"),
                 group.by = 'subtype') + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_text(colour="black", size = 12),
          panel.background = element_rect(colour = "black", size=1),
          axis.text.x = element_text(angle = 90, color = 'black')) + xlab('')
ggsave("./Figure/s1_5h.svg", s15h, dpi = 600, width = 9, height = 4.5)
marker_mono.dc.recruit <- FindMarkers(mono_dc_filter, ident.1 = c('Mono.CD14.recruit'), only.pos = T,
                                          logfc.threshold = 0.25, min.diff.pct = 0.1)

##################################################################
#                             Supp4
##################################################################
# Supp. 4c
load(file = "final/scRepertoire/TCR/cd4_tcr.rdata") # cd4_tcr
load(file = "final/scRepertoire/TCR/combined_tcr.rdata") # combined_tcr
combined_tcr.df <- rbindlist(combined_tcr)
combined_tcr.df.filter <- filter(combined_tcr.df, barcode %in% Cells(cd4_tcr))
combined_tcr.filter <- split(combined_tcr.df.filter, f = combined_tcr.df.filter$sample)
combined_tcr.filter <- addVariable(combined_tcr.filter, name = "group", 
                                   variables = c("before", "before","before","after", "before", "before", "before", "after", "before", "HC", "before",
                                                 "before", "after","before", "before", "after", "after", "HC", "HC","before", "after", "HC"))
meta <- occupiedscRepertoire(cd4_tcr, x.axis = "subtype", facet.by = "treatment",label = F,proportion = T, exportTable = T) 
s4c <- meta %>%ggplot(aes(x = subtype, y = value, fill = cloneType)) + geom_bar(stat = "identity", position = "fill") + 
    facet_grid(. ~ meta[, "treatment"]) + ylab("Proportion of Cells") + xlab("") +
    theme_bw() + labs(fill = "Clone Size") + 
    scale_fill_brewer(palette = "RdBu") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(),panel.background = element_blank(),
          axis.text = element_text(colour = "black"),
          text = element_text(size = 16))
ggsave("./Figure/s4c.svg", s4c, dpi = 600, width =14, height = 7)

##################################################################
#                             Supp6
##################################################################
# Supp. 6a
s6a <- DotPlot(pbmc_all, features = c('CCR1'))

# Supp. 6b
Idents(mono_dc_filter) <- 'subtype'
tmp <- subset(mono_dc_filter, idents = c("Macrophage","Mono.CD14","Mono.CD14.APOBEC3A+","Mono.CD14.LGALS2+","Mono.CD14.recruit","Mono.CD16"))
s6b <- VlnPlot(object = tmp, features = "CCR1",  split.by = "group", group.by = 'subtype', pt.size = 0,
               cols=c("#88a16f", "#DA9494")) + xlab('') 
ggsave("./Figure/s6b.svg", s6b, dpi = 600, width = 9, height = 4.5)
