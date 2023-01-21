library(Seurat)
library(patchwork)

# remove LL pair (test)
tmp1 <- tmp1[-2,]
tmp2 <- tmp2[-2,]

################################################################################
#
# plot the paired SLEDAI and cell type ratio
#
################################################################################

SLEDAI <- c(6,2,4,1,11,4,7,4,3,2,2,4,7,8,5,8,2,2,0,0,0,0)
names(SLEDAI) <- c('XH','XYY','XYY2','GW','WYF','WYF2','HXR','HXR2','ZPP','ZPP2',
                  'LGY','HXX','GZR','MXY','SQ','LL','LL2','WYY','QJY','ZMY1','ZH','ZS')
SLEDAI <- as.data.frame(SLEDAI)
meta <- unique(pbmc_all@meta.data[,c('treatment','orig.ident','pair','group')]) %>% remove_rownames() %>%column_to_rownames('orig.ident')
meta <-left_join(rownames_to_column(meta),rownames_to_column(SLEDAI),
                         by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
# unpaired
meta %>% filter(!treatment == 'HC') %>%
    ggboxplot( x = 'treatment' ,y = 'SLEDAI', color   = "treatment", 
                    line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(label = "p.signif", label.x = 1.5, label.y = 10) + ylab('SLEDAI index')
# paired 
paired_sledai <- data.frame(colnames= c('HXR','LL','WYF','XYY','ZPP') ,
                            before =c(7,8,11,6,3), after = c(4,2,4,2,2))
paired_sledai %>% ggpaired(cond1 = 'before', cond2 ='after', color = 'condition',
                           palette = 'npg', line.color = "gray", line.size = 0.5) +
    ylab('SELDAI index') + xlab('') +
    stat_compare_means(label = "p.signif", label.x = 1.5, label.y = 10, paired = T,method = 't.test')

# plasmablast 
tmp1<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='treated') %>%
    filter(!orig.ident == 'XYY2')
tmp2<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='untreated')
p1 <- data.frame(sample = tmp1$orig.ident, before=tmp2$Freq, after=tmp1$Freq ) %>%
    ggpaired( cond1 = 'before', cond2 = 'after', color = "condition", 
              line.color = "gray", line.size = 0.4, palette = "npg") + 
    ylab('Plasmablast / Plasma (%)') + NoLegend() + xlab('') +
    stat_compare_means(paired = F, method = 't.test',label.x = 1.5,label = "p.signif")

# proliferation T cells
tmp1<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
T_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                               before=tmp2$Ratio, after=tmp1$Ratio )
p2 <- T_pair_df %>% filter(subtype == 'T.prolife_T') %>%
    ggpaired(T_pair_df, cond1 = 'before', cond2 = 'after',
         color = "condition", line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(paired = T, method = 't.test',label.x = 1.5,label = "p.signif")+
    ylab('Prolifing T cell / T cell (%)') + NoLegend() + xlab('')

# macrophage
tmp1<- mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- mono_dc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
mono_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                           before=tmp2$Ratio, after=tmp1$Ratio )
p3 <- mono_pair_df %>% filter(subtype == 'Macrophage') %>%
ggpaired(mono_pair_df, cond1 = 'before', cond2 = 'after',
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(aes(label = paste0('p =', ..p.format..)),paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Macrophgae / Myeloid cell ') + NoLegend() + xlab('')

# CD4 Treg cells
tmp1<- cd4_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- cd4_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
cd4T_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                           before=tmp2$Ratio, after=tmp1$Ratio )
p4 <- cd4T_pair_df %>% filter(subtype == 'T.CD4.Treg') %>%
    ggpaired(cd4T_pair_df, cond1 = 'before', cond2 = 'after',legend = 'right',
             color = "condition", line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(aes(label = paste0('p =', ..p.format..)),paired = TRUE, method = 't.test',label.x = 1.5)+
    ylab('Treg cell / T cell (%)') + xlab('')

# memory B cells
tmp1<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
bcell_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                           before=tmp2$Ratio, after=tmp1$Ratio )
p5 <- bcell_pair_df %>% filter(subtype == 'B.mem') %>%
    ggpaired(bcell_pair_df, cond1 = 'before', cond2 = 'after',legend = 'right',
             color = "condition", line.color = "gray", line.size = 0.4, palette = "npg") +  
    stat_compare_means(aes(label = paste0('p =', ..p.format..)),paired = TRUE, method = 't.test',label.x = 1.5)+
    ylab('memory B cell / B cell (%)') + xlab('')

p1+p2+p4+plot_layout(ncol =3)

################################################################################
#
# colrelation plot
#
###############################################################################
t_prolife_ratio <-t_cell_merge@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='T.prolife_T') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname'))

plasmablast_ratio <-plasma_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(plasma_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='plasmablast') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname'))

macrophage_ratio <-mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='Macrophage') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname'))

treg_ratio <-cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='T.CD4.Treg') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname'))

sledai_colrelation_df <- data.frame(row.names = rownames(t_prolife_ratio), SLEDAI = t_prolife_ratio$SLEDAI,
                                    `T prolifing` = t_prolife_ratio$Ratio, Plasmablast = plasmablast_ratio$Ratio,
                                    treatment =t_prolife_ratio$treatment)
scatterplotMatrix(sledai_colrelation_df[,1:3],groups =sledai_colrelation_df$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='Correlation with SLEDAI', by.groups = F,
                  pch = c(16,16,16),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend =F , cex = 1, cex.labels = 1.5, cex.axis = 1.5,  text.width =0.01)

scatterplotMatrix(t_prolife_ratio[,c('Ratio','SLEDAI')],groups =t_prolife_ratio$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='Profiling T cell and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )

scatterplotMatrix(plasmablast_ratio[,c('Ratio','SLEDAI')],groups =plasmablast_ratio$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='Plasmablast and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )

scatterplotMatrix(macrophage_ratio[,c('Ratio','SLEDAI')],groups =macrophage_ratio$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='Plasmablast and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )

scatterplotMatrix(treg_ratio[,c('Ratio','SLEDAI')],groups =treg_ratio$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='Plasmablast and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )


################################################################################
#
# compare plasma before and after treatment 
#
################################################################################
Idents(plasma_filter) <- 'treatment'
plasma_filter_treat <- subset(plasma_filter, idents = 'treated')
plasma_filter_untreat <- subset(plasma_filter, idents = 'untreated')

marker_plasma_treat <- FindMarkers(plasma_filter, ident.1 = 'treated', ident.2 = 'untreated')
marker_plasma_treat %>% arrange(avg_log2FC) %>% head(20)
marker_plasma_treat %>% arrange(avg_log2FC) %>% tail(20)
marker_plasma_treat_filter <- marker_plasma_treat[!startsWith(rownames(marker_plasma_treat),'MT'),]
marker_plasma_treat_filter <- marker_plasma_treat_filter[!startsWith(rownames(marker_plasma_treat_filter),'RP'),]
EnhancedVolcano(marker_plasma_treat_filter,
                lab = rownames(marker_plasma_treat_filter),
                x = 'avg_log2FC', pointSize = 4, xlim = c(-3,3),
                y = 'p_val_adj', FCcutoff = 0.7,legendPosition = 'right', title = 'Plasmablast')


################################################################################
#
# marker gene of merge T cells
#
################################################################################
VlnPlot(t_cell_merge, stack = T,features = c('CD4','CD8A','CCR7','MKI67','FOXP3','ISG15'), group.by = 'subtype') + xlab('') + ylab('') + NoLegend()
DotPlot2(t_cell_merge, marker_list = c('CD4','CD8A','CCR7','MKI67','FOXP3','ISG15'), group.by = 'subtype') 




