################################################################################
#
# Figure 6A:UMAP of the Myeloid subtype  
#
################################################################################



################################################################################
#
# Figure 6B: Marker gene of the subtype 
#
################################################################################
mono_dc_filter$subtype %>% unique()
mono_dc_filter$subtype <- factor(mono_dc_filter$subtype, levels = c("Mono.CD14","Mono.CD14.LGALS2+","Mono.CD14.APOBEC3A+",
                                                                    "Mono.CD14.recruit","Mono.CD16","Macrophage",
                                                                    "cDC1","cDC2","pDC","HPSC"))

myeloid_marker <- c('CD14','LGALS2','APOBEC3A','CXCL10','FCGR3A','C1QC','CD1C','CLEC9A','JCHAIN','CD34')
DotPlot2(mono_dc_filter, marker_list = myeloid_marker, group.by = 'subtype')



################################################################################
#
# Figure 6C: Heatmap of TF activity
#
################################################################################
tf_mono <- fread('./scripts/SCENIC/pDC/treatment.csv') %>% as_data_frame() %>% column_to_rownames('V1')


t_tf_mono <- transpose(tf_mono)
# get row and colnames in order
colnames(t_tf_mono) <- rownames(tf_mono)
rownames(t_tf_mono) <- colnames(tf_mono)
t_tf_mono.mat <- as.matrix(t_tf_mono)
# do column scale and than do row scale
t_tf_mono.mat <- scale(t_tf_mono.mat, center = T, scale=T)
t_tf_mono.mat <- t(scale(t(t_tf_mono.mat), center = T, scale=T))
t_tf_mono.mat.filter <- rbind(t_tf_mono.mat %>% head(20), 
                              t_tf_mono.mat[rownames(t_tf_mono.mat) %in%c('IRF1(+)','IRF3(+)','IRF4(+)','IRF9(+)','STAT1(+)','STAT3(+)',
                                                                              'PAX5(+)','TCF4(+)','XBP1(+)'),] )
pheatmap::pheatmap( t_tf_mono.mat.filter, color=colorRampPalette(c("#430D54","#25878D","#F7E620"))(100), 
                   breaks=seq(-1.2, 1.2, length.out = 100),treeheight_row=10, treeheight_col=10, 
                   border_color=NA,fontsize_row = 8)

library(viridis)
t_tf_mono %>% rownames_to_column('TFs') %>% arrange(-(untreated)) %>% 
    filter(TFs %in% c('IRF1(+)','IRF3(+)','IRF4(+)','IRF9(+)','STAT1(+)','STAT3(+)',
                      'PAX5(+)','TCF4(+)','XBP1(+)')) %>% 
    reshape2::melt(id.vars = 'TFs') %>% 
    ggplot( aes(TFs, variable, fill= value)) + 
        geom_tile() +
        scale_fill_viridis(discrete=FALSE) +
        theme_minimal()
tmp <- t_tf_mono %>% rownames_to_column('TFs') %>% reshape2::melt(id.vars = 'TFs') %>% group_by(variable) %>%
    slice_max(value,n = 5) %>% pull(TFs)
t_tf_mono %>% rownames_to_column('TFs') %>% reshape2::melt(id.vars = 'TFs') %>% 
    filter(TFs %in% tmp) %>%
    ggplot( aes(TFs, variable, fill= value)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    theme_minimal()
