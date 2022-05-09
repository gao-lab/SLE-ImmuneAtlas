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
# Figure 6C: Ratio of Mono.recruit
#
################################################################################
mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(!treatment == 'treated') %>% filter(subtype == 'Mono.CD14.recruit') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(~subtype,scales = "free",nrow = 1) + stat_compare_means(label = 'p.signif',label.x = 1.5)+
    theme_cowplot() + ylab('Ratio of Myeloid cells') + xlab('')





################################################################################
#
# Figure 6D: Marker gene of the subtype 
#
################################################################################



################################################################################
#
# Figure 6E: Heatmap of TF activity
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
                              t_tf_mono.mat[rownames(t_tf_mono.mat) %in%c('IRF1(+)','IRF3(+)','IRF7(+)','IRF9(+)','STAT1(+)','STAT2(+)',
                                                                              'XBP1(+)','JUN(+)','PAX5(+)'),] )
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


################################################################################
#
# Extend Data Fig 6 A: Ratio of macrophage
#
################################################################################
mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(!treatment == 'treated') %>% filter(subtype == 'Macrophage') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(~subtype,scales = "free",nrow = 1) + stat_compare_means(label = 'p.format',label.x = 1.5)+
    theme_cowplot() + ylab('Ratio of Myeloid cells') + xlab('')


################################################################################
#
# Extend Data Fig 6 B : Ratio of macrophage after treatment 
#
################################################################################
# mono_dc_filter@meta.data  %>% 
#     group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
#     mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
#     left_join(mono_dc_harm@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
#     ggpubr::ggboxplot(x='group',y='Ratio', color = 'treatment',
#                       palette =c("#00AFBB", "#E7B800", "#FC4E07"))+ 
#     facet_wrap(~subtype,scales = "free") 

# Box plot :(paired)cell subtype ratio between treatment
tmp1<- mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>% filter(subtype == 'Macrophage') %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
              fill  = "condition", line.color = "gray", line.size = 0.4,
              palette = c('#DA9494','#9FB1D4')) +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free",ncol = 5 ) +
    theme_cowplot()
rm(tmp1,tmp2)


################################################################################
#
# Extend Data Fig 6 C : Ratio of macrophage after treatment 
#
################################################################################
mono_dc_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(mono_dc_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(!treatment == 'treated') %>% filter(subtype == 'Mono.CD14.APOBEC3A+') %>% 
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#DA9494", "#B4D493"))+ 
    facet_wrap(~subtype,scales = "free",nrow = 1) + stat_compare_means(label = 'p.format',label.x = 1.5)+
    theme_cowplot() + ylab('Ratio of Myeloid cells') + xlab('')


################################################################################
#
# Extend Data Fig 6D and Fig 6F: CEMiTool
#
################################################################################
library(CEMiTool)
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
# use SCENIC 
# int_df <- read_csv('./scripts/SCENIC/pDC/pdc_high_confidence_regulate.csv',col_names = T)
# int_df <- int_df[,c(2,3)]

sample_annot <- data.frame(matrix(data = NA, ncol = 2, nrow = 430) )
colnames(sample_annot) <- c('SampleName', 'Class')
sample_annot$SampleName <- Cells(pdc_filter)
sample_annot$Class <- pdc_filter$treatment
cem <- cemitool(pdc_filter@assays$RNA@data %>%as.data.frame(), sample_annot, gmt_in,interactions = int_df,
                filter=TRUE, plot=TRUE, verbose=TRUE)
# back_run(func = cemitool,out_name = 'cem',job_name = 'cemitool',
#          pdc_filter@assays$RNA@data %>%as.data.frame(), sample_annot, gmt_in, 
#          filter=TRUE, plot=TRUE, verbose=TRUE)

generate_report(cem, directory="./tmp/CEMiTool_pDC_Report", force = T)
tmp <- show_plot(cem, "ora")
show_plot(cem, "interaction")





