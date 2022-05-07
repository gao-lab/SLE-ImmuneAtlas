
################################################################################
#
# Figure 4A and Extend 4A
#
################################################################################

# UMAP of T cell main type
DimPlot(t_cell_merge, group.by = 't_meta',raster=FALSE, cols = get_color(6,set = 'Paired',set_len = 6)) + NoAxes()
DimPlot(t_cell_merge, group.by = 'subtype',raster=FALSE, cols = get_color(18,set = 'Paired',set_len = 10),label = F) + NoAxes()



################################################################################
#
# Figure 4B
#
################################################################################
# Set T cell main meta
t_cell_merge$t_meta <- 'unknown'
t_cell_merge$t_meta[grepl(pattern = 'NK',x =t_cell_merge$subtype)] <- 'NK'
t_cell_merge$t_meta[grepl(pattern = 'CD4',x =t_cell_merge$subtype)] <- 'CD4 Tcell'
t_cell_merge$t_meta[grepl(pattern = 'CD8',x =t_cell_merge$subtype)] <- 'CD8 Tcell'
t_cell_merge$t_meta[grepl(pattern = 'prolife',x =t_cell_merge$subtype)] <- 'prolife T'
t_cell_merge$t_meta[grepl(pattern = 'MAIT',x =t_cell_merge$subtype)] <- 'MAIT'
t_cell_merge$t_meta[grepl(pattern = 'yd',x =t_cell_merge$subtype)] <- 'γδ Tcell'

# pub Fig 4b up
t_cell_merge@meta.data  %>% 
    group_by(orig.ident,t_meta) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4,5,6)]  %>%  distinct() ) %>%
    filter(group != 'treated') %>% filter(t_meta %in% c('CD8 Tcell','MAIT','prolife T')) %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group',
                      palette =c("#E3292C" ,"#4CAB45"))+ 
    facet_wrap(~t_meta,scales = "free",ncol = 3) + xlab('') + ylab('Ratio of T cells') +
    stat_compare_means(comparisons = list(c("SLE", "HC")),label = "p.signif" ) +  theme(legend.position="right")

# pub Fig 4b middle 
# please follow  analysis_celltypist.R to process the 'obs' object
dim(obs)
# Set T cell main meta
obs$t_meta <- 'unknown'
obs$t_meta[grepl(pattern = 'NK',x =obs$label)] <- 'NK'
obs$t_meta[grepl(pattern = 'CD4',x =obs$label)] <- 'CD4 Tcell'
obs$t_meta[grepl(pattern = 'CD8',x =obs$label)] <- 'CD8 Tcell'
obs$t_meta[grepl(pattern = 'prolife',x =obs$label)] <- 'prolife T'
obs$t_meta[grepl(pattern = 'MAIT',x =obs$label)] <- 'MAIT'
obs$t_meta[grepl(pattern = 'yd',x =obs$label)] <- 'γδ Tcell'
table(obs$t_meta)

obs %>% group_by(sample,t_meta) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(obs[,c(5,28)]  %>%  distinct(),by = c('sample' = 'sample') ) %>%
    filter(t_meta %in% c('CD8 Tcell','MAIT','prolife T')) %>%
    filter(!group %in% c('IFNbeta_stim','sle_treated') )  %>%
    ggpubr::ggboxplot(x='group',y='Ratio', fill = 'group', palette = 'npg')+ 
    facet_wrap(~t_meta,scales = "free",ncol = 3) + xlab('') + ylab('Ratio of T cells') +
    stat_compare_means(label = "p.signif" ,comparisons = list(c("sle", "sle_flare"),c('hc_child','sle_child'),c('sle','hc'),c('sle_flare','hc')))+
    theme(axis.text.x = element_text(angle = 30, hjust = 1),legend.position="right")

# stat_compare_means(comparisons = list(c("SLE", "HC")) ) 

# pub Fig 4b down 
# Box plot :(paired)cell t_meta ratio between treatment
tmp1<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,t_meta) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- t_cell_merge@meta.data  %>% 
    group_by(orig.ident,t_meta) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(t_cell_merge@meta.data[,c(1,4,5,6)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
data.frame(sample = tmp1$orig.ident, t_meta= tmp1$t_meta, 
           before=tmp2$Ratio, after=tmp1$Ratio ) %>% filter(t_meta %in% c('CD8 Tcell','MAIT','prolife T')) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',fill  = "condition", 
              line.color = "gray", line.size = 0.4, palette = 'npg') +
    stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Prolife Plasma ratio') + facet_wrap(~t_meta,scales= "free",ncol = 3 ) +
    theme(legend.position="right")


################################################################################
#
# Figure 4C: UMAP of CD8T cell subtype
#
################################################################################
DimPlot(cd8_filter, group.by = 'subtype',raster=FALSE, cols = 
            get_color(5,set = 'Paired',set_len = 5)) + NoAxes()


################################################################################
#
# Figure 4D: Marker of the IFN response 
#
################################################################################
DotPlot2(cd8_filter, marker_list =  head(rownames(marker_cd8_c10)), group.by = 'subtype')
VlnPlot(cd8_filter, features = head(rownames(marker_cd8_c10)),
        group.by = 'subtype', stack = T, cols = get_color(6))
