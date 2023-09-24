# plasma

setwd('/data/sle')
source('./scripts/function_R/utils.R')

plasma <- subset(bcell_harm, idents = c(6, 8))
plasma_blast <- subset(all_cell, idents = 20)

plasma <- merge(plasma, plasma_blast)
DotPlot2(plasma, marker_list = marker_list)
VlnPlot(plasma, features = feats[1:3])
DimPlot(plasma, group.by = 'orig.ident') + DimPlot(plasma, label = T)

#-------------------------- Filter and re-cluster ------------------------------
# plasma <- do_seurat(plasma, res = c(0.6))
plasma_harm <- do_harmony(plasma, harmony_slot = 'orig.ident', 
                          max.iter = 30, res =c(0.8), from_begin =T)

plasma_harm <- subset(plasma_harm, idents = c(0, 8, 10, 13, 14), invert = T)
plasma_harm <- do_harmony(plasma_harm, harmony_slot = 'orig.ident', 
                          max.iter = 30, res =c(0.8), from_begin =T)
VlnPlot(plasma_harm, features = marker_list, stack = T)
plasma_filter <- subset(plasma_harm, idents = c(9, 10, 11, 13), invert = T)
save(plasma_filter, file = './final/addtional/seurat/plasma/02-plasma_harm_filter.rdata')

DimPlot(plasma_filter, group.by = 'orig.ident') + DimPlot(plasma_filter, label = T)

#----------------------------- Compare key genes -------------------------------
VlnPlot(plasma_filter, stack = T,
        features = grep(rownames(plasma),pattern = '^IGH[AEGDM]',value = T))
DotPlot2(plasma_filter ,group.by = 'treatment',
         marker_list = grep(rownames(plasma),pattern = '^IGH[AEGDM]',value = T))
DotPlot2(plasma_filter, marker_list = c(plasma_marker, c('AICDA')))
VlnPlot(plasma_filter, features = feats)
# IFN related
VlnPlot(plasma_filter, features =c('IFIT1','IFIT3','ISG15','IFI27','TLR7') , stack = T)

# Heatmap: marker of clusters
back_run(FindAllMarkers, out_name = 'marker_all_plasma', job_name = 'marker_all_plasma',
         plasma_filter,logfc.threshold = 0.25, min.pct = 0.2,min.diff.pct = 0.1,only.pos = T)
marker_all_plasma %<>% filter(p_val_adj < 0.05)
marker_all_plasma %>% group_by(cluster) %>% slice_max(avg_log2FC ,n = 10) %>% View()
seu_plot_heatmap(plasma_filter, marker_all_plasma,n = 5,
                 sort_var =  c('seurat_clusters','group','treatment'),row_font_size = 8,
                 anno_var = c('seurat_clusters','group','treatment','orig.ident'),
                 anno_colors = list('Set2','Set3','Set2','Set1'))

#------------------------------- Anno subtype ----------------------------------
plasma_filter$subtype <- 'plasma'
plasma_filter$subtype[which(plasma_filter$seurat_clusters %in% c(3, 4))] <- 'plasmablast'
plasma_filter$subtype[which(plasma_filter$seurat_clusters %in% c(7))] <- 'plasma.IgA'
plasma_filter$subtype[which(plasma_filter$seurat_clusters %in% c(1, 2, 6, 8))] <- 'plasma'
plasma_filter$subtype[which(plasma_filter$seurat_clusters %in% c(0, 5, 12))] <- 'plasma.IgG'
table(plasma_filter$subtype)
DimPlot(plasma_filter, group.by = 'subtype',cols = get_color(len = 4,set = 'Set1',set_len = 4)) & NoAxes()

DotPlot2(plasma_filter, marker_list = c('IFI27','IGHA1','IGHA2','IGHG1','MKI67'),group.by = 'subtype')
VlnPlot(plasma_filter, features =c('IFIT1','IFIT3','ISG15','IFI27','IFI44L') , stack = T,group.by = 'subtype')

plasma_filter@meta.data  %>% 
    group_by(orig.ident, subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num / sample_num*100) %>%
    left_join(plasma_filter@meta.data[, c(1,4,5,6)] %>% distinct()) %>%
    filter(!treatment == 'treated') %>%
    ggpubr::ggboxplot(x = 'group', y = 'Ratio', fill = 'group', palette = c("#DA9494", "#B4D493")) + 
    facet_wrap(~subtype, scales = "free",nrow = 2) + stat_compare_means(label.x = 1.2)+
    theme_cowplot()

#------------------------- Vis Marker Genes After Filter -----------------------
DimPlot(plasma_filter, group.by = 'orig.ident') + DimPlot(plasma_filter, label = T)
plot_scdata(plasma_filter, pal_setup = 'Set2') + plot_scdata(plasma_filter, color_by = "group", pal_setup = 'Set2')
plot_scdata(plasma_filter, color_by = "treatment", pal_setup = 'Set1')+ plot_scdata(plasma_filter, color_by = "subtype", pal_setup = 'Set1')

VlnPlot(plasma_filter, features = c('IHGM','IGHD','IGHG1','IGHG2','IGHG3','IGHG4',
                                       'IGHA1','IGHA2','IGHE'), stack = T, group.by = 'group')

# Heatmap
Idents(plasma_filter) <- 'subtype'
marker_plasma_filter_subtype <- FindAllMarkers(plasma_filter,logfc.threshold = 0.25,
                                               min.pct = 0.2,min.diff.pct = 0.1)
seu_plot_heatmap(plasma_filter, marker_plasma_filter_subtype
                 ,sort_var =  c('subtype','disease','treatment'),row_font_size = 8,
                 anno_var = c('subtype','disease','treatment'),anno_colors = list('Set2','Set3','Set2'))

# ------------------------ Ratio of Plasma -------------------------------------
# Boxplot: plasma/pbmc ratio
my_comparisons <- list(c("SLE", "HC"))
as.data.frame((table(plasma_filter$orig.ident)/table(all_cell$orig.ident))* 100) %>% 
    rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% 
    # filter(!treatment == 'treated') %>%
    ggboxplot(x ='group',y ='Freq', fill = 'group',
              palette  = c("#00AFBB","#E7B800","#FC4E07"), add="jitter")+
    stat_compare_means(comparisons = my_comparisons ,method = 't.test') + ylab('Plasma / PBMC (%)')

# Boxplot: plasmablast/pbmc ratio (only before treatment )
plasmablast_meta <- plasma_filter@meta.data %>% filter(subtype == 'plasmablast')
as.data.frame((table(plasmablast_meta$orig.ident)/table(pbmc_all$orig.ident))* 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasmablast_meta[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique()  %>% mutate(oppo = 100-Freq) %>%
    # filter(!treatment == 'treated') %>%
    ggboxplot(x ='group',y='Freq',fill = 'group',
              palette  = c("#00AFBB","#E7B800","#FC4E07"), add="jitter")+ ylab('Plasmablast / PBMC (%)')+
    stat_compare_means(label.x = 1.5)

# Boxplot: plasmablast/plasma ratio (only before treatment )
as.data.frame((table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident))* 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasmablast_meta[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique()  %>% mutate(oppo = 100-Freq) %>%
    # filter(!treatment == 'treated') %>%
    ggboxplot(x ='group',y='Freq',fill = 'group',
              palette  = c("#00AFBB","#E7B800","#FC4E07"), add="jitter")+ ylab('Plasmablast / PBMC (%)')+
    stat_compare_means(comparisons = my_comparisons,  method = "t.test" )

# Boxplot: (paired)plasmablast/plasma ratio 
tmp1<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='treated') %>%
    filter(!orig.ident == 'XYY2')
tmp2<- (table(plasmablast_meta$orig.ident)/table(plasma_filter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='untreated')
data.frame(sample = tmp1$orig.ident, before=tmp2$Freq, after=tmp1$Freq ) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
             fill  = "condition", line.color = "gray", line.size = 0.4,
             palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
        ylab('Plasmablast / Plasma (%)')
    
# Boxplot: (paired)plasmablast/pbmc ratio 
tmp1<- (table(plasmablast_meta$orig.ident)/table(pbmc_all$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='treated') %>%
    filter(!orig.ident == 'XYY2')
tmp2<- (table(plasmablast_meta$orig.ident)/table(pbmc_all$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_filter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='untreated')
data.frame(sample = tmp1$orig.ident, before=tmp2$Freq, after=tmp1$Freq ) %>%
    ggpaired( cond1 = 'before', cond2 = 'after',
             fill  = "condition", line.color = "gray", line.size = 0.4,
             palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
        ylab('Plasmablast / PBMC (%)')

# Idents(plasma_filter) <- 'subtype'
# plot_stat(plasma_filter, plot_type = "prop_multi", group_by = "group", 
#           pal_setup = c("sienna","bisque3"),)

# ------------------------------ Health vs SLE----------------------------------

marker_plasma_hc.sle <- FindMarkers(plasma_filter,ident.1 = 'HC',ident.2 = 'SLE',group.by = 'group')
marker_plasma_hc.sle %<>% filter(p_val_adj <0.05) %>% arrange(avg_log2FC)

VlnPlot(plasma_filter, features = c('PLAAT4','PLAAT2'),split.by = 'pair',split.plot = T)
VlnPlot(all_cell,features = c('PLAAT4','PLAAT2'), group.by = 'main_type', split.by = 'group')
VlnPlot(plasma_filter, features = c('PLAAT4'),split.by = 'pair',split.plot = T)
VlnPlot(all_cell,features = c('PLAAT4'), group.by = 'main_type', split.by = 'group')
DotPlot(plasma_filter, features = c('PLAAT4','PLAAT2'), group.by = 'group')
plot_all_cluster_go(marker_plasma_filter_disease, org = "human", ont = "BP")


# ------------------------------ Treatment vs none -----------------------------
Idents(plasma_filter) <- 'treatment'
marker_plasma_filter_treatment <- FindAllMarkers(plasma_filter,logfc.threshold = 0.25,
                                               min.pct = 0.2,min.diff.pct = 0.1,only.pos = T)
seu_plot_heatmap(plasma_filter, marker_plasma_filter_treatment, n = 15,
                 sort_var =  c('disease','treatment','subtype','seurat_clusters'),row_font_size = 8,
                 anno_var = c('treatment','disease','subtype','seurat_clusters'),
                 anno_colors = list('Set2','Set3','Set2','Set3'))
Idents(plasma_filter) <- 'pair'
plasma_filter_treatment <- subset(plasma_filter, idents = 'unpaired', invert =T)
VlnPlot(plasma_filter_treatment, features = c('PLAAT2'),group.by = 'pair',split.by = 'treatment')
VlnPlot(plasma_filter_treatment, features = c('PLAAT2'),group.by = 'treatment')
DotPlot(plasma_filter_treatment, features = c('PLAAT2','PLAAT4'),group.by = 'treatment')
DotPlot(plasma_filter_treatment, features = c('PLAAT2','PLAAT4','IFI27'),group.by = 'orig.ident')

# pathway analysis
marker_plasma_filter_treatment_filter <- marker_plasma_filter_treatment %>% 
plot_all_cluster_go(marker_plasma_filter_treatment, org = "human", ont = "BP")

DotPlot2(plasma_filter, marker_list =c('IFI27','IFIT1','IFI44','IFI44L','PLAAT2'
                                       ,'PLAAT4','CD44','CD52',
                                       'S100A4','HLA-DRB5','NFKBIA') , group.by = 'treatment')

# ------------------------ Plasma harmony(abondon)------------------------------
# I find use harmony may loss some info so I decide not use harmony in plasma
# run Harmony ------
plasma_filter_harmony <- plasma_filter %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

plasma_filter_harmony <- plasma_filter_harmony %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
# save(plasma_filter_harmony, file = paste0(output_path, 'plasma_filter_harmony.rdata'))

DimPlot(plasma_filter_harmony,label = T) + DimPlot(plasma_filter_harmony, group.by = 'orig.ident', label = T) 
DimPlot(plasma_filter_harmony, group.by = 'Phase')
VlnPlot(plasma_filter_harmony, features = c('IHGM','IGHD','IGHG1','IGHG2','IGHG3','IGHG4',
                                    'IGHA1','IGHA2','IGHE'), stack = T)
marker_plasma_filter_all <- FindAllMarkers(plasma_filter_harmony,logfc.threshold = 0.25,min.pct = 0.2,min.diff.pct = 0.1)
seu_plot_heatmap(plasma_filter_harmony, marker_plasma_filter_all,sort_var =  c('seurat_clusters','group'),row_font_size = 8,
                 ,anno_var = c('seurat_clusters','group'), anno_colors = list('Set2','Set3'))



#------------------------ Save the Anno Objcet ---------------------------------
save(plasma_filter, file = './final/addtional/seurat/plasma/03-plasma_anno_filter_harm.rdata')

