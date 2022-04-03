# plasma

setwd('/data/sle')
output_path <- './output_file/seurat/plasma/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/b_cell/b_cell_harmony_5iter_anno.rdata')

plasma <- subset(bcell_harmony, idents = c(5,11))
plasma <- do_seurat(plasma)
# save the seurat result without harmony directly
save(plasma, file = paste0(output_path,'plasma_re_clustering.rdata'))

DimPlot(plasma, group.by = 'orig.ident') + DimPlot(plasma, label = T)


#---------------------------- Vis Marker Genes ---------------------------------
DotPlot2(plasma, marker_list = b_marker_list)
VlnPlot(plasma, features = feats[1:3])


#------------------------------ Compare groups ---------------------------------
VlnPlot(plasma, group.by = 'group', stack = T,
        features = grep(rownames(plasma),pattern = '^IGH[AEGDM]',value = T))
VlnPlot(plasma, group.by = 'orig.ident', features = 'IGHM')


#------------------------- Filter Out the Doublets -----------------------------
plasma_fliter <- subset(plasma, idents = c(8,10,11), invert =T)
plasma_pic <- subset(plasma, idents = c(8,10,11))
save(plasma_pic, file = paste0(output_path, 'pic_plasma_recluster_no_harmony.rdata'))
Idents(plasma_fliter) <- 'group'
plasma_fliter <- subset(plasma_fliter, idents = 'pSS_pah', invert = T)
do_seurat_back(plasma_fliter,name = 'plasma_fliter')
plasma_fliter <- subset(plasma_fliter, idents = c(8,9), invert = T)

# add meta
plasma_fliter$disease <- 'SLE'
plasma_fliter$disease[which(plasma_fliter$group == 'HC')] <- 'HC'

save(plasma_fliter, file = paste0(output_path,'plasma_fliter_remove_ss_recluster_without_harmony.rdata'))

#------------------------- Vis Marker Genes After Filter -----------------------
DimPlot(plasma_fliter, group.by = 'orig.ident') + DimPlot(plasma_fliter, label = T)
plot_scdata(plasma_fliter, pal_setup = pal) + plot_scdata(plasma_fliter, color_by = "group", pal_setup = pal)
plot_scdata(plasma_fliter, color_by = "treatment", pal_setup = 'Set1')+ plot_scdata(plasma_fliter, color_by = "subtype", pal_setup = 'Set1')
DotPlot2(plasma_fliter, marker_list = b_marker_list)
VlnPlot(plasma_fliter, features = c('IHGM','IGHD','IGHG1','IGHG2','IGHG3','IGHG4',
                                       'IGHA1','IGHA2','IGHE'), stack = T, group.by = 'group')
marker_plasma_fliter_all <- FindAllMarkers(plasma_fliter,logfc.threshold = 0.25,min.pct = 0.2,min.diff.pct = 0.1)
seu_plot_heatmap(plasma_fliter, marker_plasma_fliter_all,sort_var =  c('seurat_clusters','group'),row_font_size = 8,
                 anno_var = c('seurat_clusters','group'),anno_colors = list('Set2','Set3'))
VlnPlot(plasma_fliter, features = c('IFI27'))

Idents(plasma_fliter) <- 'subtype'
marker_plasma_fliter_subtype <- FindAllMarkers(plasma_fliter,logfc.threshold = 0.25,
                                               min.pct = 0.2,min.diff.pct = 0.1)
seu_plot_heatmap(plasma_fliter, marker_plasma_fliter_subtype
                 ,sort_var =  c('subtype','disease','treatment'),row_font_size = 8,
                 anno_var = c('subtype','disease','treatment'),anno_colors = list('Set2','Set3','Set2'))
# ------------------------ Ratio of Plasma -------------------------------------
# overall plasma ratio
my_comparisons <- list(c("SLE", "HC"))
(table(plasma_fliter$orig.ident)/table(all_cell$orig.ident)[-c(10,11)] * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_fliter@meta.data[c('orig.ident', 'group', 'treatment', 'pair','disease')]) %>% 
    unique() %>% filter(!treatment == 'treated') %>%
    ggboxplot(x ='disease',y='Freq', color = 'disease',
              palette  = c("#00AFBB","#E7B800","#FC4E07"), add="jitter")+
    stat_compare_means(comparisons = my_comparisons,  method = "t.test" )

# plasmablast ratio 
plasmablast_meta <- plasma_fliter@meta.data %>% filter(subtype == 'Plasma.prolife')
(table(plasmablast_meta$orig.ident)/table(plasma_fliter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_fliter@meta.data[c('orig.ident', 'group', 'treatment', 'pair','disease')]) %>% 
    unique() %>% filter(!treatment == 'treated') %>% mutate(oppo = 100-Freq) %>%
    ggboxplot(x ='disease',y='Freq',fill = 'disease',
              palette  = c("#00AFBB","#E7B800","#FC4E07"), add="jitter")+ ylab('Prolife Plasma ratio')+
    stat_compare_means(comparisons = my_comparisons,  method = "t.test" )

# paired test
tmp1<- (table(plasmablast_meta$orig.ident)/table(plasma_fliter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_fliter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='treated') %>%
    filter(!orig.ident == 'XYY2')
tmp2<- (table(plasmablast_meta$orig.ident)/table(plasma_fliter$orig.ident) * 100) %>% 
    data.frame() %>% rename(orig.ident = Var1 ) %>% 
    left_join(plasma_fliter@meta.data[c('orig.ident', 'group', 'treatment', 'pair')]) %>% 
    unique() %>% filter(!pair == 'unpaired') %>% filter(treatment =='untreated')
plasma_pair_df <- data.frame(sample = tmp1$orig.ident, before=tmp2$Freq, after=tmp1$Freq )
ggpaired(plasma_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Prolife Plasma ratio')
    
back_run(load,out_name = 'all',job_name = 'loading all cell','./output_file/seurat/all_pbmc/all_cell_filter_anno.rdata')

# plasma subset ratio
plasma_fliter@meta.data %>% filter(group != 'pSS_pah') %>% 
ggplot( aes(x = group , fill = seurat_clusters))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(plasma_fliter$seurat_clusters))) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

plasma_fliter@meta.data %>% filter(group != 'pSS_pah') %>% 
    ggplot( aes(x = seurat_clusters , fill = group))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(plasma_fliter@meta.data$group))) + 
    scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

plasma_fliter@meta.data %>% filter(group != 'pSS_pah') %>% 
    ggplot( aes(x = seurat_clusters , fill = orig.ident))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(plasma_fliter$orig.ident))) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()
plasma_fliter@meta.data %>% filter(group != 'pSS_pah') %>% 
    ggplot( aes(x = orig.ident , fill = seurat_clusters))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(plasma_fliter$orig.ident))) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

plot_stat(plasma_fliter, plot_type = "prop_multi", group_by = "group", pal_setup = c("sienna","bisque3"))

# ------------------------------ Health vs SLE----------------------------------
Idents(plasma_fliter) <- 'disease'
marker_plasma_fliter_disease <- FindAllMarkers(plasma_fliter,logfc.threshold = 0.25,
                                               min.pct = 0.2,min.diff.pct = 0.1,only.pos = T)
marker_plasma_fliter_disease %<>% filter(p_val_adj < 0.05)
seu_plot_heatmap(plasma_fliter, marker_plasma_fliter_disease,n = 15,
                 sort_var =  c('disease','treatment','subtype'),row_font_size = 8,
                 anno_var = c('disease','treatment','subtype','orig.ident'),anno_colors = list('Set2','Set3','Set2','Set3'))
VlnPlot(plasma_fliter, features = c('PLAAT4','PLAAT2'),split.by = 'pair',split.plot = T)
DotPlot(plasma_fliter, features = c('PLAAT4','PLAAT2'), group.by = 'disease')
plot_all_cluster_go(marker_plasma_fliter_disease, org = "human", ont = "BP")


# ------------------------------ Treatment vs none -----------------------------
Idents(plasma_fliter) <- 'treatment'
marker_plasma_fliter_treatment <- FindAllMarkers(plasma_fliter,logfc.threshold = 0.25,
                                               min.pct = 0.2,min.diff.pct = 0.1,only.pos = T)
seu_plot_heatmap(plasma_fliter, marker_plasma_fliter_treatment,n = 15,
                 sort_var =  c('disease','treatment','subtype','seurat_clusters'),row_font_size = 8,
                 anno_var = c('treatment','disease','subtype','seurat_clusters'),anno_colors = list('Set2','Set3','Set2','Set3'))
Idents(plasma_fliter) <- 'pair'
plasma_fliter_treatment <- subset(plasma_fliter, idents = 'unpaired', invert =T)
VlnPlot(plasma_fliter_treatment, features = c('PLAAT2'),group.by = 'pair',split.by = 'treatment')
VlnPlot(plasma_fliter_treatment, features = c('PLAAT2'),group.by = 'treatment')
DotPlot(plasma_fliter_treatment, features = c('PLAAT2','PLAAT4'),group.by = 'treatment')
DotPlot(plasma_fliter_treatment, features = c('PLAAT2','PLAAT4','IFI27'),group.by = 'orig.ident')

# pathway analysis
marker_plasma_fliter_treatment_filter <- marker_plasma_fliter_treatment %>% 
plot_all_cluster_go(marker_plasma_fliter_treatment, org = "human", ont = "BP")

DotPlot2(plasma_fliter, marker_list =c('IFI27','IFIT1','IFI44','IFI44L','PLAAT2'
                                       ,'PLAAT4','CD44','CD52',
                                       'S100A4','HLA-DRB5','NFKBIA') , group.by = 'treatment' )

# ------------------------ Plasma harmony(abondon)------------------------------
# I find use harmony may loss some info so I decide not use harmony in plasma
# run Harmony ------
plasma_fliter_harmony <- plasma_fliter %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

plasma_fliter_harmony <- plasma_fliter_harmony %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
save(plasma_fliter_harmony, file = paste0(output_path, 'plasma_fliter_harmony.rdata'))

DimPlot(plasma_fliter_harmony,label = T) + DimPlot(plasma_fliter_harmony, group.by = 'subtype', label = T) 
DimPlot(plasma_fliter_harmony, group.by = 'Phase')
VlnPlot(plasma_fliter_harmony, features = c('IHGM','IGHD','IGHG1','IGHG2','IGHG3','IGHG4',
                                    'IGHA1','IGHA2','IGHE'), stack = T)
marker_plasma_fliter_all <- FindAllMarkers(plasma_fliter_harmony,logfc.threshold = 0.25,min.pct = 0.2,min.diff.pct = 0.1)
seu_plot_heatmap(plasma_fliter_harmony, marker_plasma_fliter_all,sort_var =  c('seurat_clusters','group'),row_font_size = 8,
                 ,anno_var = c('seurat_clusters','group'),anno_colors = list('Set2','Set3'))




