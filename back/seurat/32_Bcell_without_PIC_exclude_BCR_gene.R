

setwd('/data/sle')
output_path <- './output_file/seurat/b_cell/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/b_cell/b_cell_harmony_5iter_anno.rdata')

#------------------------------ Re Analysis ------------------------------------
bcell_filter <- subset(bcell_harmony, idents = c(0,1,2,3,4,8,10,12))
bcell_filter <- do_seurat(bcell_filter)

#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
bcell_filter <- bcell_filter %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

bcell_filter <- bcell_filter %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
save(bcell_filter, file = paste0(output_path, 'bcell_filter_harmony_exclude_BCR_gene.rdata'))

load(file = paste0(output_path, 'bcell_filter_harmony_exclude_BCR_gene.rdata'))
# DimPlot(bcell_filter,label = T) + DimPlot(bcell_filter, group.by = 'subtype', label = T) 
# DimPlot(bcell_filter, group.by = 'Phase')
bcell_filter <- subset(bcell_filter, idents = c(9,10),invert = T)

# add meta
bcell_filter$disease <- 'SLE'
bcell_filter$disease[which(bcell_filter$group == 'HC')] <- 'HC'

load('./output_file/seurat/b_cell/final_bcell.rdata')

#----------------------------------- Vis marker --------------------------------
plot_scdata(bcell_filter, color_by = "treatment", pal_setup = 'Set2') +
    plot_scdata(bcell_filter, color_by = "subtype", pal_setup = pal)

DimPlot(bcell_filter,label = T) + DimPlot(bcell_filter, group.by = 'subtype', label = T) 
DimPlot(bcell_filter, group.by = 'Phase')
DotPlot2(bcell_filter, marker_list = c(b_marker_list,'IL21R' ,'BACH2','BCL6'))
DotPlot2(bcell_filter, marker_list = c('AICDA','IRF4' ,'PRDM1','PAX5','BCL6','BACH2','IFNAR1','IFNAR2'),group.by = 'treatment')

DotPlot2(bcell_filter, group.by = 'subtype', marker_list = c('PAX5','IL4R','TCL1A','IFIT1','IFIT3','ISG15','TLR7','IFI27','CD5','CD9','MME','IGHA2','IGHG1','AICDA'))
DotPlot2(subset(bcell_filter, idents = 5),marker_list = c('IFIT1','IFIT3','ISG15','IFI27','TLR7'), group.by = 'treatment')
DotPlot2(bcell_filter,marker_list = c('IFIT1','IFIT3','ISG15','IFI27','TLR7'), group.by = 'treatment')
DotPlot2(bcell_filter,marker_list = c('IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHA3','IGHE','IFIT1','IFIT3','ISG15'), group.by = 'treatment')
VlnPlot(bcell_filter, features =  c('IGHG1','IGHG2','IGHG3','IGHG4','IGHA1','IGHA2','IGHA3','IGHE','IFIT1','IFIT3','ISG15'), stack = T,group.by = 'treatment')

#----------------------------------- Ratio -------------------------------------
# bar plot
bcell_filter@meta.data %>% filter(group != 'pSS_pah') %>% 
    ggplot( aes(x = seurat_clusters, fill = orig.ident))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(tmp_df$group))) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

# By Disease
bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    ggpubr::ggboxplot(x='disease',y='Ratio', fill = 'disease',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free") + stat_compare_means( method = 't.test',label.x = 1.2)

# By treatment
bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>%
    ggpubr::ggboxplot(x='treatment',y='Ratio', color = 'treatment',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free") + stat_compare_means( method = 't.test',label.x = 1.2)

# paired test
tmp1<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
plasma_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                             before=tmp2$Ratio, after=tmp1$Ratio )
ggpaired(plasma_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Prolife Plasma ratio') + facet_wrap(~subtype,scales= "free" )

#----------------------------------- Vis marker --------------------------------
bcell_filter$subtype <- 'unknown'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(0,1))] <- 'B.naive'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(2))] <- 'B.transition'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(3,8))] <- 'B.memory.IGHA'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(6,7))] <- 'B.memory.IGHG'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(4))] <- 'B.GC-like'
bcell_filter$subtype[which(bcell_filter$seurat_clusters %in% c(5))] <- 'B.IFN-response'
# bcell_filter$subtype[which(bcell_filter$old.ident %in% c())] <- 'B.IGHV4+'
# bcell_filter$subtype[which(bcell_filter$old.ident %in% c(6))] <- 'B.MZB'
DimPlot(bcell_filter,label = T) + DimPlot(bcell_filter, group.by = 'subtype', label = T) 

marker_b_filter_all <- FindAllMarkers(bcell_filter, min.pct = 0.2, logfc.threshold = 0.25, only.pos = T)

seu_plot_heatmap(bcell_filter, marker_b_filter_all
                 ,sort_var =  c('seurat_clusters','disease','treatment'),row_font_size = 8,
                 ,anno_var = c('seurat_clusters','disease','treatment'),anno_colors = list('Set2','Set3','Set2'))

#----------------------------------- Save File ---------------------------------
save(bcell_filter, file = './output_file/seurat/b_cell/final_bcell.rdata')
