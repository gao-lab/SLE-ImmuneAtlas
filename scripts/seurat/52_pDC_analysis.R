library(Scillus)
library(CEMiTool)

setwd('/data/sle')
output_path <- './output_file/seurat/mono_dc/'
source('./scripts/function_R/utils.R')

Idents(mono_dc_filter) <- 'subtype'
pdc <- subset(mono_dc_filter, idents = 'pDC')
pdc <- do_seurat(pdc, scale_all = T)
DimPlot(pdc) + DimPlot(pdc, group.by = 'orig.ident')


#---------------------------- Harmony Batch Remove -----------------------------
# run Harmony ------
pdc_harmony <- pdc %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20) %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
save(pdc_harmony, file = 'final/seurat/pDC/01-pDC_raw_harm.rdata')
pdc_filter <- subset(pdc_harmony, idents = c(4,6,7,8),invert =T)

DimPlot(pdc_filter,group.by = 'orig.ident') + 
    DimPlot(pdc_filter, label = T) 
plot_scdata(pdc_filter, color_by = "treatment", pal_setup = 'Set2') + 
    plot_scdata(pdc_filter, color_by = "seurat_clusters", pal_setup = 'Set1')

# FeaturePlot(pdc, features = c("CST3","IRF7","ITGAX","CD86","LILRA4","CLEC4C","CD1C","FCER1A"))


#----------------------------- Check Marker Gene -------------------------------
DotPlot2(pdc_filter, marker_list = marker_list)
VlnPlot(pdc_filter, features = feats)
DotPlot(pdc_filter, features =  grep(rownames(pdc), pattern = '^TLR', value = T))
VlnPlot(pdc_filter, features =  grep(rownames(pdc), pattern = '^TLR', value = T)[-10], stack = T)

FeaturePlot(pdc_filter, features = c('AXL','FCER1A','CST3','CD1C'))
FeaturePlot(pdc_filter, features = c('CD3D','CD79A','CD14','PPBP'))

DotPlot(pdc_filter, features = c('IRF9','STAT1','STAT2'), 
        split.by = 'treatment', cols = c("red", "blue", "green"))

DotPlot2(pdc_filter,marker_list  =c('IRF9','STAT1','STAT2'), group.by = 'treatment')


#----------------------------- Find Marker Gene --------------------------------
# marker_pdc_c6 <- FindMarkers(pdc_filter, ident.1 = 6, only.pos = T)
marker_pdc_sle <-  FindMarkers(pdc_filter, ident.1 = 'untreated', only.pos = F,group.by = 'treatment') %>%
    filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))

marker_pdc_all <- FindAllMarkers(pdc_filter, only.pos = T)
marker_pdc_all %<>% filter(p_val_adj < 0.05)  %>% top_n(avg_log2FC  , n = 10) %>% View()
pdc_marker <- marker_pdc_sle %>% filter(p_val_adj < 0.05)  %>% top_n(avg_log2FC  , n = 10) %>% rownames()
VlnPlot(pdc_filter, features = c('IRF9','STAT1','STAT2'))

DotPlot2(pdc_filter, marker_list = pdc_marker, group.by = 'group') + theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))

seu_plot_heatmap(pdc_filter, marker_pdc_sle %>% filter(p_val_adj < 0.05)  %>% top_n(avg_log2FC  , n = 10) %>% rownames(),
                 sort_var =  c('group','treatment'),row_font_size = 8,
                 anno_var = c('group','treatment'),anno_colors = list('Set2','Set3'))

Idents(pdc_filter) <- 'group'
gsea_reslut <-plot_gsea(pdc_filter, group_by = 'group', focus = 'SLE',
                        title = 'Platelet GSEA enrichment',category = 'H')

gsea_reslut

#----------------------------- Pathway enrichment ------------------------------
# gsva -----
gsva_result  <- analyse_sc_clusters(pdc_filter, verbose = TRUE)
pathway_expression <- pathways(gsva_result )
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

plot_gsva_pathway(gsva_result , pathway_id = rownames(max_difference)[11])
plot_gsva_heatmap(gsva_result , max_pathways =20,
                  # pathway_ids = relevant_pathways, # limit to these pathways
                  margins = c(6,20), # adapt the figure margins in heatmap.2
                  dendrogram = "col", # only plot column dendrogram
                  scale = "row", # scale for each pathway
                  key = FALSE, # don't display the color key
                  lwid=c(0.1,4)) # remove the white space on the left

plot_gsva_pca(gsva_result)
#------------------------------ Anno pDC Subtype -------------------------------
pdc$subtype <- 'unknown'

#----------------------------- Check the Ratio ---------------------------------


# bar plot
plot_stat(pdc_filter, plot_type = "prop_multi", pal_setup = "Set2",group_by = 'treatment')
pdc_filter@meta.data %>% filter(group != 'pSS_pah') %>% 
    ggplot( aes(x = seurat_clusters, fill = orig.ident))+
    geom_bar(stat = 'count',position = 'fill')+
    labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pdc_filter$orig.ident))) + 
    # scale_fill_manual(values =c('#90EE90','#DC143C','#F4A460')) +
    labs(fill="group") + coord_flip()

# By Disease
pdc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pdc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    ggpubr::ggboxplot(x='disease',y='Ratio', fill = 'disease',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~seurat_clusters,scales = "free",ncol = 4) + stat_compare_means( method = 't.test',label.x = 1.2)

# By treatment
pdc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pdc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>%
    ggpubr::ggboxplot(x='treatment',y='Ratio', color = 'treatment',
                      palette =c("#FC4E07", "#00AFBB"))+ 
    facet_wrap(~subtype,scales = "free") + stat_compare_means( method = 't.test',label.x = 1.2)

# paired test
tmp1<- pdc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pdc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='treated')%>% 
    filter(!orig.ident == 'XYY2')
tmp2<- pdc_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,seurat_clusters) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(pdc_filter@meta.data[,c(1,12,13,14,18)]  %>%  distinct()) %>%
    filter(!pair == 'unpaired') %>% filter(treatment =='untreated')%>% 
    filter(!orig.ident == 'XYY2')
pdc_pair_df <- data.frame(sample = tmp1$orig.ident, subtype= tmp1$subtype, 
                             before=tmp2$Ratio, after=tmp1$Ratio )
ggpaired(pdc_pair_df, cond1 = 'before', cond2 = 'after',
         fill  = "condition", line.color = "gray", line.size = 0.4,
         palette = "npg") +  stat_compare_means(paired = TRUE, method = 't.test',label.x = 1.4)+
    ylab('Prolife Plasma ratio') + facet_wrap(~seurat_clusters,scales= "free" )







# 
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

DoHeatmap2(pdc_filter,marker_list = c('APBBIIP','ABCE1','A1BG','COMMD1'))


library('igraph')
int_df <- int_df[int_df$importance > 10,]
links <- int_df[,-1]
nodes <- rbind(data.frame(gene = unique(links$TF),  attr = rep('TF',length(unique(links$TF)))),
              data.frame(gene = unique(links$target),  attr = rep('target',length(unique(links$target))))
              ) 
nodes <- nodes[!duplicated(nodes[,1]),]
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T) 
l <- layout_with_fr(net)

deg <- degree(net, mode="all")
V(net)$deg <- deg
V(net)$size <- log2(deg) + 1
sub_net <- subgraph(net, V(net)$deg > 2)
# V(net)$color <- ifelse(answers[V(net), 2] == 1, "blue", "red")
# E(net)$width <- E(net)$/6
plot(sub_net, layout=layout_with_fr,
     edge.arrow.size=.1,edge.color="#FFE4B5",
     vertex.size=V(net)$size ,vertex.label.cex=0.5,
     vertex.color="#40E0D0", vertex.label.color="#696969",vertex.frame.color="#ffffff")

