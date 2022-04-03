library(RColorBrewer)
library(scRepertoire)
library(Seurat)
library(tidyverse)

bcr_df <- rbindlist(combined_bcr)
# load the bcell seurat object and combined_bcr
# bcell_bcr <- combineExpression(combined_bcr, bcell_filter, 
#                             cloneCall="aa", group.by = "sample", proportion = FALSE, 
#                             cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
scRepertoire_barcode <- bcell_filter@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 

bcell_filter <- RenameCells(bcell_filter, new.names = scRepertoire_barcode$scRepertoire)
intersect(Cells(bcell_filter),bcr_df$barcode) %>% length()

# debug(bcell_bcr <- combineExpression(combined_bcr, bcell_filter, 
#                                      cloneCall="aa", group.by = "sample", proportion = FALSE, 
#                                      cloneTypes=c(Single=1, Small=10, Medium=10, Large=1000)))
bcell_bcr <- combineExpression(combined_bcr, bcell_filter, 
                                  cloneCall="aa", group.by = "sample", proportion = FALSE, 
                                  cloneTypes=c(Single=1, Small=10, Medium=10, Large=1000))


quantContig(bcell_bcr, cloneCall="aa", scale = T, chain = "both",split.by  = 'orig.ident',) + 
    stat_compare_means(comparisons =  list(c("treated", "untreated"),c('untreated','HC'),c('treated','HC')),  method = "t.test",
                       bracket.size = 0.5)

clonalDiversity(bcell_bcr, cloneCall = "aa",group.by = 'orig.ident') # not difference
clonalDiversity(bcell_bcr, cloneCall = "aa") # not difference
clonalProportion(bcell_bcr, cloneCall = "aa",split  = 'orig.ident', exportTable = T)



# --------------- overlap in B cell subtype between diff groups ---------------
Idents(bcell_bcr) <- 'treatment'
hc_bcr <- subset(bcell_bcr,idents = 'HC')
Idents(hc_bcr) <- 'subtype'
hc_bcr <- subset(hc_bcr, idents = 'B.IFN-response',invert =T)
clonalOverlap(hc_bcr, cloneCall="aa", method="jaccard",exportTable = T) %>% reshape2::melt('names') %>% 
ggplot(aes(x = names,y = variable,fill = value))+
    geom_tile()+theme_bw()+
    theme_minimal()+ # 设置主题为无边框
    scale_fill_gradient(low = "white",high = "red", na.value = "white",limit = c(0,0.005),space = "Lab",name = "Jaccard Index")+
    labs(title = "Correlation Heatmap")+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,size = 10,hjust = 1),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 12),
          axis.title = element_blank(), # 去除横纵坐标轴标题标签
          panel.grid.major = element_blank(), # 去除背景网格线
          panel.background = element_blank(), # 去除背景颜色
          legend.justification = c(1,0),
          legend.position = c(1.3,0.3), # 图例位置
          legend.direction = "horizontal")+ # vertical(垂直)/horizontal(水平)
    coord_fixed()+ # 确保x轴一个单位与y轴一个单位长度相同
    guides(fill = guide_colorbar(barwidth = 8,barheight = 1.5, # 图例长宽
                                 title.position = "top",title.hjust = 0.5))
untreat_bcr <- subset(bcell_bcr,idents = 'untreated')
Idents(untreat_bcr) <- 'subtype'
# untreat_bcr <- subset(untreat_bcr, idents = 'B.IFN-response',invert =T)
clonalOverlap(untreat_bcr, cloneCall="aa", method="jaccard",exportTable = T) %>% reshape2::melt('names') %>% 
    ggplot(aes(x = names,y = variable,fill = value))+
    geom_tile()+theme_bw()+
    theme_minimal()+ # 设置主题为无边框
    scale_fill_gradient(low = "white",high = "red", na.value = "white",limit = c(0,0.005),space = "Lab",name = "Jaccard Index")+
    labs(title = "Correlation Heatmap")+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,size = 10,hjust = 1),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 12),
          axis.title = element_blank(), # 去除横纵坐标轴标题标签
          panel.grid.major = element_blank(), # 去除背景网格线
          panel.background = element_blank(), # 去除背景颜色
          legend.justification = c(1,0),
          legend.position = c(1.3,0.3), # 图例位置
          legend.direction = "horizontal")+ # vertical(垂直)/horizontal(水平)
    coord_fixed()+ # 确保x轴一个单位与y轴一个单位长度相同
    guides(fill = guide_colorbar(barwidth = 8,barheight = 1.5, # 图例长宽
                                 title.position = "top",title.hjust = 0.5))

treat_bcr <- subset(bcell_bcr,idents = 'treated')
Idents(treat_bcr) <- 'subtype'
clonalOverlap(treat_bcr, cloneCall="aa", method="jaccard",exportTable = T) %>% reshape2::melt('names') %>% 
    ggplot(aes(x = names,y = variable,fill = value))+
    geom_tile()+theme_bw()+
    theme_minimal()+ # 设置主题为无边框
    scale_fill_gradient(low = "white",high = "red", na.value = "white",limit = c(0,0.005),space = "Lab",name = "Jaccard Index")+
    labs(title = "Correlation Heatmap")+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,size = 10,hjust = 1),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 12),
          axis.title = element_blank(), # 去除横纵坐标轴标题标签
          panel.grid.major = element_blank(), # 去除背景网格线
          panel.background = element_blank(), # 去除背景颜色
          legend.justification = c(1,0),
          legend.position = c(1.3,0.3), # 图例位置
          legend.direction = "horizontal")+ # vertical(垂直)/horizontal(水平)
    coord_fixed()+ # 确保x轴一个单位与y轴一个单位长度相同
    guides(fill = guide_colorbar(barwidth = 8,barheight = 1.5, # 图例长宽
                                 title.position = "top",title.hjust = 0.5))

# -------------------------- VDJ connection preference -------------------------
bcr_df_filter <- bcr_df %>% filter(barcode %in% Cells(bcell_bcr))
bcr_df_filter.list <- split(bcr_df_filter, f=   bcr_df_filter$sample)
bcr_df_filter.list <- addVariable(bcr_df_filter.list, name = "group", 
                             variables = c("before", "before","before","after", "before", "before", "before", "after", "before", "HC", "before",
                                           "before", "after","before", "before", "after", "after", "HC", "HC","before", "after", "HC"))

bcr_df_filter.list.hc <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("HC"))
vizGenes(bcr_df_filter.list.hc, gene = "V", chain = "IGH", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")
vizGenes(bcr_df_filter.list.hc, gene = "V", chain = "IGL", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")

bcr_df_filter.list.untreated <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("before"))
vizGenes(bcr_df_filter.list.untreated, gene = "V", chain = "IGH", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")
vizGenes(bcr_df_filter.list.untreated, gene = "V", chain = "IGL", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")

bcr_df_filter.list.treated <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("after"))
vizGenes(bcr_df_filter.list.treated, gene = "V", chain = "IGH", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")
vizGenes(bcr_df_filter.list.treated, gene = "V", chain = "IGL", y.axis = "J", 
         plot = "heatmap", scale = TRUE, order = "gene") +  
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_gradient(low = "#C0C0C0",high = "#FF0000", na.value = "white",limit = c(0,0.08),space = "Lab",name = "Preference")


# ------------------------ paired --------------------------------------------
# not use combined_bcr and use bcr_df_filter.list,
debug(scatterClonotype)
undebug(scatterClonotype)
tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="aa",
                 x.axis = "LL", 
                 y.axis = "LL2",
                 dot.size = "total",
                 graph = "count",exportTable = T)

tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="aa",
                 x.axis = "ZPP", 
                 y.axis = "ZPP2",
                 dot.size = "total",
                 graph = "count",exportTable = T)
tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="aa",
                 x.axis = "WYF", 
                 y.axis = "WYF2",
                 dot.size = "total",
                 graph = "count",exportTable = T)
tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="aa",
                 x.axis = "HXR", 
                 y.axis = "HXR2",
                 dot.size = "total",
                 graph = "count",exportTable = T)

tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="aa",
                        x.axis = "XYY", 
                        y.axis = "XYY2",
                        dot.size = "total",
                        graph = "count",exportTable = T)

tmp <- data.frame(matrix(NA,ncol = 2,nrow = 2,dimnames = list(c('unexpand','expand'), c('unshare', 'share'))))
tmp[1,] <- c(775,0);tmp[2,] <- c(2,0)
pheatmap(tmp,display_numbers = T,number_color='white',fontsize_number = 15,number_format = '%.0f',
         border_color ='black', color = c('#0000CD','red'),
         cluster_rows =F,cluster_cols = F,scale = 'col')
