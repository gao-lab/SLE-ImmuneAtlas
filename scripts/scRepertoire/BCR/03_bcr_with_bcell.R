library(RColorBrewer)
library(scRepertoire)
library(Seurat)
library(tidyverse)

# --------------------------- Read and Filter ----------------------------------
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
                                  cloneCall="gene+nt", group.by = "sample", proportion = FALSE, 
                                  cloneTypes=c(Single=1, Small=10, Medium=10, Large=1000))
bcr_df_filter <- bcr_df %>% filter(barcode %in% Cells(bcell_bcr))

quantContig(bcell_bcr, cloneCall="gene+nt", scale = T, chain = "both",split.by  = 'orig.ident',) + 
    stat_compare_means(comparisons =  list(c("treated", "untreated"),c('untreated','HC'),c('treated','HC')),  method = "t.test",
                       bracket.size = 0.5)

clonalDiversity(bcell_bcr, cloneCall = "gene+nt",group.by = 'orig.ident') # not difference
clonalDiversity(bcell_bcr, cloneCall = "gene+nt") # not difference
# useless because different sample has different number of BCR
# clonalProportion(combined_bcr, cloneCall = "gene+nt",split = c(10, 100, 10000), exportTable = F)


# ---------------- overlap in B cell subtype between diff groups ---------------
Idents(bcell_bcr) <- 'treatment'
hc_bcr <- subset(bcell_bcr,idents = 'HC')
Idents(hc_bcr) <- 'subtype'
clonalOverlap2(hc_bcr, cloneCall="gene+nt", method="jaccard",title = 'HC',limit = c(0,0.045)) 
scale_fill_gradient(low = "white",high = "red", na.value = "white",limit = c(0,0.05),space = "Lab",name = "Jaccard Index")

untreat_bcr <- subset(bcell_bcr,idents = 'untreated')
Idents(untreat_bcr) <- 'subtype'
clonalOverlap2(untreat_bcr, cloneCall="gene+nt", method="jaccard",title = 'untreat_bcr',limit = c(0,0.045))
clonalOverlap(untreat_bcr,cloneCall="gene+nt", method="jaccard") + 
    scale_fill_gradient(low = "white",high = "red", na.value = "white",limit = c(0,0.05),space = "Lab",name = "Jaccard Index")


treat_bcr <- subset(bcell_bcr,idents = 'treated')
Idents(treat_bcr) <- 'subtype'
clonalOverlap2(treat_bcr, cloneCall="gene+nt", method="jaccard",title = 'untreat_bcr',limit = c(0,0.045))
clonalOverlap(treat_bcr,cloneCall="gene+nt", method="jaccard") + 
    scale_fill_gradient(low = "blue",high = "red", na.value = "white",limit = c(0,0.05),space = "Lab",name = "Jaccard Index")


# -------------------------- VDJ connection preference -------------------------
tmp <- bcr_df_filter[(grepl('IGHV3-23',bcr_df_filter$IGH) + grepl('IGHJ4',bcr_df_filter$IGH)) == 2,]
table(tmp$sample)
rm(tmp)
bcr_df_filter.list <- split(bcr_df_filter, f=   bcr_df_filter$sample)
bcr_df_filter.list <- addVariable(bcr_df_filter.list, name = "group", 
                             variables = c("before", "before","before","after", "before", "before", "before", "after", "before", "HC", "before",
                                           "before", "after","before", "before", "after", "after", "HC", "HC","before", "after", "HC"))

bcr_df_filter.list.hc <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("HC"))
vizGenes2(bcr_df_filter.list.hc, gene = 'V',chain = 'IGH', y.axis = "J")
vizGenes2(bcr_df_filter.list.hc, gene = "V", chain = "IGL", y.axis = "J") 

bcr_df_filter.list.untreated <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("before"))
vizGenes2(bcr_df_filter.list.untreated, gene = "V", chain = "IGH", y.axis = "J")
vizGenes2(bcr_df_filter.list.untreated, gene = "V", chain = "IGL", y.axis = "J")

bcr_df_filter.list.treated <-  subsetContig(bcr_df_filter.list, name = "group", variables = c("after"))
vizGenes2(bcr_df_filter.list.treated, gene = "V", chain = "IGH", y.axis = "J")
vizGenes2(bcr_df_filter.list.treated, gene = "V", chain = "IGL", y.axis = "J") 


# ------------------------ show paired sample correlation ----------------------
# 
# not use combined_bcr and use bcr_df_filter.list,
debug(scatterClonotype)
undebug(scatterClonotype)
# do not use this 
if(F){
    tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = "LL", 
                            y.axis = "LL2",
                            dot.size = "total",
                            graph = "count",exportTable = F)
    
    tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = "ZPP", 
                            y.axis = "ZPP2",
                            dot.size = "total",
                            graph = "count",exportTable = F)
    tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = "WYF", 
                            y.axis = "WYF2",
                            dot.size = "total",
                            graph = "count",exportTable = F)
    tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = "HXR", 
                            y.axis = "HXR2",
                            dot.size = "total",
                            graph = "count",exportTable = F)
    
    tmp <- scatterClonotype(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = "XYY", 
                            y.axis = "XYY2",
                            dot.size = "total",
                            graph = "count",exportTable = F)
    
}


p <- list()
paired_samples <- c('XYY','HXR','WYF','LL','ZPP')
for(paired_sample in paired_samples){
    paired_sample2 <-paste0(paired_sample,'2')
    tmp <- treat_overlap(bcr_df_filter.list, cloneCall ="gene+nt",
                            x.axis = paired_sample, 
                            y.axis = paired_sample2,
                            dot.size = "total")
    single_not_shared <- (tmp[[paired_sample]] ==1 & tmp[[paired_sample2]] ==0) %>% sum() +
        (tmp[[paired_sample]] ==0 & tmp[[paired_sample2]] ==1) %>% sum()
    single_shared <- (tmp[[paired_sample]]==1 & tmp[[paired_sample2]] ==1) %>% sum()
    expand_not_shared <- (tmp[[paired_sample]] >1 & tmp[[paired_sample2]] <2) %>% sum() +
        (tmp[[paired_sample]] <2 & tmp[[paired_sample2]] >1) %>% sum()
    expand_shared <-  (tmp[[paired_sample]] >1 & tmp[[paired_sample2]] >1) %>% sum()
    plot_df <- data.frame(matrix(NA,ncol = 2,nrow = 2,dimnames = list(c('unexpand','expand'), c('unshare', 'share'))))
    plot_df[1,] <- c(single_not_shared,single_shared);plot_df[2,] <- c(expand_not_shared,expand_shared)
    plot<- plot_df %>% as.matrix() %>% reshape2::melt() %>% 
        ggplot( aes(x=Var1, y=Var2, fill =value),) +   geom_tile() +  
        geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
        scale_fill_viridis( na.value = "white",space = "Lab",name = "cells")+
        theme_classic() +  coord_fixed()
    p[[paired_sample]] <- plot
}
do.call(ggarrange, p)

