library(tidyverse)
library(scRepertoire)
setwd('/data/sle/')
source('./scripts/function_R/utils.R')

################################################################################
#
# Create seurat object 
#
################################################################################
# NOTE: we use some shared object in the 03_CD8_TCR.R (combined_tcr.df)
load('./final/seurat/t_cell/04-CD4_Tcell_filter_anno.rdata')

# -------------- build seurat object with TCR -------------
scRepertoire_barcode <- cd4_filter@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 

cd4_tcr <- RenameCells(cd4_filter, new.names = scRepertoire_barcode$scRepertoire)
# cd4_filter : 52169 
# combined_tcr.df: 151666
# intersect 49375 
intersect(Cells(cd4_tcr),combined_tcr.df$barcode) %>% length()

cd4_tcr <- combineExpression(combined_tcr, cd4_tcr, cloneCall="gene+nt", 
                             group.by = "sample", proportion = FALSE, 
                             cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000, Super = 5000))
# modify show order
cd4_tcr$cloneType <- factor(cd4_tcr$cloneType, 
                            levels = c("Large (100 < X <= 1000)", "Medium (10 < X <= 100)", 
                                       "Small (1 < X <= 10)", "Single (0 < X <= 1)", NA))
save(cd4_tcr, file = 'final/scRepertoire/TCR/cd4_tcr.rdata')

load('final/scRepertoire/TCR/cd4_tcr.rdata')
################################################################################
#
# Clone expansion analysis 
#
################################################################################
table(cd4_tcr$cloneType)

DotPlot(cd4_filter,features = c('CD8A','CD8B','CD4'))
DimPlot(cd4_tcr, group.by = "cloneType",pt.size = 0.25)  + 
    scale_color_brewer(palette = "RdBu", direction = 1) + NoAxes()
DimPlot(cd4_tcr, group.by = "cloneType", split.by = 'treatment',pt.size = 0.5)  + 
    scale_color_brewer(palette = "RdBu", direction = 1) + NoAxes()
occupiedscRepertoire(cd4_tcr, x.axis = "subtype", facet.by = 'treatment',label = F,proportion = F) + 
    scale_fill_brewer(palette = 'RdBu') +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    labs(fill = "Clone Size")

# clone size corresponding occupied space(we need subset the subtype)
clonalProportion(combined_tcr.filter, cloneCall = "gene+nt",split = c(10, 100, 1000, 100000)) + 
    scale_fill_brewer(palette = 'RdBu',name = "Clone Size", labels = c("Top10", "Top100", "Top1000",'Other')) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab('')




################################################################################
#
# Clone overlap  analysis 
#
################################################################################
# subtype overlap 
Idents(cd4_tcr) <- 'treatment'
hc_cd4_tcr <- subset(cd4_tcr,idents = 'HC'); Idents(hc_cd4_tcr) <- 'subtype'
untreat_cd4_tcr <- subset(cd4_tcr,idents = 'untreated'); Idents(untreat_cd4_tcr) <- 'subtype'
treat_cd4_tcr <- subset(cd4_tcr,idents = 'treated'); Idents(treat_cd4_tcr) <- 'subtype'

clonalOverlap2(hc_cd4_tcr, cloneCall="gene+nt", method="jaccard",title = 'Health Control',limit = c(0,0.012)) 
clonalOverlap2(hc_cd4_tcr, cloneCall="gene+nt", method="overlap",title = 'Health Control',limit = c(0,0.6)) 

clonalOverlap2(untreat_cd4_tcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Untreated',limit = c(0,0.2))
clonalOverlap2(untreat_cd4_tcr, cloneCall="gene+nt", method="overlap",title = 'SLE Untreated',limit = c(0,0.6)) 

clonalOverlap2(treat_cd4_tcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Treated',limit = c(0,0.2))
clonalOverlap2(treat_cd4_tcr, cloneCall="gene+nt", method="overlap",title = 'SLE Treated',limit = c(0,0.6))

rm(hc_cd4_tcr,untreat_cd4_tcr,treat_cd4_tcr);gc()

# sample overlap before and after treatment
scatterClonotype(combined_tcr.filter, cloneCall ="gene+nt", x.axis = "ZPP", y.axis = "ZPP2",
                 dot.size = "total",graph = "count")





