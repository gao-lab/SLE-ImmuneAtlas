library(tidyverse)
library(scRepertoire)
setwd('/data/sle/')

################################################################################
#
# Create seurat object 
#
################################################################################
# output is in final/scRepertoire/TCR/

load('./final/scRepertoire/TCR/combined_tcr.rdata')
combined_tcr.df <- rbindlist(combined_tcr)

# -------------- bulid seurat object with TCR -------------
scRepertoire_barcode <- cd8_filter@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 

cd8_tcr <- RenameCells(cd8_filter, new.names = scRepertoire_barcode$scRepertoire)
# cd8_filter : 62491 
# combined_tcr.df: 151666
# intersect 55795 
intersect(Cells(cd8_tcr),combined_tcr.df$barcode) %>% length()

cd8_tcr <- combineExpression(combined_tcr, cd8_tcr,  cloneCall="gene+nt", 
                               group.by = "sample", proportion = FALSE, 
                               cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000, Super = 5000))
# modify show order
cd8_tcr$cloneType <- factor(cd8_tcr$cloneType, 
                            levels = c("Large (100 < X <= 1000)", "Medium (10 < X <= 100)", 
                                       "Small (1 < X <= 10)", "Single (0 < X <= 1)", NA))
save(cd8_tcr, file = 'final/scRepertoire/TCR/cd8_tcr.rdata')

#------- build scRepertoire object with filter TCR -------
combined_tcr.df.filter <-  filter(combined_tcr.df, barcode %in% Cells(cd8_tcr))
combined_tcr.filter <- split(combined_tcr.df.filter, f = combined_tcr.df.filter$sample)
combined_tcr.filter <- addVariable(combined_tcr.filter, name = "group", 
                                  variables = c("before", "before","before","after", "before", "before", "before", "after", "before", "HC", "before",
                                                "before", "after","before", "before", "after", "after", "HC", "HC","before", "after", "HC"))


################################################################################
#
# Clone expansion analysis 
#
################################################################################
table(cd8_tcr$cloneType)

DimPlot(cd8_tcr, group.by = "cloneType",pt.size = 0.25)  + 
    scale_color_brewer(palette = "RdBu", direction = 1) + NoAxes()
DimPlot(cd8_tcr, group.by = "cloneType", split.by = 'treatment',pt.size = 0.5)  + 
    scale_color_brewer(palette = "RdBu", direction = 1) + NoAxes()
occupiedscRepertoire(cd8_tcr, x.axis = "subtype", facet.by = 'treatment',label = F,proportion = F) + 
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
Idents(cd8_tcr) <- 'treatment'
hc_cd8_tcr <- subset(cd8_tcr,idents = 'HC'); Idents(hc_cd8_tcr) <- 'subtype'
untreat_cd8_tcr <- subset(cd8_tcr,idents = 'untreated'); Idents(untreat_cd8_tcr) <- 'subtype'
treat_cd8_tcr <- subset(cd8_tcr,idents = 'treated'); Idents(treat_cd8_tcr) <- 'subtype'

clonalOverlap2(hc_cd8_tcr, cloneCall="gene+nt", method="jaccard",title = 'Health Control',limit = c(0,0.2)) 
clonalOverlap2(hc_cd8_tcr, cloneCall="gene+nt", method="overlap",title = 'Health Control',limit = c(0,0.6)) 

clonalOverlap2(untreat_cd8_tcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Untreated',limit = c(0,0.2))
clonalOverlap2(untreat_cd8_tcr, cloneCall="gene+nt", method="overlap",title = 'SLE Untreated',limit = c(0,0.6)) 

clonalOverlap2(treat_cd8_tcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Treated',limit = c(0,0.2))
clonalOverlap2(treat_cd8_tcr, cloneCall="gene+nt", method="overlap",title = 'SLE Treated',limit = c(0,0.6))

rm(hc_cd8_tcr,untreat_cd8_tcr,treat_cd8_tcr);gc()

# sample overlap before and after treatment
scatterClonotype(combined_tcr.filter, cloneCall ="gene+nt", x.axis = "ZPP", y.axis = "ZPP2",
                 dot.size = "total",graph = "count")


