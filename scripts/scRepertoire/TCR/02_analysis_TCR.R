library(scRepertoire)
library(Seurat)
library(tidyverse)
setwd('/data/sle/')

################################################################################
#
# Create seurat object 
#
################################################################################
# output is in final/scRepertoire/TCR/
# NOTE : only contain ab T cell !

load('./final/scRepertoire/TCR/combined_tcr.rdata')
combined_tcr_filter <- rbindlist(combined_tcr)

scRepertoire_barcode <- t_cell_merge@meta.data %>% rownames_to_column('barcode')%>% 
    mutate(new_barcode = str_split_fixed(barcode ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 

t_cell_merge <- RenameCells(t_cell_merge, new.names = scRepertoire_barcode$scRepertoire)
# t_cell_merge : 145716 
# combined_tcr_filter: 151666
# intersect 114129 
intersect(Cells(t_cell_merge),combined_tcr_filter$barcode) %>% length()

t_cell_AB <- combineExpression(combined_tcr, t_cell_merge,  cloneCall="gene+nt", 
                                  group.by = "sample", proportion = FALSE, 
                                  cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000))
save(t_cell_AB, file = 'final/scRepertoire/TCR/t_cell_AB.rdata')
rm(t_cell_AB)

# we focus on CD8 T cell please see 03_CD8_TCR.R


