# 2022.04.24
# consider analysis Plasma with B cell and just divided into plasma, plasmablast

# --------------------------- Read and Filter ----------------------------------
load('./final/seurat/plasma/03-plasma_anno_filter_No_harm.rdata')
load('./final/seurat/b_cell/03-b_cell_anno_filter_harm.rdata')
load('./final/scRepertoire/BCR/combined_bcr.rdata')

bcr_df <- rbindlist(combined_bcr)

# plasma cells
plasma_filter %>% Cells() %>% head()
scRepertoire_barcode <- plasma_filter@meta.data %>% rownames_to_column('barcode_bcr')%>% 
    mutate(new_barcode = str_split_fixed(barcode_bcr ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 
plasma_filter_bcr <- RenameCells(plasma_filter, new.names = scRepertoire_barcode$scRepertoire)


# bcell(not need change seurat cells bracode because done )
# bcell_filter

all_b <- merge(x= bcell_filter, y = plasma_filter_bcr)
intersect(Cells(all_b),bcr_df$barcode) %>% length()

all_b_bcr <- combineExpression(combined_bcr, all_b, 
                               cloneCall="gene+nt", group.by = "sample", proportion = FALSE, 
                               cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000))
bcr_df_filter <- bcr_df %>% filter(barcode %in% Cells(all_b_bcr))
save(all_b_bcr,file = 'final/scRepertoire/BCR/all_b_bcr.rdata')


# some check
quantContig(all_b_bcr, cloneCall="gene+nt", scale = T, chain = "both",split.by  = 'treatment',) + 
    stat_compare_means(comparisons =  list(c("treated", "untreated"),c('untreated','HC'),c('treated','HC')),  method = "t.test",
                       bracket.size = 0.5)
# clonalDiversity(plasma_filter_bcr, cloneCall = "gene+nt",group.by = 'orig.ident') # not difference
# clonalDiversity(plasma_filter_bcr, cloneCall = "gene+nt")

# -------- overlap in B cell and Plamsa subtype between diff treatment ---------
Idents(all_b_bcr) <- 'treatment'
hc_bcr <- subset(all_b_bcr,idents = 'HC')
Idents(hc_bcr) <- 'subtype'
clonalOverlap2(hc_bcr, cloneCall="gene+nt", method="jaccard",title = 'Health Control',limit = c(0,0.045)) 
clonalOverlap2(hc_bcr, cloneCall="gene+nt", method="overlap",title = 'Health Control',limit = c(0,0.2)) 

untreat_bcr <- subset(all_b_bcr,idents = 'untreated')
Idents(untreat_bcr) <- 'subtype'
clonalOverlap2(untreat_bcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Untreated',limit = c(0,0.045))
clonalOverlap2(untreat_bcr, cloneCall="gene+nt", method="overlap",title = 'SLE Untreated',limit = c(0,0.2)) 

treat_bcr <- subset(all_b_bcr,idents = 'treated')
Idents(treat_bcr) <- 'subtype'
clonalOverlap2(treat_bcr, cloneCall="gene+nt", method="jaccard",title = 'SLE Treated',limit = c(0,0.045))
clonalOverlap2(treat_bcr, cloneCall="gene+nt", method="overlap",title = 'SLE Treated',limit = c(0,0.2))

rm(hc_bcr,untreat_bcr,treat_bcr);gc()
