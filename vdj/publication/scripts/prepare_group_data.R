library(scRepertoire)
# prepare the BCR data for circle plot

# read in data
bcr_sample <- list.files('./data/10x_bcr/',,pattern = 'csv$')

bcr_list <- list()
for(i in bcr_sample){
    tmp <- as.data.frame(read.csv(paste0('./data/10x_bcr/',i)))
    tmp$sample <- str_split_fixed(i,pattern = '_',2)[,1]
    assign( i, tmp)
    rm(tmp)
}

contig_list_bcr <- list(GW_bcr.csv,GZR_bcr.csv,HXR_bcr.csv,HXR2_bcr.csv,HXX_bcr.csv,
                        LGY_bcr.csv,LL_bcr.csv,LL2_bcr.csv,MXY_bcr.csv,QJY_bcr.csv,
                        SQ_bcr.csv,WYF_bcr.csv, WYF2_bcr.csv,WYY_bcr.csv,XH_bcr.csv,
                        XYY_bcr.csv,XYY2_bcr.csv,ZH_bcr.csv,ZMY1_bcr.csv,ZPP_bcr.csv,
                        ZPP2_bcr.csv,ZS_bcr.csv)
contig_list_bcr <- rbindlist(contig_list_bcr)
head(contig_list_bcr)

# add the meta
meta <- unique(pbmc_all@meta.data[,c('treatment','orig.ident', 'group')]) %>% remove_rownames() %>%column_to_rownames('orig.ident')
contig_list_bcr <-left_join(contig_list_bcr,rownames_to_column(meta),
                                by = c('sample' = 'rowname')) 
contig_list_bcr$barcode <- paste0(contig_list_bcr$sample,'_',contig_list_bcr$group,'_',contig_list_bcr$barcode)

# filter the data
contig_barcode <- Cells(bcell_filter)

# b_cell_merge <- RenameCells(b_cell_merge, new.names = contig_barcode$contig_barcode)
intersect(contig_barcode,contig_list_bcr$barcode) %>% length()

contig_list_bcr.filter <- contig_list_bcr[contig_list_bcr$barcode %in%  contig_barcode,]

# group by group and chain
for (  i in contig_list_bcr.filter$treatment %>% unique() ){
    for ( bcr_chain in  contig_list_bcr.filter$chain %>% unique()){
        print( i);print(bcr_chain)
        tmp <- contig_list_bcr.filter %>% filter(treatment ==  i & chain == bcr_chain)
        print(dim(tmp))
        write_csv(tmp, file = paste0('./vdj/publication/input/', i,"_",bcr_chain,'_bcr.csv'),col_names = T)
    }
}

# group only by group
for ( i in contig_list_bcr.filter$treatment %>% unique() ){
    print( i)
    tmp <- contig_list_bcr.filter %>% filter(treatment ==  i )
    print(dim(tmp))
    write_csv(tmp, file = paste0('./vdj/publication/input/', i,'_bcr.csv'),col_names = T)
    
}
