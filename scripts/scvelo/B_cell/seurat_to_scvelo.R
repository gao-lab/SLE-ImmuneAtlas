
load('./final/seurat/b_cell/03-b_cell_anno_filter_harm.rdata') #name 'bcell_filter'

bcell_filter_meta <- merge(bcell_filter@meta.data, Embeddings(bcell_filter,'umap'),by.x = 0, by.y=0,all =T)

# I change seurat barcode in the scRepertoire,here need to recover it
sample_index.df <- cd4_filter@meta.data %>% rownames_to_column('rowname') %>% 
    select( rowname,orig.ident) %>% separate('rowname',c('front','index'),'_') %>%
    select(index,orig.ident) %>% distinct()

bcell_filter_meta.index <- str_split_fixed(bcell_filter_meta$Row.names,'_',3)[,1]  %>% 
    data.frame(sample = .) %>%
    left_join(sample_index.df ,by = c('sample' = 'orig.ident'))

bcell_filter_meta$Row.names <- paste0(str_split_fixed(bcell_filter_meta$Row.names,'_',3)[,3] , '_' , bcell_filter_meta.index$index)


write.csv(bcell_filter_meta, file = './scripts/scvelo/B_cell/bcell_filter_meta.csv')
rm(bcell_filter_meta)
