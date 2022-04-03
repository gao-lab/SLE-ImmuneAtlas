
setwd('/data/sle')
source('./scripts/function_R/utils.R')
library(reshape2)
library(pheatmap)

bk =  c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
# ----------------------------------- IFN --------------------------------------
allegs = get('GO:0034340', org.Hs.egGO2ALLEGS)
genes = unlist(mget(allegs,org.Hs.egSYMBOL)) %>% unique()
# IFN_genes
IFN_genes.2 <- intersect(genes, row.names(pbmc_final))
pbmc_final <- AddModuleScore(pbmc_final,features = list(IFN_genes.2),name = 'IFN_score')

IFN_logFC <- c()


pbmc_final@meta.data %>% 
    select('disease','treatment','subtype','orig.ident','IFN_score1') %>% 
    group_by(treatment,subtype) %>% summarise(mean_IFN = mean(IFN_score1)) %>%
    spread(treatment,mean_IFN) %>% mutate(treated_HC = log2(treated/HC), 
                                          untreated_HC = log2(untreated/HC) ,
                                          untreated_treated = log2(untreated/treated)) %>%
    column_to_rownames('subtype') %>% select(4,5,6) %>%
    pheatmap(scale ='none',cluster_rows =F, cluster_cols = F,color =  c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
             breaks = bk)

pbmc_final@meta.data %>% 
    select('disease','treatment','subtype','orig.ident','IFN_score1') %>% 
    group_by(treatment,subtype) %>% summarise(mean_IFN= mean(IFN_score1)) %>%
    spread(treatment,mean_IFN) %>% mutate(treated_HC = 2^(treated)-2^(HC), 
                                          untreated_HC = 2^(untreated)-2^(HC) ,
                                          untreated_treated = 2^(untreated)-2^(treated)) %>%
    column_to_rownames('subtype') %>% 
    # select(1,2,3) %>%
    select(5,4) %>%
    pheatmap(scale ='none',cluster_rows =F, cluster_cols = F, color = c(colorRampPalette(colors = c("blue",'white'))(25),
                                                                        colorRampPalette(colors = c('white',"red"))(25)),
             breaks = c(seq(-0.25,0.25,by=0.01)))

# ----------------------------------- JAK --------------------------------------
bk.2 =  c(seq(0,5,by=0.5)
JAK_genes.2 <- intersect(JAK_genes, row.names(pbmc_final))
pbmc_final <- AddModuleScore(pbmc_final,features = list(JAK_genes.2),name = 'JAK_score')

pbmc_final@meta.data %>% 
    select('disease','treatment','subtype','orig.ident','JAK_score1') %>% 
    group_by(treatment,subtype) %>% summarise(mean_IFN = mean(JAK_score1)) %>%
    spread(treatment,mean_IFN) %>% mutate(treated_HC = 2^(treated)-2^(HC), 
                                          untreated_HC = 2^(untreated)-2^(HC) ,
                                          untreated_treated = 2^(untreated)-2^(treated)) %>%
    column_to_rownames('subtype') %>% 
    select(1,2,3) %>%
    # select(4,5,6) %>%
    pheatmap(scale ='none',cluster_rows =F, cluster_cols = F, color = c(colorRampPalette(colors = c("blue",'white'))(11),
                                                                        colorRampPalette(colors = c('white',"red"))(15)),
             breaks = c(seq(-0.1,0.15,by=0.01)))

pbmc_final@meta.data %>% 
    select('disease','treatment','subtype','orig.ident','JAK_score1') %>% 
    group_by(treatment,subtype) %>% summarise(mean_JAK= mean(JAK_score1)) %>%
    spread(treatment,mean_JAK) %>% mutate(treated_HC = 2^(treated)-2^(HC), 
                                          untreated_HC = 2^(untreated)-2^(HC) ,
                                          untreated_treated = 2^(untreated)-2^(treated)) %>%
    column_to_rownames('subtype') %>% 
    # select(1,2,3) %>%
    select(5,4) %>%
    pheatmap(scale ='none',cluster_rows =F, cluster_cols = F, color = c(colorRampPalette(colors = c("blue",'white'))(7),
                                                                        colorRampPalette(colors = c('white',"red"))(7)),
             breaks = c(seq(-0.06,0.06,by=0.01)))

# ------------------------------------- mTOR -----------------------------------
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
mTOR_gene <- getBM(attributes=c('hgnc_symbol'),
                   filters = 'go', values = 'GO:0003723', mart = ensembl)

mTOR_gene.2 <- intersect(mTOR_gene$hgnc_symbol, row.names(pbmc_final))
pbmc_final <- AddModuleScore(pbmc_final,features = list(mTOR_gene.2),name = 'mTOR_score')

pbmc_final@meta.data %>% 
    dplyr::select('disease','treatment','subtype','orig.ident','mTOR_score1') %>% 
    group_by(treatment,subtype) %>% summarise(mean_mTOR= mean(mTOR_score1)) %>%
    spread(treatment,mean_mTOR) %>% mutate(treated_HC = 2^(treated)-2^(HC), 
                                          untreated_HC = 2^(untreated)-2^(HC) ,
                                          untreated_treated = 2^(untreated)-2^(treated)) %>%
    column_to_rownames('subtype') %>% 
    # select(1,2,3) %>%
    dplyr::select(5,4) %>%
    pheatmap(scale ='none',cluster_rows =F, cluster_cols = F, color = c(colorRampPalette(colors = c("blue",'white'))(10),
                                                                        colorRampPalette(colors = c('white',"red"))(10)),
             breaks = c(seq(-0.1,0.1,by=0.01)))

