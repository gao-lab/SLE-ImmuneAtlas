library(Seurat)
library(tidyverse)
library(scRepertoire)
library(stringr)
library(data.table)
library(magrittr)
library(ggpubr)
library(patchwork)
library(viridis)
library(ggrepel)
set.seed(520)
setwd('/data/sle')
source('scripts/function_R/utils.R')

## Input
# combined_bcr                : combined cellranger output
# meta                        : meta info of sample
# bcell_immacantation_SHM.csv : SHM result from Ig blast
# 10X_clone-pass_germ-pass.tsv: clonotype info result from immacantation
# 04-final-b_cell_filter.rdata: filtered b cell rna data
# seurat_bcell_filter_meta.csv: meta of B cells

## Output 
# Fig.2a,b,e(left); Fig.5a,b,c,d,e; Extend Data Fig.5

################################################################################
#
# Read in data
#
################################################################################
# read in data
bcr_sample <- list.files('data/10x_bcr/',pattern='csv$')
bcr_sample_name <- str_split_fixed(bcr_sample,pattern='_',2)[,1]
load('final/seurat/b_cell/04-final-b_cell_filter.rdata')

bcr_list <- list()
for(i in bcr_sample){
    tmp <- as.data.frame(read.csv(paste0('./data/10x_bcr/',i)))
    tmp$sample <- str_split_fixed(i,pattern='_',2)[,1]
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
meta <- read.csv('data/meta.csv', row.names=1)
contig_list_bcr <-left_join(contig_list_bcr,rownames_to_column(meta),
                            by=c('sample'='rowname')) 
contig_list_bcr$barcode <- paste0(contig_list_bcr$sample,'_',contig_list_bcr$group,'_',contig_list_bcr$barcode)


################################################################################
#
# BCR repertoire analysis via scRepertoire
#
################################################################################

## STEP1: build the combined_bcr object
load('final/scRepertoire/BCR/combined_bcr.rdata')

## STEP2: filter the BCR repertoire
combined_bcr.df <- rbindlist(combined_bcr)
combined_bcr.df.filter <- combined_bcr.df %>% filter(barcode %in% Cells(bcell_filter))
combined_bcr.filter <- split(combined_bcr.df.filter, f=combined_bcr.df.filter$sample)


## STEP3: BCR repertoire diversity(Fig.2e left)
clonalHomeostasis(combined_bcr, cloneCall="gene+nt",
                  cloneTypes=c(Rare=0.005, Small =0.05,
                                 Large=0.1, Hyperexpanded=1)) # just test
tmp <- clonalHomeostasis(combined_bcr, cloneCall="gene+nt",
                         cloneTypes=c(Rare=0.01, Small=0.05,
                                        Large=0.1, Hyperexpanded=1),exportTable=T) %>% 
    melt () %>% left_join(meta %>% rownames_to_column(),by=c('Var1'='rowname' )) %>%
    arrange(treatment,Var1)
# tmp$Var1 %<>% unique() %>% factor(.,levels=.)
tmp$Var1 <- fct_inorder(tmp$Var1)
tmp$sample <- paste0('Sample', c(1:22) %>% rep(4)  %>% sort())
tmp$sample <- fct_inorder(tmp$sample)
ggplot(tmp, aes(x=sample, y=value, fill=Var2)) + 
    geom_bar(stat="identity", position="fill", color="black", lwd=0.25) + 
    scale_fill_manual(name="Clonotype Group", values=c('#7BC6FF','#C7FDED','#FFB2AD','#E7525A') ) + 
    xlab("") +  ylab("Relative Abundance") + theme_classic() + 
    theme(axis.text.x=element_text(angle=30, hjust=1))


################################################################################
#
# B cell and plasma with RNA 
#
################################################################################
### Fig.5a and Fig.5b
### B cell
bcell_bcr<- combineExpression(combined_bcr, bcell_filter, 
                                  cloneCall="gene+nt", group.by="sample", proportion=FALSE, 
                                  cloneTypes=c(Single=1, Small=10, Medium=10, Large=100))
Idents(bcell_bcr) <- 'cloneType'
bcell_bcr.filter <- subset(bcell_bcr, idents=c('Single (0 < X <= 1)','Small (1 < X <= 10)'))
DimPlot(bcell_bcr.filter, group.by="cloneType",cols=c('#D6604D','#92C5DE') %>% rev(), 
        split.by='treatment',pt.size=1.5) + NoAxes()
ggplot(data=bcell_bcr.filter@meta.data, aes(x=treatment, fill=cloneType))+
    geom_bar(stat='count',position='fill') + labs(y='proportions', x="") + 
    scale_fill_discrete(labels= names(table(bcell_bcr.filter$cloneType))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values= c('#92C5DE','#D6604D')) + theme_classic() 

### plasma
load('./final/seurat/plasma/03-plasma_anno_filter_No_harm.rdata')
scRepertoire_barcode <- plasma_filter@meta.data %>% rownames_to_column('barcode_bcr')%>% 
    mutate(new_barcode=str_split_fixed(barcode_bcr ,'_',2)[,1])%>% 
    mutate(scRepertoire=paste0(orig.ident,'_',group,'_',new_barcode) )%>% select(scRepertoire) 
plasma_bcr <- RenameCells(plasma_filter, new.names=scRepertoire_barcode$scRepertoire)
plasma_bcr <- combineExpression(combined_bcr, plasma_bcr, addLabel =T,
                                   cloneCall="gene+nt", group.by="sample", proportion=FALSE, 
                                   cloneTypes=c(Single=1, Small=10, Medium=10, Large=100))
Idents(plasma_bcr) <- 'cloneType'
plasma_bcr.filter <- subset(plasma_bcr, idents=c('Single (0 < X <= 1)','Small (1 < X <= 10)'))
DimPlot(plasma_bcr.filter, group.by="cloneType",cols=c('#92C5DE','#D6604D'),
        split.by='treatment',pt.size=1.5) + NoAxes()
ggplot(data = plasma_bcr.filter@meta.data, aes(x = treatment, fill = cloneType))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(plasma_filter$cloneType))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values= rev(c('#92C5DE','#D6604D'))) + theme_classic() 


################################################################################
#
# BCR antibody isotype ratio
#
################################################################################
## STEP1: read in the data
files <- dir('vdj/immcatation/data/')
meta$name <- rownames(meta)
HC_df <- data.frame(matrix(NA, ncol=dim(tmp)[2] + 1, nrow=0))
SLE_df <- data.frame(matrix(NA, ncol=dim(tmp)[2] + 1, nrow=0))
for (file in files){
    name=str_split(file,'_')[[1]][1]
    tmp <-  read.csv(paste0('vdj/immcatation/data/',file,'/10X_clone-pass_germ-pass.tsv'),header=T,sep='\t')
    tmp %<>% mutate(sample=name)
    if(meta$group[which(meta$name == name)] == 'HC' ){
        HC_df <- rbind(HC_df,tmp)
    }else if(meta$group[which(meta$name == name)] == 'SLE'){
        SLE_df <- rbind(SLE_df, tmp)
    }
}
dim(HC_df)
dim(SLE_df)
SLE_df$disease <- 'SLE'
HC_df$disease <- 'HC'
all_df <-rbind(SLE_df, HC_df)

all_df$barcode_match <- paste0(all_df$sample,'_',all_df$disease,'_',str_split_fixed(all_df$sequence_id, '_',n=3)[,1])
bcell_rna_meta <- read.csv('vdj/immcatation/analysis/seurat_bcell_filter_meta.csv')
plasma_rna_meta <- read.csv('vdj/immcatation/analysis/seurat_plasma_filter_meta.csv')
bcell_rna_meta <- bcell_rna_meta[,c('X','orig.ident','nCount_RNA','nFeature_RNA','group','treatment','pair','percent_mito','percent_ribo','S.Score','G2M.Score','Phase','old.ident','RNA_snn_res.0.8','RNA_snn_res.1','seurat_clusters','main_type','subtype')]
plasma_rna_meta <- plasma_rna_meta[,c('X','orig.ident','nCount_RNA','nFeature_RNA','group','treatment','pair','percent_mito','percent_ribo','S.Score','G2M.Score','Phase','old.ident','RNA_snn_res.0.8','RNA_snn_res.1','seurat_clusters','main_type','subtype')]
rna_meta <- rbind(bcell_rna_meta, plasma_rna_meta)

intersect(all_df$barcode_match, rna_meta$X ) %>% length()

## STEP2: filter via rna
filter_all_df <- all_df[which(all_df$barcode_match %in% rna_meta$X),] 
filter_all_df %<>% left_join( rna_meta, by=c('barcode_match'='X')) %>%
    filter(!c_call %in% c('','IGHG4'))
v_j_combine <- filter_all_df %>% select(c('v_call','j_call','disease','subtype','orig.ident')) %>% mutate(combine=(grepl('IGHV3-23',.[,1]) + grepl('IGHJ4',.[,2])) ) %>% filter(combine == 2)
v_j_combine[1:5,]
filter_all_df_sle <- filter_all_df[which(filter_all_df$disease == 'SLE'),]
filter_all_df_hc <- filter_all_df[which(filter_all_df$disease == 'HC'),]

## STEP3: plot
### total plot
#### sle
sle_c_pie <- table(filter_all_df_sle$c_call , filter_all_df_sle$main_type) %>% data.frame()
sle_c_pie_bcell <- sle_c_pie %>% filter(Var2 == 'Bcell') %>% filter(Var1 != '')
colnames(sle_c_pie_bcell) <- c('C gene','main_type','count')
a <- pie_plot(sle_c_pie_bcell,'count','C gene')

sle_c_pie_plasma <- sle_c_pie %>% filter(Var2 == 'Plasma') %>% filter(Var1 != '')
colnames(sle_c_pie_plasma) <- c('C gene','main_type','count')
b <- pie_plot(sle_c_pie_plasma,'count','C gene')
##### hc
hc_c_pie <- table(filter_all_df_hc$c_call , filter_all_df_hc$main_type) %>% data.frame()
hc_c_pie_bcell <- hc_c_pie %>% filter(Var2 == 'Bcell') %>% filter(Var1 != '')
colnames(hc_c_pie_bcell) <- c('C gene','main_type','count')
c <- pie_plot(hc_c_pie_bcell,'count','C gene')

hc_c_pie_plasma <- hc_c_pie %>% filter(Var2 == 'Plasma') %>% filter(Var1 != '')
colnames(hc_c_pie_plasma) <- c('C gene','main_type','count')
d <- pie_plot(hc_c_pie_plasma,'count','C gene')
options(repr.plot.width=14, repr.plot.height=14)
ggarrange(a,b,c,d)

### by B cell subtype
filter_all_df$subtype[grepl(filter_all_df$subtype,pattern='plasma')] <- 'plasma'
filter_all_df$subtype %>% table()
sub_filter_all_df <- filter_all_df %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) 
p2_list <- list();j=1
for (treatment in c('HC','treated','untreated')) {
    # print(treatment)
    sub <- sub_filter_all_df[which(sub_filter_all_df$treatment == treatment),]
    p_list <- list()
    i=1
    for (celltype in c('B.mem','B.mem.CD27-','B.mem.CXCR3+','B.mem.IGHM+','plasma')){
        # print(celltype)
        sub2 <- sub[which(sub$subtype == celltype),]
        # plot
        plot_tab <- table(sub2$c_call , sub2$subtype) %>% data.frame() %>% filter(Var1 != '')
        colnames(plot_tab) <- c('C gene','subtype','count')
        p <- pie_plot(plot_tab,'count','C gene')
        p_list[[i]] <- p
        # names(p_list[[i]]) <- paste0(treatment,celltype)
        i=i+1
    }
    p2 <- do.call(patchwork::wrap_plots,c(p_list,ncol =1)) 
    p2_list[[j]] <- p2; j=j+1
}

do.call(patchwork::wrap_plots,c(p2_list,ncol =3))

################################################################################
#
# BCR SHM ratio
#
################################################################################
# bcell_immacantation_SHM.csv is the result from SHM analysis
SHM_info <- read.csv('./scripts/immcantaion/bcell_immacantation_SHM.csv')
SHM_info$barcode <- str_split_fixed(SHM_info$sequence_id, pattern='_',n=2)[,1] %>% 
    paste0('_',SHM_info$sample)

bcell_bcr$barcode <- str_split_fixed(Cells(bcell_bcr), pattern='_',n=3)[,3]  %>% 
    paste0('_',bcell_bcr$orig.ident)
# 20892/23322=90%: 90% of b cell detect heavy chain in immcantation results
intersect(SHM_info$barcode , bcell_bcr$barcode) %>% length()

## Add SHM to seurat meta
bcell_bcr@meta.data %<>%  left_join(SHM_info, by=c('barcode'='barcode'))
bcell_bcr$mu_freq_seq_r[is.na(bcell_bcr$mu_freq_seq_r)] <- 0

# we only access the heavy chain SHM ratio 
SHM_info %>% group_by(barcode) %>% summarise(max=max(mu_freq_seq_r)) %>% dim()
SHM_cell <- bcell_filter@meta.data %>% left_join(SHM_info, by=c('barcode'='barcode')) %>%
    drop_na(mu_freq_seq_r)
ggboxplot(SHM_cell, "subtype", "mu_freq_seq_r", color="disease",
          palette=c("#4169E1", "#FF8C00")) + ggtitle("Total mutations(replace)") +
    xlab("Isotype") + ylab("Mutation frequency") +
    theme_bw() + stat_compare_means(aes(group=disease), label="p.signif") +
    facet_wrap(~c_call,scales="free")

SHM_cell %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    ggboxplot( "subtype", "mu_freq_seq_r", color="disease",
               palette=c("#4169E1", "#FF8C00")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() + stat_compare_means(aes(group=disease), label="p.signif") +
    facet_wrap(~c_call,ncol=4) + theme(axis.text.x=element_text(angle=30, hjust=1))

compare_group <-  list(c('HC','untreated'),c('treated','untreated'),c('HC','treated'))
SHM_cell$treatment <- factor(SHM_cell$treatment,levels=c('HC','untreated','treated')) 
# Fig.5c
SHM_cell %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHD')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", fill="treatment",
               palette=c("#B4D493", "#DA9494","#9FB1D4")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    facet_grid(subtype ~ c_call)  +
    stat_compare_means(mapping=aes(treatment),label="p.signif",hide.ns=F,
                       comparisons =compare_group,label.y =c(0.1,0.12,0.14)) +
    xlab('') + theme_cowplot()+ theme(axis.text.x=element_text(angle=30, hjust=1))

## overall
### Total mutations(replace)
SHM_info %<>% left_join(meta %>% rownames_to_column() %>% select('rowname','treatment'),
                       by = c('sample'='rowname'))
ggboxplot(SHM_info, "c_call", "mu_freq_seq_r", fill = "treatment",
          palette =  c("#B4D493", "#DA9494","#9FB1D4")) + ggtitle("Replace mutations") +
    xlab("Isotype") + ylab("Mutation frequency")  +
    theme_pubr() + stat_compare_means(aes(group = treatment ), label = "p.signif", method = 'anova')

ggboxplot(SHM_info, "c_call", "mu_freq_seq_s", fill = "treatment",
          palette =  c("#B4D493", "#DA9494","#9FB1D4")) + ggtitle("Silent mutations") +
    xlab("Isotype") + ylab("Mutation frequency")  +
    theme_pubr() + stat_compare_means(aes(group = treatment ), label = "p.signif", method = 'anova')


################################################################################
#
# B Cell subtype BCR overlap
#
################################################################################
all_b <- merge(x= bcell_filter, y=plasma_filter_bcr)
intersect(Cells(all_b),combined_bcr.df$barcode) %>% length()

all_b_bcr <- combineExpression(combined_bcr, all_b, 
                               cloneCall="gene+nt", group.by="sample", proportion=FALSE, 
                               cloneTypes=c(Single=1, Small=10, Medium=100, Large=1000))
combined_bcr.df.filter <- combined_bcr.df %>% filter(barcode %in% Cells(all_b_bcr))

# Fig.5e
Idents(all_b_bcr) <- 'treatment'
hc_bcr <- subset(all_b_bcr,idents='HC')
Idents(hc_bcr) <- 'subtype'
clonalOverlap2(hc_bcr, cloneCall="gene+nt", method="jaccard",title='Health Control',limit=c(0,0.045)) 

untreat_bcr <- subset(all_b_bcr,idents='untreated')
Idents(untreat_bcr) <- 'subtype'
clonalOverlap2(untreat_bcr, cloneCall="gene+nt", method="jaccard",title='SLE Untreated',limit=c(0,0.045))

treat_bcr <- subset(all_b_bcr,idents='treated')
Idents(treat_bcr) <- 'subtype'
clonalOverlap2(treat_bcr, cloneCall="gene+nt", method="jaccard",title='SLE Treated',limit=c(0,0.045))


################################################################################
#
# BCR preparation for VDJ tools
#
################################################################################
contig_list_bcr.filter <- contig_list_bcr[contig_list_bcr$barcode %in%  Cells(all_b),]

# group by group and chain
for (  i in contig_list_bcr.filter$treatment %>% unique() ){
    for ( bcr_chain in  contig_list_bcr.filter$chain %>% unique()){
        print( i);print(bcr_chain)
        tmp <- contig_list_bcr.filter %>% filter(treatment ==  i & chain == bcr_chain)
        print(dim(tmp))
        write_csv(tmp, file=paste0('vdj/publication/input/', i,"_",bcr_chain,'_bcr.csv'),col_names=T)
    }
}

# group only by group
for ( i in contig_list_bcr.filter$treatment %>% unique() ){
    print( i)
    tmp <- contig_list_bcr.filter %>% filter(treatment ==  i )
    print(dim(tmp))
    write_csv(tmp, file=paste0('vdj/publication/input/', i,'_bcr.csv'),col_names=T)
}

# run 'snakemake -j24 -p' in corresponding dirs

################################################################################
#
# finish
#
################################################################################
sessionInfo()
rm(hc_bcr,untreat_bcr,treat_bcr,all_b,all_b_bcr);gc()