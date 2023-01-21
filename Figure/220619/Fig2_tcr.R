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


################################################################################
#
# TCR repertoire analysis via scRepertoire
#
################################################################################

## STEP1: build the combined_tcr object
load('final/scRepertoire/TCR/combined_tcr.rdata')
load('final/seurat/t_cell/05-tcell_merge_final.rdata')
meta <- read.csv('data/meta.csv', row.names=1)

## STEP2: filter the TCR repertoire
combined_tcr.df <- rbindlist(combined_tcr)
index <- paste0(t_cell_merge$orig.ident, '_', t_cell_merge$group, '_', 
                (str_split_fixed(Cells(t_cell_merge),pattern = '_',n=3))[,3])
combined_tcr.df.filter <- combined_tcr.df %>% filter(barcode %in% index)
combined_tcr.filter <- split(combined_tcr.df.filter, f=combined_tcr.df.filter$sample)


## STEP3: TCR repertoire diversity(Fig.2e left)
# clonalHomeostasis(combined_tcr.filter, cloneCall = "gene+nt",
#                   cloneTypes = c(Rare = 0.01, Small = 0.05,
#                                  Large = 0.1, Hyperexpanded = 1))
tmp <- clonalHomeostasis(combined_tcr.filter, cloneCall = "gene+nt",
                         cloneTypes = c(Rare = 0.01, Small = 0.05,
                                        Large = 0.1, Hyperexpanded = 1),exportTable = T) %>% 
  melt () %>% left_join(meta %>% rownames_to_column(),by = c('Var1'='rowname' )) %>%
  arrange(treatment,Var1)
# tmp$Var1 %<>% unique() %>% factor(.,levels = .)
tmp$Var1 <- fct_inorder(tmp$Var1)
tmp$sample <- paste0('Sample', c(1:22) %>% rep(4)  %>% sort())
tmp$sample <- fct_inorder(tmp$sample)
ggplot(tmp, aes(x = sample, y = value,  fill = Var2)) + 
  geom_bar(stat = "identity", position = "fill", color = "black", lwd = 0.25) + 
  scale_fill_manual(name = "Clonotype Group", values = c('#7BC6FF','#C7FDED','#FFB2AD','#E7525A') ) + 
  xlab("") +  ylab("Relative Abundance") + theme_classic() + 
  theme(axis.text.x=element_text(angle=30, hjust=1))

quantContig(combined_tcr.filter, cloneCall="gene", group = "ID", scale = TRUE) + 
  stat_compare_means()
################################################################################
#
# prepare VDJtools TCR input
#
################################################################################
# please use  tcr_vdj/publication/scripts/prepare_group_data.R to generate
# input contig_list_tcr_filter is TCR filter by barcode
dim(contig_list_tcr_filter)
contig_list_tcr_filter$treatment %>% table()

for ( group in contig_list_tcr_filter$treatment %>% unique() ){
  tmp <- contig_list_tcr_filter %>% filter(treatment == group )
  print(dim(tmp))
  write_csv(tmp, file = paste0('./tcr_vdj/publication/input/',group,'_tcr.csv'),col_names = T)
}


################################################################################
#
# TCR repertoire expansion
#
################################################################################
# input combined_tcr.filter is TCR filter by barcode

clonalHomeostasis(combined_tcr.filter, cloneCall = "gene+nt",
                  cloneTypes = c(Rare = 0.01, Small = 0.05,
                                 Large = 0.1, Hyperexpanded = 1))
tmp <- clonalHomeostasis(combined_tcr.filter, cloneCall = "gene+nt",
                         cloneTypes = c(Rare = 0.01, Small = 0.05,
                                        Large = 0.1, Hyperexpanded = 1),exportTable = T) %>% 
  melt () %>% left_join(meta %>% rownames_to_column(),by = c('Var1'='rowname' )) %>%
  arrange(treatment,Var1)
# tmp$Var1 %<>% unique() %>% factor(.,levels = .)
tmp$Var1 <- fct_inorder(tmp$Var1)
tmp$sample <- paste0('Sample', c(1:22) %>% rep(4)  %>% sort())
tmp$sample <- fct_inorder(tmp$sample)
ggplot(tmp, aes(x = sample, y = value,  fill = Var2)) + 
  geom_bar(stat = "identity", position = "fill", color = "black", lwd = 0.25) + 
  scale_fill_manual(name = "Clonotype Group", values = c('#7BC6FF','#C7FDED','#FFB2AD','#E7525A') ) + 
  xlab("") +  ylab("Relative Abundance") + theme_classic() + 
  theme(axis.text.x=element_text(angle=30, hjust=1))

# clonalProportion(combined_tcr.filter, cloneCall = "gene+nt")
# abundanceContig(combined_tcr.filter, cloneCall = "gene+nt", group = "group", scale = FALSE)
# 
# quantContig(combined_tcr.filter, cloneCall="gene+nt", scale = TRUE, group.by = 'group')

