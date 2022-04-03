# by chenrui
## input  1) vector of target cell type 2) vector of background
##        3) x lab name 
## output 1) boxplot with p value 

comprare_cell_ratio <- function(target_cell, background, x_lab_name){
 
  library(ggpubr) 
  library(RColorBrewer)
  library(ggplot2)
  
  theme_boxplot <- theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
                         axis.line=element_line(colour="black",size=0.25),
                         axis.title=element_text(size=13,face="plain",color="black"),
                         axis.text = element_text(size=12,face="plain",color="black"),
                         legend.position="none")
  
  ## here we focus on TLB ratio between 2 groups 
  
  tlb_pbmc_ratio <- table(target_cell)/table(child_all$orig.ident)*100
  tlb_b_ratio <- table(target_cell)/table(child_b$orig.ident)*100
  
  tlb_ratio_df <- data.frame(group = rep('unknown', length(tlb_pbmc_ratio)), 
                             ratio = tlb_pbmc_ratio,ratio2 = tlb_b_ratio, row.names = 2)
  tlb_ratio_df$group[which(rownames(tlb_ratio_df) %in% cHD_sample_list)] <- x_lab_name[1]
  tlb_ratio_df$group[which(rownames(tlb_ratio_df) %in% cSLE_sample_list)] <- x_lab_name[2]
  table(tlb_ratio_df$group)
  
  
  compaired <- list(x_lab_name)
  col_list <- c(brewer.pal(7, "Set2")[c(1, 2, 4)])
  # wilcox.test
  ggboxplot(tlb_ratio_df, x = "group", y = "ratio.Freq",
            fill = "group", palette = col_list,
            add = "jitter", size=0.5) +
    stat_compare_means(comparisons = compaired, method = "wilcox.test") + 
    theme_classic() + theme_boxplot
  ggboxplot(tlb_ratio_df, x = "group", y = "ratio2.Freq",
            fill = "group", palette = col_list,
            add = "jitter", size=0.5) +
    stat_compare_means(comparisons = compaired, method = "wilcox.test") + 
    theme_classic() + theme_boxplot

}
