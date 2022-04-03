library(cowplot)
library(harmony)
library(tidyverse)
# library(RISC)
library(Seurat)
library(future)

setwd('/data/sle')
output_path <- './output_file/seurat/b_cell/'
source('./scripts/function_R/utils.R')

load('./output_file/seurat/b_cell/bcell_filter_harmony.rdata')
all_cell_num_df <- as.data.frame(table(all_cell$orig.ident)) 
colnames(all_cell_num_df) <- c('orig.ident','all_cell')

plot_df <- bcell_filter@meta.data %>% filter(group !='pSS_pah') %>% 
    group_by(orig.ident,subtype) %>% summarise(celltype = n()) %>% 
    mutate(b_cell = sum(celltype)) %>% mutate(B_ratio = celltype/b_cell*100) %>%
    left_join(bcell_filter@meta.data[,c(1,12,13,14)]  %>%  distinct() )  %>% 
    filter(subtype=='B.IFN-response') %>% 
    filter(!pair == c('unpaired')) %>% left_join(all_cell_num_df) %>% 
    mutate(All_ratio=celltype/all_cell*1000) %>% filter(!orig.ident =='XYY' )
plot_df$treatment <- factor(plot_df$treatment,levels =c('untreated','treated'))

ggplot(data = plot_df,aes(x=treatment,y=All_ratio,group = pair ,color=pair ,shape=pair ))+
    geom_point()+geom_line()+
    xlab("")+ ylab("Ratio")+theme_bw() 
    # theme(panel.grid.major=element_line(colour=NA),
    #       panel.background = element_rect(fill = "transparent",colour = NA),
    #       plot.background = element_rect(fill = "transparent",colour = NA),
    #       panel.grid.minor = element_blank(),#以上theme中代码用于去除网格线且保留坐标轴边框
    #       text = element_text(family = "STXihei"),#设置中文字体的显示
    #       # legend.position = c(.075,.915),#更改图例的位置，放至图内部的左上角
    #       legend.box.background = element_rect(color="black"))#为图例田间边框线

ggplot(data = plot_df,aes(x=treatment,y=B_ratio,group = pair ,color=pair ,shape=pair ))+
    geom_point()+geom_line()+
    xlab("")+ylab("Ratio")+
    theme_bw() # remove background grey
