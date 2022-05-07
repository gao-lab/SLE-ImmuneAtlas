library(ggplot2)
library(tidyverse)
library(ggpubr)

mip1_beta <- read.csv('/data/sle/data/blood_cytokine/MCP-1.csv')
g1 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL3', fill = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test') + xlab('') +
  ylab('CCL3 in blood')
g2 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL4', fill = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test') + xlab('') +
  ylab('CCL4 in blood')
g3 <- mip1_beta %>% filter(!group == 'other') %>% ggboxplot('group','CCL5', fill = c('#DA9494','#9FB1D4','#B4D493')) + 
  stat_compare_means(comparisons = list(c('before','hc'),c('after','hc')), method = 't.test') + xlab('') +
  ylab('CCL5 in blood')

g1 + g2 + g3
