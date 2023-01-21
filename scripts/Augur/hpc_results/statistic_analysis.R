library(tidyverse)

pbmc_subtype <- c(bcell_filter$subtype %>% unique(), t_cell_merge$subtype %>% unique(),
                    mono_dc_filter$subtype %>% unique(),
                  plasma_filter$subtype %>% unique(), platelet_filter$subtype %>% unique())
score <- c(0.648,0.652,0.605,0.693,0.593,0.642,0.672,
           0.647,0.706,0.708,0.635,0.740,0.625,0.670,0.672,
           0.637,0.693,0.603,0.707,0.570,0.718,0.698,
           0.709,0.682,0.556,
           0.547,0.673,0.632,0.664,0.674,0.601,0.709,0.626,0.590,0.565,
           0.784,0.730,0,0,0.573,0.749)

pbmc_augur_subtype_score <- data.frame(subtype = pbmc_subtype,score = score) %>% arrange(desc(score))

# c('naive','memory','effector','depress','other')
cell_group <- c('effector','other','effector','effector','effector','effector','other','effector', # 8
                'effector','effector','effector','memory','effector','effector','other','other', #16
                'memory','effector','effector','other','naive','memory','naive','memory','memory',
                'depress','other','other','depress','naive','naive','other','naive','other','other','depress',
                'other','naive','other','unknown','unknown')
pbmc_augur_subtype_score$cell_group <- cell_group 
pbmc_augur_subtype_score %<>% filter(!cell_group %in% c('unknown','other'))
main_group <- c('B','T','B','T','T','T','T','T','T','B',
                'T','T','B','T','T','B','B','T','B','T',
                'T','T','B','T','B','T','T')
pbmc_augur_subtype_score$main_group <- main_group 
t_b_subtype_score <- pbmc_augur_subtype_score %>% filter(!cell_group %in% c('unknown','other','depress'))
# pbmc_augur_subtype_score$subtype


ggboxplot(t_b_subtype_score,x = 'cell_group', y='score') + ylim(c(0.5,0.8))

t_b_subtype_score$cell_group <- factor(t_b_subtype_score$cell_group,levels = c('naive','memory','effector'))
ggline(t_b_subtype_score, x = "cell_group", y = "score", color = "main_group",point.size = 2,
       add = c("mean_se", "jitter"), palette = c("#DA9494", "#9FB1D4",'red','blue')) +
    xlab('Cell Status') + ylab('Augur score')+ stat_compare_means(method = "anova",label.x = 1.8) + 
    theme(legend.position="right",legend.title = element_text()) + labs(col="Cell Type")
