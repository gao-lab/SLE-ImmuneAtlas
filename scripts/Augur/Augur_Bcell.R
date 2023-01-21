library(Augur)
library(ggplot2)
library(ggpubr)

# ----------------------  HC and SLE none treatment ----------------------------
Idents(bcell_filter) <- 'treatment'
b_tmp1 <- subset(bcell_filter,idents = c('HC','untreated'))
# augur.B = calculate_auc(input = b_tmp1@assays$RNA@counts %>% as.matrix(),
#                       meta = b_tmp1@meta.data,min_cells = 20,
#                       label_col = 'treatment',
#                       cell_type_col = 'subtype', n_threads=32)
# VERY cpu consume
back_run(func = calculate_auc,out_name = 'augur',job_name ='augur' ,
         input = b_tmp1@assays$RNA@counts %>% as.matrix(),
         meta = b_tmp1@meta.data,min_cells = 20,
         label_col = 'treatment',show_progress = T,
         cell_type_col = 'subtype', n_threads=32)


plot_lollipop(augur)

plot_umap(augur = augur, sc = b_tmp1)
plot_umap(augur, pbmc_final_downsample, mode = 'rank')
plot_scatterplot()


# ------------------------  HC and SLE treatment ------------------------------
Idents(bcell_filter) <- 'treatment'
b_tmp2 <- subset(bcell_filter,idents = c('HC','treated'))
# augur_treat = calculate_auc(input = b_tmp2@assays$RNA@counts %>% as.matrix(),
#                         meta = b_tmp2@meta.data,min_cells = 20,
#                         label_col = 'treatment',
#                         cell_type_col = 'subtype', n_threads=32)
# VERY cpu consume
back_run(func = calculate_auc,out_name = 'augur_treat',job_name ='augur_treat' ,
         input = b_tmp2@assays$RNA@counts %>% as.matrix(),
         meta = b_tmp2@meta.data,min_cells = 20,show_progress = T,
         label_col = 'treatment',
         cell_type_col = 'subtype', n_threads=32)

plot_lollipop(augur_treat)

plot_umap(augur_treat, b_tmp2)
plot_scatterplot(augur , augur_treat)

# compare the treated and not treated result
augur_df <- data.frame(augur$AUC) %>%  mutate(group = 'HC VS untreated') %>% 
    rbind(augur_treat$AUC %>%  mutate(group = 'HC VS treated'))
augur_df$group <- factor(augur_df$group,levels = c('HC VS untreated','HC VS treated'))
augur_df$auc <- round(augur_df$auc, digits = 2)
ggbarplot(augur_df, "cell_type", "auc",
          fill = "group", color = "group", palette = "Paired",
          label = TRUE, digits =2,
          position = position_dodge(0.9))

# ------------------------  sample level  compare ------------------------------
pairs <- names(bcell_filter$pair %>% table())[-3]
Idents(bcell_filter) <- 'group'
bcell_hc <- subset(bcell_filter, idents = 'HC')
Idents(bcell_filter) <- 'pair'
augur_list <- lapply(pairs,FUN = function(pair){
                        tmp <- subset(bcell_filter, idents = pair)
                        Idents(tmp) <- 'treatment'
                        before <- subset(tmp, idents = 'untreated')
                        after <- subset(tmp ,idents = 'treated')
                        before <- merge(before, bcell_hc)
                        after <- merge(after, bcell_hc)
                        before <- calculate_auc(input = before@assays$RNA@counts %>% as.matrix(),
                                                meta = before@meta.data,min_cells = 20,
                                                label_col = 'treatment',
                                                cell_type_col = 'subtype', n_threads=4)
                        after <- calculate_auc(input = after@assays$RNA@counts %>% as.matrix(),
                                                meta = after@meta.data,min_cells = 20,
                                                label_col = 'treatment',
                                                cell_type_col = 'subtype', n_threads=4)
                        augur_df <- data.frame(before$AUC) %>%  mutate(group = 'HC VS untreated') %>% 
                            rbind(after$AUC %>%  mutate(group = 'HC VS treated'))
                        augur_df$group <- factor(augur_df$group,levels = c('HC VS untreated','HC VS treated'))
                        return(augur_df)
                        
})



# -------------------  HC and SLE  untreated(down sample) -----------------------


# -------------------  HC and SLE  treated(down sample) -----------------------


# ------------------------------- Save Files -----------------------------------
save(augur,file = './output_file/Augur/')
save(augur,file = './output_file/Augur/')

# --------------------------- Run in HPC(abandon) ------------------------------
# 
# augur_sub = calculate_auc(input = expr_mat,
#                           meta = pbmc_final@meta.data,
#                           label_col = 'disease',
#                           cell_type_col = 'subtype', n_threads=24)
# save(augur_sub, file = './augur_sub.rdata')
# 
# augur_main = calculate_auc(input = expr_mat,
#                            meta = pbmc_final@meta.data,
#                            label_col = 'disease',
#                            cell_type_col = 'compare_meta', n_threads=24)
# save(augur_main, file = './augur_main.rdata')
# 
# pdf(file = './Augur_pics.pdf')
# print(plot_lollipop(augur_sub))
# print(plot_lollipop(augur_main))
# print(plot_umap(augur_sub, pbmc_final))
# print(plot_umap(augur_sub, pbmc_final, mode = 'rank'))
# print(plot_umap(augur_main, pbmc_final))
# print(plot_umap(augur_main, pbmc_final, mode = 'rank'))
# dev.off()
# 
# 
