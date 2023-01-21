# NOTE: need to run b_cell_SHM.R and plasma_SHM.R in advance 
#       so we have object 'SHM_cell' and 'SHM_plasma'

dim(SHM_plasma);dim(SHM_cell)
tmp <- intersect(colnames(SHM_cell),colnames(SHM_plasma))

SHM_all <- rbind(SHM_cell[,tmp], SHM_plasma[,tmp])
dim(SHM_all)


################################################################################
#
# Fig3 
#
################################################################################
compare_group <-  list(c('HC','untreated'),c('treated','untreated'),c('HC','treated'))
SHM_all$treatment <- factor(SHM_all$treatment,levels = c('HC','untreated','treated')) 
SHM_all$subtype[grepl(SHM_all$subtype,pattern = 'plasma')] <- 'plasma'
SHM_all %>% filter(!subtype %in% c('B.naive','B.IFN-response','B.transition')) %>%
    filter(!c_call %in% c('IGHG4','IGHD')) %>% 
    ggboxplot( 'treatment', "mu_freq_seq_r", fill = "treatment",
               palette = c("#B4D493", "#DA9494","#9FB1D4")) + ggtitle("Total mutations(replace)") +
    xlab("B cell subtype") + ylab("Mutation frequency") +
    theme_bw() + 
    facet_grid(subtype ~ c_call)  +
    stat_compare_means(mapping = aes(treatment),label = "p.signif",hide.ns = F,
                       comparisons =compare_group,label.y =c(0.1,0.12,0.14)) +
    xlab('') + theme_cowplot()+ theme(axis.text.x=element_text(angle=30, hjust=1))

for (variable in vector) {
    
}