

macro <- subset(mono_dc_filter,idents = 'Macrophage')
VlnPlot(macro , features = c('STAT1 ','CD86','CD163','STAT3','TLR7'),group.by = 'treatment',stack = T)
marker_macro_sle <- FindMarkers(macro, ident.1 = 'untreated', ident.2 = 'HC', group.by = 'treatment')
# treatment before and after do not have significant difference
marker_macro_sle.treat <- FindMarkers(macro, ident.1 = 'untreated', ident.2 = 'treated', group.by = 'treatment') %>%
     filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
marker_macro_sle %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
