

Idents(bcell_filter) <- 'treatment'
bcell_treat <- subset(bcell_filter, idents = c('HC'), invert = T)
Idents(bcell_treat) <- 'subtype'

marer_b_treat <- FindMarkers(bcell_treat,ident.1 = 'treated', ident.2 = 'untreated', only.pos = F)
marer_b_treat %<>% filter(p_val_adj <0.05) %>% arrange(desc(avg_log2FC))

bcell_treat$subtype <- factor(bcell_treat$subtype)
VlnPlot(bcell_treat,features =c('CXCR3','CXCR4','CXCR5'), split.by = 'treatment', split.plot = T ,group.by = 'subtype' )

DotPlot(bcell_treat, features = c('TLR7','TLR9'), split.by = 'subtype', cols = get_color(8))

subset(bcell_treat, idents = 'B.transition') %>% DotPlot2(marker_list = c('TLR9','ISG15','TLR7','IL10RA'), group.by = 'treatment')
