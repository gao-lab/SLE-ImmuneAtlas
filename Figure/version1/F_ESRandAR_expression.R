
load('./final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')
DotPlot(pbmc_all, features = c('ESR1','AR','OR5AR1','PGR'), group.by = 'subtype')+ 
    scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
    # scale_size(breaks = c(0, 25, 50, 75)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
    theme_bw() +theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))

Idents(pbmc_all) <- 'subtype'
# pDC
DotPlot(pbmc_all %>% subset(idents = 'pDC'), features = c('ESR1','AR'), group.by = 'treatment')+
    scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
    theme_bw() +theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))

# plasma
DotPlot(pbmc_all %>% subset(idents = c('plasmablast','plasma','plasma.IgG','plasma.IgA')), 
        features = c('ESR1','AR'), group.by = 'treatment')+
    scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
    theme_bw() +theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
    theme(axis.text.x = element_text(angle = 90))

p <- list()
for (celltype in  c('plasmablast','plasma','plasma.IgG','plasma.IgA')) {
    p[[celltype]] <- DotPlot(pbmc_all %>% subset(idents = celltype), 
                             features = c('ESR1','AR'), group.by = 'treatment')+
        scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
        theme_bw() +theme(panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
        theme(axis.text.x = element_text(angle = 90))
}
do.call(grid.arrange,p)
