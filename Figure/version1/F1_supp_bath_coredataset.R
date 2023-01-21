

DimPlot(pbmc_all, cols = c('#386CB0','#BF5B17','#BEAED4'), group.by = 'treatment', raster = F, pt.size = 0.1)

# change the sample name
sample_name <- paste0(rep('sample',22),c(1:22))


# plot the sample level batch 
DimPlot(pbmc_all, cols = rev(get_color(22)), group.by = 'sample_name', raster = F, pt.size = 0.1)

# plot the sample cell type distribution
pbmc_all$orig.ident <- factor(pbmc_all$orig.ident, 
                              levels = c('QJY','ZH','ZMY1','ZS',
                                         'GW','GZR','HXX','LGY','MXY','SQ','WYY','HXR','LL','WYF','XH','ZPP',
                                         'HXR2','LL2','WYF2','XYY','XYY2','ZPP2'))
pbmc_all$sample_name <- 'unknown'
for (i in c(1: length(unique(pbmc_all$orig.ident)))) {
    tmp <-  c('QJY','ZH','ZMY1','ZS',
              'GW','GZR','HXX','LGY','MXY','SQ','WYY','HXR','LL','WYF','XH','ZPP',
              'HXR2','LL2','WYF2','XYY','XYY2','ZPP2')[i]
    pbmc_all$sample_name[which(pbmc_all$orig.ident == tmp)] <- paste0(pbmc_all$treatment[which(pbmc_all$orig.ident == tmp)],i)
}
table(pbmc_all$sample_name)
pbmc_all$sample_name <- factor(pbmc_all$sample_name, 
                              levels =names(table(pbmc_all$sample_name))[c(c(1:4), c(18:22),c(11:17), c(5:10))])

ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$orig.ident, 
                                      fill =pbmc_all$main_type))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pbmc_all$main_type))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values=get_color(14,set = 'Set1',set_len = 9)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$sample_name, 
                                      fill =pbmc_all$main_type))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pbmc_all$main_type))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values=get_color(14,set = 'Set1',set_len = 9)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
