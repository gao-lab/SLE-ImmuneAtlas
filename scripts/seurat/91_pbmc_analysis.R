



load('./final/seurat/pbmc/03-pbmc_all_concat.rdata')
# modify the meta
pbmc_all$subtype[which(pbmc_all$subtype == 'Tfh')] <- 'T.CD4.Tfh'
pbmc_all$subtype[which(pbmc_all$subtype == 'T0_naive_ifn')] <- 'T.CD4.IFN-response'
pbmc_all$subtype[which(pbmc_all$subtype == 'Th1')] <- 'T.CD4.Th1'
pbmc_all$subtype[which(pbmc_all$subtype == 'Th1.cxcr3')] <- 'T.CD4.Th1.cxcr3'
pbmc_all$subtype[which(pbmc_all$subtype == 'Th17')] <- 'T.CD4.Th17'
pbmc_all$subtype[which(pbmc_all$subtype == 'Th2')] <- 'T.CD4.Th2'
pbmc_all$subtype[which(pbmc_all$subtype == 'Treg')] <- 'T.CD4.Treg'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD8.IFN-response')] <- 'T.CD8.IFN-response'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD8.mem')] <- 'T.CD8.mem'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD8.Naive')] <- 'T.CD8.Naive'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD8.Teff')] <- 'T.CD8.Teff'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD8.Tex')] <- 'T.CD8.Tex'
pbmc_all$subtype[which(pbmc_all$subtype == 'CD16.Mono')] <- 'Mono.CD16'


pbmc_all$main_type[which(pbmc_all$subtype == 'HPSC')] <- 'HPSC'

pbmc_all$main_type[which(pbmc_all$main_type == 'mDC')] <- 'cDC'

# pub
DimPlot(pbmc_all, group.by = 'main_type',label = T, cols = get_color(12), raster = F) + NoAxes()
DimPlot(pbmc_all, group.by = 'subtype',label = T, cols = get_color(41), raster = F, repel = T) + NoAxes()

# -------------------------------- Vis the marker ------------------------------

Idents(pbmc_all) <- 'subtype'
marker_pbmc <- c('CD3D','CCR7','IL7R','GZMA','NKG7','MKI67','KLRC1','CD79A','MZB1','LILRA4',
                 'ITGAX','CD86','CD14','FCGR3A','CD34','PPBP')
table(pbmc_all$main_type)
pbmc_all$main_type <- factor(pbmc_all$main_type, 
                               levels = rev(c("T.naive","T.cyto","T.prolife","NK",
                                              "Bcell","Plasma",'pDC',"cDC","Mono.CD14","Mono.CD16",'HPSC',
                                              'Platelet')))
# pub
p <- VlnPlot(pbmc_all, features = marker_pbmc, stack = T,group.by = 'main_type',cols = get_color(16)) + 
    NoLegend() + RotatedAxis() + NoGrid() + ylab('') + xlab('') 
p+  theme(axis.text.y=element_text(angle=300, hjust=1), axis.text.x=element_text(angle=270, hjust=1))

# pub
pbmc_all$subtype  %>% table()
pbmc_all$subtype <- factor(pbmc_all$subtype, 
                             levels = rev(c("T.CD4.naive","T.CD4.IFN-response","T.CD4.Tfh","T.CD4.Th1","T.CD4.Th1.cxcr3","T.CD4.Th17","T.CD4.Th2","T.CD4.Treg",
                                            "T.CD8.Naive","T.CD8.mem","T.CD8.Teff","T.CD8.IFN-response","T.CD8.Tex",
                                            "T.MAIT","T.yd",'T.prolife',"NK","NK.cd56+",
                                            "B.transition","B.naive","B.IFN-response","B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-",
                                            "plasmablast","plasma","plasma.IgG","plasma.IgA","pDC","cDC1","cDC2","Mono.CD14","Mono.CD14.APOBEC3A+",
                                            "Mono.CD14.LGALS2+","Mono.CD14.recruit","Mono.CD16","Macrophage","HPSC",
                                            "Platelet","Platelet.IFITM1+")))
marker_pbmc_detailed <- c('CD3D','CD4','IL7R','ICOS','TBX21','CXCR3','RORC','GATA3','FOXP3',
                          'CD8A','EOMES','GZMB','GZMA','HAVCR2','TRAV1-2','TRGV9','MKI67','KLRC1','NCAM1',
                          'CD79A','CD19','IGHD','IFIT1','CD27','MS4A1','TOP2A','MZB1','IGHG1','IGHA1','LILRA4','FCER1A','CD1C','CLEC9A','CD14',
                          'APOBEC3A','LGALS2','FCGR3A','C1QC','CD34',
                          'PPBP','IFITM1') %>% unique()
p <- VlnPlot(pbmc_all, features = marker_pbmc_detailed, stack = T,group.by = 'subtype',
             cols = get_color(length(marker_pbmc_detailed))) + 
    NoLegend() + RotatedAxis() + NoGrid() + ylab('') + xlab('') +
 theme(axis.text.y=element_text(angle=360, hjust=1), axis.text.x=element_text(angle=270, hjust=1))

DotPlot2(pbmc_all,marker_list = marker_pbmc_detailed, group.by = 'subtype')


pbmc_all$main_type[is.na(pbmc_all$main_type)] <- 'HPSC'

pbmc_all$subtype[which(pbmc_all$main_type == 'T.prolife')] <- 'T.prolife'
pbmc_all$subtype[is.na(pbmc_all$subtype) ] <- 'T.CD4.IFN-response'

# -------------------------------- Bar Plot  -----------------------------------
pbmc_all$orig.ident <- factor(pbmc_all$orig.ident, 
                           levels = c('QJY','ZH','ZMY1','ZS',
                                      'GW','GZR','HXR','HXX','LGY','LL','MXY','SQ','WYF','WYY','XH','ZPP',
                                      'HXR2','LL2','WYF2','XYY','XYY2','ZPP2'))
# pub
# Bar plot :main cell type distribution in PBMC
ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$orig.ident, 
                                      fill =pbmc_all$main_type))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pbmc_all$main_type))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values=get_color(14,set = 'Set1',set_len = 9)) +
     theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# pub
# Bar plot :subclusters distribution in PBMC
ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$orig.ident, 
                                          fill =pbmc_all$subtype))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pbmc_all$subtype))) + xlab('')+
    labs(fill="") +  scale_fill_manual(values=get_color(41,set = 'Set1',set_len = 11)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# -------------------------------- Save the data -------------------------------
save(pbmc_all, file = 'final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')
