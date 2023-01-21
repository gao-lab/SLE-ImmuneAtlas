library(Seurat)

################################################################################
#
# patching the pictures
#
################################################################################
Idents(cd4_filter) <- 'subtype'
Idents(cd8_filter) <- 'subtype'
p11 <- DimPlot(bcell_filter, cells.highlight = Cells(subset(bcell_filter, idents ='B.IFN-response'))) + NoLegend() + NoAxes() + ggtitle('B cell')
p12 <- DimPlot(cd4_filter,  cells.highlight =  Cells(subset(cd4_filter, idents ='T.CD4.IFN-response')))+ NoLegend()+ NoAxes()+ ggtitle('CD4 T cell')
p13 <- DimPlot(cd8_filter,  cells.highlight =  Cells(subset(cd8_filter, idents ='T.CD8.IFN-response' )))+ NoLegend()+ NoAxes()+ ggtitle('CD8 T cell')

p21 <- FeaturePlot(bcell_filter, features = c('IFI6'), min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p22 <- FeaturePlot(cd4_filter, features = c('IFI6'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p23 <- FeaturePlot(cd8_filter, features = c('IFI6'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()

p31 <- FeaturePlot(bcell_filter, features = c('MX1'), min.cutoff = 2, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p32 <- FeaturePlot(cd4_filter, features = c('MX1'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p33 <- FeaturePlot(cd8_filter, features = c('MX1'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()

p41 <- FeaturePlot(bcell_filter, features = c('ISG15'), min.cutoff = 2, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p42 <- FeaturePlot(cd4_filter, features = c('ISG15'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()
p43 <- FeaturePlot(cd8_filter, features = c('ISG15'),min.cutoff = 1, cols = c('grey','red'),pt.size = 0.6)+ NoLegend()+ NoAxes()

# merge
(p11 | p12 | p13 ) /(p21 |p22 |p23) /(p31 |p32 |p33) /(p41 |p42 |p43)

# DotPlot2(pbmc_all, marker_list  =c('IFITM3','IFI27') )
