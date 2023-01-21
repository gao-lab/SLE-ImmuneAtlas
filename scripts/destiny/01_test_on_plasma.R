library(destiny)
library(Seurat) 
library(SingleCellExperiment)

# README
# This test show the destiny runs slow in the default mode and not work when use 
# n_pcs parameters so I decide use DPT in scanpy instead

sce <- as.SingleCellExperiment(plasma_filter)
#this has the cell classification
table(sce$ident)

dm <- DiffusionMap(sce, verbose = TRUE,n_pcs = 50,n_threads =24)


cellLabels <- sce$ident
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Samples = cellLabels)
# pdf("./DC1_DC2.pdf", w=11, h=8.5)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Samples)) +
    geom_point()  + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
# dev.off()


sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # rank cells by their dpt dm1
sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # rank cells by their dpt dm2
sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # rank cells by their dpt dm1 reverse order
sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # rank cells by their dpt dm2 reverse order

SortedDM1 <- data.frame(DM1Sort = as.data.frame(colData(sce))$pseud_dm1,
                        Samples = as.data.frame(colData(sce))$ident)
SortedDM2 <- data.frame(DM2Sort = as.data.frame(colData(sce))$pseud_dm2,
                        Samples = as.data.frame(colData(sce))$ident)
SortedDM1R <- data.frame(DM1SortR = as.data.frame(colData(sce))$pseud_dm1R,
                         Samples = as.data.frame(colData(sce))$ident)
SortedDM2R <- data.frame(DM2SortR = as.data.frame(colData(sce))$pseud_dm2R,
                         Samples = as.data.frame(colData(sce))$ident)

ggplot(SortedDM1, aes(x=SortedDM1[,1], y=Samples,color=Samples)) +
    geom_jitter() + xlab("Diffusion component 1 (DC1)") + ylab("Samples") +
    ggtitle("Cells ordered by DC1")
ggplot(SortedDM2, aes(x=SortedDM2[,1], y=Samples,color=Samples)) +
    geom_jitter() + xlab("Diffusion component 2 (DC2)") + ylab("Samples") +
    ggtitle("Cells ordered by DC2")


library(plotly)

#interactive 2D
p2 = plot_ly(x=tmp$DC1, y=tmp$DC2, type="scatter", mode="markers", color=tmp$Samples, marker.size = 0.5)
htmlwidgets::saveWidget(as_widget(p2), "Interactive2D_DiffM.html", title = "Diffusion map")

#interactive 3D
p = plot_ly(x=tmp$DC1, y=tmp$DC2, z=tmp$DC3, type="scatter3d", mode="markers", color=tmp$Samples, marker = list(size = 2 ))
htmlwidgets::saveWidget(as_widget(p), "Interactive3D.html", title = "Diffusion map")

ggplot(SortedDM1R, aes(x=SortedDM1R[,1], y=Samples,color=Samples)) +
    geom_jitter() + xlab("Minus Diffusion component 1 (DC1)") + ylab("Samples") +
    ggtitle("Cells ordered by reversed DC1")
ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Samples,color=Samples)) +
    geom_jitter() + xlab("Minus Diffusion component 2 (DC2)") + ylab("Samples") +
    ggtitle("Cells ordered by reversed DC2")
