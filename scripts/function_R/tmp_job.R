
#------------------------------ Install RISC -----------------------------------
library(devtools)
install.packages(c("Matrix", "matrixStats", "irlba", "doParallel", "foreach", 
                   "data.table", "Rtsne", "umap", "MASS", "pbmcapply", "Rcpp", 
                   "RcppEigen", "densityClust", "FNN", "igraph", "RColorBrewer", 
                   "ggplot2", "gridExtra", "pheatmap"))
install_github("https://github.com/bioinfoDZ/RISC.git")
library(RISC)
# install.packages("/Path/to/RISC_1.0.tar.gz", repos = NULL, type = "source")


#------------------------------ -----------------------------------