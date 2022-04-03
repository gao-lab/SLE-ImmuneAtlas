# install.packages("./vdj/Platypus_3.3.2.tar.gz", repos = NULL, type="source")

library(Platypus)
library(Seurat)
library(tidyverse)
library(utils)
library(scales)
library(pals)
library(ggrepel)
library(igraph)

VDJ.out.directory.list_B <- list()
VDJ.out.directory.list_B[[1]] <- c("/data/teamwork/sle/data/xiehe/batch6_7_27/raw_data/vdj_wyf_bcr/outs/") # WYF
VDJ.out.directory.list_B[[2]] <- c("/data/teamwork/sle/data/xiehe/batch6_7_27/raw_data/vdj_zs_bcr/outs/") # ZS

GEX.out.directory.list <- list()
GEX.out.directory.list[[1]] <- c("/data/teamwork/sle/data/xiehe/batch6_7_27/raw_data/count_wyf_count/outs/") # WYF 
GEX.out.directory.list[[2]] <- c("/data/teamwork/sle/data/xiehe/batch6_7_27/raw_data/count_zs_count/outs/") # ZS

vgm_b <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list_B,                         
                        GEX.out.directory.list = GEX.out.directory.list,
                        exclude.on.cell.state.markers = c("CD3G+", "CD3E+", "CD4+",
                                                          "CD8A+","CD8B1+","C1QA+",
                                                          "C1QB+","RORA+","NKG7+"),
                        #Strict filtering to achieve best possible clustering of B cells
                        group.id = c("WYF","ZS"), trim.and.align = T,
                        parallel.processing  ='mclapply') 
saveRDS(vgm_b,"EAE_vgm_b.rds")


# ------------------- test ---------------------------------
EAE_Bcells <- PlatypusDB_fetch(
    PlatypusDB.links = c("kreiner2021b/ALL/ALL"),
    save.to.disk = F, 
    load.to.enviroment = F, 
    load.to.list = T, 
    combine.objects = T)
