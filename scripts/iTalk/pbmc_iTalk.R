
library(Seurat)
library(tidyverse)
library(iTALK)
library(data.table)

setwd('/data/sle')
source('./scripts/function_R/utils.R')
output_path <- 'output_file/seurat/all_pbmc/'
load('./output_file/seurat/all_pbmc/downsample_final_concat_pbmc_without_pSS.rdata')

#------------------------------ Build iTALK data -------------------------------
# use downsample data first
downsample_count_mat <- sparse_to_dense(dgc_matrix = pbmc_final_downsample@assays$RNA@counts,formart = 'DF')
iTalk_data <- transpose(downsample_count_mat)
colnames(iTalk_data) <- rownames(downsample_count_mat)
rownames(iTalk_data) <- colnames(downsample_count_mat)
iTalk_data <- cbind(data.frame(cell_type = pbmc_final_downsample$subtype), iTalk_data)
iTalk_data <- cbind(data.frame(compare_group = pbmc_final_downsample$disease), iTalk_data)
unique(iTalk_data$cell_type)
unique(iTalk_data$compare_group)


#-------------------------- iTALK analysis pipeline ----------------------------
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
# interaction type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(get_color(len =length(cell_types)), names=cell_types)

# par(mfrow=c(1,2))
iTalk_res <- NULL
for(comm_type in comm_list){
    res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
    # res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    # plot by ligand category
    ## 1.overall network plot
    # NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    ## 2.top 20 ligand-receptor pairs
    # LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
    # title(comm_type)
    iTalk_res <- rbind(iTalk_res, res_cat)
}

iTalk_res_filter <-iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:50,]
NetView(iTalk_res_filter,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(iTalk_res_filter[1:20,],datatype='mean count',cell_col=cell_col, 
       link.arr.lwd=iTalk_res_filter$cell_from_mean_exprs[1:20], 
       link.arr.width=iTalk_res_filter$cell_to_mean_exprs[1:20])


#-------------------------- Disease Vs Healthy ---------------------------------
## significant ligand-receptor pairs between compare groups

# randomly assign the compare group to each sample
data<-data %>% mutate(compare_group=sample(2,nrow(data),replace=TRUE))
# find DEGenes of regulatory T cells and NK cells between these 2 groups
deg_t<-DEG(data %>% filter(cell_type=='regulatory_t'),method='Wilcox',contrast=c(2,1))
deg_nk<-DEG(data %>% filter(cell_type=='cd56_nk'),method='Wilcox',contrast=c(2,1))
# find significant ligand-receptor pairs and do the plotting
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
    res_cat<-FindLR(deg_t,deg_nk,datatype='DEG',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
    #plot by ligand category
    if(nrow(res_cat)==0){
        next
    }else if(nrow(res_cat>=20)){
        LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
    }else{
        LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
    }
    NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    title(comm_type)
    res<-rbind(res,res_cat)
}
if(is.null(res)){
    print('No significant pairs found')
}else if(nrow(res)>=20){
    res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
    NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
}else{
    NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
    LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
}
# I just randomly assigned the compare group to samples which has no biological difference for showing how to use the package.
# So there should be no significant genes to be expected. 