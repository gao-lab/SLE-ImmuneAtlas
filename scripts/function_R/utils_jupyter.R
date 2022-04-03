# -------------------------------- Initialize ----------------------------------
setwd('/rd2/user/xiacr/sle')
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

library('tidyverse')

# -------------------------------- Functions -----------------------------------
# Convert seurat object to iTALK data
# @seu_obj       : seurat object
# @cell_type     : one column of seurat meta, identify cell type annotation
# @compare_group : one column of seurat meta, identify compare group
# @top_genes     : parametter in rawParse(), set 1-100
# >>> iTalk_data and highly exp genes for iTALK analysis
seu_to_iTALK <- function(seu_obj,cell_type = 'subtype',compare_group = 'disease',
                         top_genes =50 ){
  downsample_count_mat <- sparse_to_dense(dgc_matrix = seu_obj@assays$RNA@counts,formart = 'DF')
  iTalk_data <- transpose(downsample_count_mat)
  colnames(iTalk_data) <- rownames(downsample_count_mat)
  rownames(iTalk_data) <- colnames(downsample_count_mat)
  rm(downsample_count_mat)
  iTalk_data <- cbind(data.frame(cell_type = do.call(`$`,args = list(seu_obj,cell_type))), iTalk_data)
  iTalk_data <- cbind(data.frame(compare_group = do.call(`$`,args = list(seu_obj,compare_group))), iTalk_data)
  
  return(iTalk_data)
  # highly_exprs_genes <- rawParse(iTalk_data, top_genes=top_genes, stats="mean")
  # return(list(iTalk_data, highly_exprs_genes))
}


# Convert sparse matrix to data.table or matrix(solve as.Matrix() not work)
# @dgc_matrix    : sparse matrix in R
# @formart       : which format should function return(data.table, data.frame or matrix)
# @file_name     : tmp file name  
# @file_path     : tmp file path 
# @keep_row_name : if keep colname (gene names in Seurat)
# @threads       : threads to accelerate the fread and fwrite
# >>> default: data.table but can be matrix or data.frame
sparse_to_dense <- function(dgc_matrix, formart = 'data.table', file_name = NULL,
                            file_path = '/rd2/user/xiacr/sle/tmp/', keep_row_name = T, threads = 12){
  # library(Matrix)
  library(data.table)
  library(stringr)
  print(dim(dgc_matrix))
  n <- dim(dgc_matrix)[1]%/%1000
  m <- dim(dgc_matrix)[1]%%1000
  print(n)
  
  if(is.null(file_name)){
    file_name <- str_replace_all(paste0(date(),'.txt'), ' ','_')
    file_name <- paste0(file_path,file_name)
  }
  
  for(i in (1:n)){
    up_limit <- 1000*i
    down_limit <- up_limit - 999
    print('-------------------')
    print(down_limit)
    print("to")
    print(up_limit)
    tmp_mat <- as.matrix(dgc_matrix[down_limit:up_limit,])
    fwrite(tmp_mat, file_name, row.names = F, col.names = F, nThread = threads,
           append = T)
  }
  down_limit <- 1000*n+1
  up_limit <- dim(dgc_matrix)[1]
  tmp_mat <- as.matrix(dgc_matrix[down_limit:up_limit,])
  fwrite(tmp_mat, file_name, row.names = F, col.names = F, nThread = threads,
         append = T)
  
  expr <- fread(file_name, header = F, nThread = threads)
  
  if(formart %in% c('matrix','Matrix','mat','M')){
    print('return matrix')
    expr <- as.matrix(expr)
    colnames(expr) <- colnames(dgc_matrix)
    if(keep_row_name){
      rownames(expr) <- rownames(dgc_matrix)
    }
  }else if(formart %in% c('data.frame','Data.Frame','df','DataFrame','DF')){
    print('return data.frame')
    expr <- as.data.frame(expr)
    colnames(expr) <- colnames(dgc_matrix)
    rownames(expr) <- rownames(dgc_matrix)
  }else{
    print('return data.table')
    colnames(expr) <- colnames(dgc_matrix)
    if(keep_row_name){
      expr <- cbind(rownames(dgc_matrix),expr)
    }
  }
  return(expr)
}


# Quick get color via RColorBrewer
# ! If len is Null print the RColorBrewer
# @len     : length of the colors you want
# @set     : which set of RColorBrewer do you want 
# @set_len : length of the set you want
# >>> color in 256 rgb format
get_color <-function(len = NULL,set = 'Set1',set_len = 9){
  if(is.null(len)){
    print(display.brewer.all())
    return(NULL)}
  library(RColorBrewer)
  getPalette = colorRampPalette(brewer.pal(set_len, set))
  return(getPalette(len))
}

# Run GSEA on seurat data set
# @seu_obj    : a seurat obejct which scaled and clustered
# @group_by   : which meta.data should be used
# @focus      : which sub class should we focus on @group_by 
# @title      : change the plot title
# @pvalue_cut : significance p value cut off (0.05)
# @plot_lengrh: how many top gsea results should be plot(20)
# @species    : species('Homo sapiens')
# @category   : which MisgDB database to use('C2')
# ! adapted from https://erilu.github.io/single-cell-rnaseq-analysis/ and
# ! https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
# @seu_obj  :seurat object  
# >>> None (plot in rstudio)
plot_gsea <- function(seu_obj,group_by,focus,title = 'please change title',
                    pvalue_cut= 0.05,plot_length =20,
                    species ="Homo sapiens",category = 'C2'){
  library(msigdbr)
  library(fgsea)
  library(presto)
  seu_obj.genes <- wilcoxauc(seu_obj,group_by) # fast rank for seurat object 
  m_df<- msigdbr(species = species, category = category)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  # seu_obj.genes %>%
  #   dplyr::filter(group == focus) %>%
  #   arrange(desc(logFC), desc(auc)) %>%
  #   head(n = 10)
  focus.genes<- seu_obj.genes %>%
    dplyr::filter(group == focus) %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  ranks<- deframe(focus.genes)
  fgsea_results<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  
  fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()
  
  plot <- fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = title)+
    theme(axis.text.y = element_text(size = 7), 
          plot.title = element_text(hjust = 1)) %>% 
    return(plot)

}


# Do harmony
# ! can add the funtion that save harmony convergence plot to pdf
# @seu_obj   : a seurat object
# @harmony_slot:
#             use which slot to remove the batch effect, such as 'sample'
# @from_begin: if do seurat analysi from beginning
# @max.iter  : max iter loop number in Harmony
# @scale_all : use all of the genes to scale(recommend False to save memory 
#              and speed up)m, but necessaery when using DoHeatmp()
# @res       : resolution of the clustering ,can set a vector like c(0.5,0.8,1)
# @qc        : if calculate MT and ribosome gene proportion
# @feature_exclude 
#            : genes excluded from HVGs so that their will not influence the clustering 
#              and UMAPauto remove IG TR RP MT gene family, only work when using Gene Symbol 
# @reg       : if regress the variable in @vars.reg
# @vars.reg  : only work when @reg is True
# >>> seurat object clustered and with UMAP
do_harmony <- function(seu_obj = seu_obj,harmony_slot='orig.ident',from_begin=F,max.iter=20,
                       scale_all = F, res = 0.8, qc = F,
                      feature_exclude = '^HLA-|^IG[HJKL]|^RNA|^MT-|^RP[SL]|^TR[ABDY][VJ]',
                      reg = F,
                      vars.reg = c('nCount_RNA','pct_counts_mt')){
  
  library(Seurat)
  library(harmony)
  set.seed(520)
  
  if(from_begin){
    if(qc){
      seu_obj <- PercentageFeatureSet(seu_obj, "^MT-", col.name = "percent_mito")
      seu_obj <- PercentageFeatureSet(seu_obj, "^RP[SL]", col.name = "percent_ribo")
    }
    seu_obj <- NormalizeData(object = seu_obj)
    seu_obj <- FindVariableFeatures(object = seu_obj)
    seu_obj <- CellCycleScoring(seu_obj,s.features=cc.genes$s.genes,
                                g2m.features=cc.genes$g2m.genes,set.ident=T)
    hvg <- VariableFeatures(seu_obj)
    
    feature_exclude = feature_exclude # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
    hvg = grep(feature_exclude, hvg, invert=T, value=T)
    
    if(scale_all){
      all.genes <- rownames(seu_obj)
      seu_obj <- ScaleData(object = seu_obj, features = all.genes)
    }else if(reg){
      seu_obj <- ScaleData(object = seu_obj, vars.to.regress =vars.reg , features = hvg)
    }else{
      seu_obj <- ScaleData(object = seu_obj, features = hvg)
    }
    
    seu_obj <- RunPCA(object = seu_obj, verbose = F, features = hvg )
    # seu_obj <- FindNeighbors(object = seu_obj)
    # seu_obj <- FindClusters(object = seu_obj, resolution = res)
    # seu_obj <- RunUMAP(object = seu_obj,dims = 1:30)
  }
# pdf()
  seu_obj <- seu_obj %>% 
    RunHarmony(group.by.vars = harmony_slot, plot_convergence = F, max.iter.harmony = max.iter)%>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = res) %>% 
    identity()
  return(seu_obj)
}

# Call function background mode
# ! must know and send the all the parameters of the original function 
# ! should assign the para of @func clearly, such as seu_obj = t.cell
# ! not need assign the output variable of back_fun(), it will use 'out_name'
# @func      : the name of a current function
# @out_name  : name of the export variable
# @job_name  : name of the job show on rstudio surface
# @import    : import parameters for job()
# @export    : export parameters for job()
# @...       : the list of the current function parameters
back_run <- function(func, out_name, job_name, ...){
  job::job({
    print('----start----')
    # print(ls())
    source('./scripts/function_R/utils.R')
    assign(out_name, func(...))
    print('----finish----')
    print(ls())
  },import = c(func,out_name,...), title = job_name )
  
  print('job submitted')
}

# Trans Ensemble ID to Gene Symbol
# !Can add a para to choose mouse or human
# @mat:  matrix whose row name is Ensemble ID
# >>> matrix whose row name is Gene Symbol
trans_Ensemble_Symbol <- function(mat, species = 'human'){
  library(biomaRt)
  biomart_dataset = 'hsapiens_gene_ensembl'
  out_name = 'hgnc_symbol'
  if(species == 'mouse'){
    biomart_dataset = 'mmusculus_gene_ensembl'
    out_name = 'mgi_symbol'
  }
  mart <- useDataset(biomart_dataset, useMart("ensembl"))
  trans_list <- getBM(filters= "ensembl_gene_id", 
                      attributes= c("ensembl_gene_id",out_name),
                      values=rownames(mat),mart= mart)
  trans_list <- trans_list[!duplicated(trans_list[,1]),]
  trans_list <- trans_list[!duplicated(trans_list[,2]),]
  trans_list <- trans_list[-which(trans_list[,2]==''),]
  gene_list <- merge(as.data.frame(rownames(mat)),trans_list, 
                     by.x= 'rownames(mat)', by.y='ensembl_gene_id',
                     all.x = T,sort = F)
  mat <- mat[!is.na(gene_list[,2]),]
  rownames(mat) <- gene_list[,2][!is.na(gene_list[,2])]
  return(mat)
}

# Beautify the seurat function DotPlot()
# @para seu_obj:     seurat object
# @para marker_list: marker genes list
DotPlot2 <- function(seu_obj,marker_list, group.by = NULL,split.by=NULL){
  print(
    DotPlot(seu_obj, features = marker_list, group.by = group.by, split.by =split.by ) + 
      scale_colour_gradient2(low = "#4169E1", mid = "#F5F5F5", high = "#DC143C")+
      scale_size(breaks = c(0, 25, 50, 75)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('')
  )
}

# Beautify the seurat function DoHeatmap()
# @para seu_obj:     seurat object
# @para marker_list: marker genes list
DoHeatmap2 <- function(seu_obj,marker_list){
  if(length(Cells(seu_obj)) > 15000){
    print('please downsample')
  }else{
    print(
      DoHeatmap(seu_obj, marker_list)+guides(color = FALSE) +
        scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), 
                             mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                             midpoint = 0, guide = "colourbar", aesthetics = "fill")
    )
  }
}

# Compute the TPM value from the raw count matrix
# @para counts: raw counts(RNA:gene x cell)
# @para len: length of the gene(transcript), must match with @count gene
# >>> TPM matrix(gene x cell)
counts_to_tpm  <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# Do Seurat pipeline in background mode
# @name: name of the seurat object returned. Must specify for every job!
# Other: parameters are same as the function do_seurat()
# >>> seurat object clustered and with UMAP
do_seurat_back <- function(seu_obj, name ,scale_all = F, res = 0.8, qc = F,
                           feature_exclude = '^HLA-|^IG[HJKL]|^RNA|^MT-|^RP[SL]|^TR[ABDY][VJ]',
                           reg = F,vars.reg = c('nCount_RNA','pct_counts_mt')){
  print('start')
  # rm(seu_obj)
  job::job({
    source('/rd2/user/xiacr/TLB/scripts/function_R/utils.R')
    # name <- do_seurat(seu_obj,scale_all=scale_all,res=res,qc=qc,
    # feature_exclude=feature_exclude,reg =reg,vars.reg=vars.reg)
    assign(name , do_seurat(seu_obj,scale_all,res,qc,feature_exclude,reg,vars.reg))
    # job::export()
  },import ='auto', title = paste0(name,'_job')
  )
  
}

# Do seurat pipeline automatically
# @seu_obj  : a seurat object
# @scale_all: use all of the genes to scale(recommend False to save memory 
#             and speed up)m, but necessaery when using DoHeatmp()
# @res      : resolution of the clustering ,can set a vector like c(0.5,0.8,1)
# @qc       : if calculate MT and ribosome gene proportion
# @feature_exclude: genes excluded from HVGs so that their will not influence the clustering 
#             and UMAPauto remove IG TR RP MT gene family, only work when using Gene Symbol 
# @reg      :if regress the variable in @vars.reg
# @vars.reg
#       only work when @reg is True
# >>> seurat object clustered and with UMAP
do_seurat <- function(seu_obj,scale_all = F, res = 0.8, qc = F,
                      feature_exclude = '^HLA-|^IG[HJKL]|^RNA|^MT-|^RP[SL]|^TR[ABDY][VJ]',
                      reg = F,
                      vars.reg = c('nCount_RNA','pct_counts_mt')){
  
  library(Seurat)
  library(Matrix)
  set.seed(520)
  if(qc){
    seu_obj <- PercentageFeatureSet(seu_obj, "^MT-", col.name = "percent_mito")
    seu_obj <- PercentageFeatureSet(seu_obj, "^RP[SL]", col.name = "percent_ribo")
  }
  seu_obj <- NormalizeData(object = seu_obj)
  seu_obj <- FindVariableFeatures(object = seu_obj)
  seu_obj <- CellCycleScoring(seu_obj,s.features=cc.genes$s.genes,
                              g2m.features=cc.genes$g2m.genes,set.ident=T)
  hvg <- VariableFeatures(seu_obj)
  
  feature_exclude = feature_exclude # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
  hvg = grep(feature_exclude, hvg, invert=T, value=T)
  
  if(scale_all ){
    all.genes <- rownames(seu_obj)
    seu_obj <- ScaleData(object = seu_obj, features = all.genes)
  } else if(reg ){
    seu_obj <- ScaleData(object = seu_obj, vars.to.regress =vars.reg )
  }else{
    seu_obj <- ScaleData(object = seu_obj)
  }
  
  seu_obj <- RunPCA(object = seu_obj, verbose = F, features = hvg )
  seu_obj <- FindNeighbors(object = seu_obj)
  seu_obj <- FindClusters(object = seu_obj, resolution = res)
  seu_obj <- RunUMAP(object = seu_obj,dims = 1:30)
  return(seu_obj)
}


















