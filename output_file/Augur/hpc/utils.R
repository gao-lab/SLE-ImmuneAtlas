# -------------------------------- Initialize ----------------------------------
#setwd('/data/sle')
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
#.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")
# .libPaths("/home/liny/R/x86_64-pc-linux-gnu-library/4.0")

packages <- c("Seurat", "tidyverse","cowplot")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

# source other scripts
#source('/data/sle/scripts/function_R/do_seurat.R')
#source('/data/sle/scripts/function_R/markers.R')


# ----------------------Define some markers of cell types-----------------------
# TLR_genes
TLR_genes <- paste0('TLR',c(1:10))

# IL genes
if(exists('t.cd4.filter')){
IL_genes <- c(grep(rownames(t.cd4.filter), pattern = '^IL[0-9]$',value = T) %>% sort(),
              grep(rownames(t.cd4.filter), pattern = '^IL[0-9][0-9]$',value = T) %>% sort())
}else(print('not exists t.cd4.filter, need initialise IL_genes'))

# CD8 markers 
t_cd8_effect_marker <- c('GZMK','KLRD1','GZMB','CCL5','IFNG ')
t_cd8_tex_marker <- c('PDCD1','CTLA4','LAYN','HAVCR2')
t_cd8_naive_marker <- c('LEF1','TCF7')
t_cd8_resident_marker <- c('ITGA1','CXCR6')
t_cd8_act_marker <- c('IFNG','CCL4','CCL3')
t_cd8_temra_marker <- c('PRF1','NKG7')
t_cd8_mem_marker <- c('EOMES')
t_cd8_all_marker <- c(t_cd8_effect_marker,t_cd8_tex_marker,t_cd8_naive_marker,
                      t_cd8_resident_marker,t_cd8_act_marker,t_cd8_temra_marker,
                      t_cd8_mem_marker) %>% unique()

# T cell markers
t_cd4_th1_marker <- c('CXCR3','TBX21','IFNG','RUNX2')
t_cd4_th2_marker <- c('GATA3','IL13','IL5','IL17RB','LTB4R','CCR8')
t_cd4_th17_marker <- c('RORC','CCR6','RUNX1','PALLD')
t_cd4_treg_marker <- c('FOXP3','CCL5')
t_cd4_tfh_marker <- c('CXCR5','BCL6','TOX','TOX2','BTLA','ICOS','CD40LG','ICA1')
t_cd4_all_marker <- c(t_cd4_th1_marker,t_cd4_th2_marker,t_cd4_th17_marker,
                      t_cd4_tfh_marker,t_cd4_treg_marker) %>% unique()

t_yd_marker <- c('TRDV1','TRDV2','TRGV9')
t_MAIT_marker <- c('TRAV1-2','TRAJ33','TRAJ20','TRAJ12')
t_invarient_NKT_marker <- c('TRAV10','TRAJ18','CD1D')
t_rare_marker <- c(t_yd_marker,t_MAIT_marker,t_invarient_NKT_marker)
t_sub_marker <- c('MK167','TOP2A', # prolife
                  'LAG3','PDCD1','CD38','TIGIT','HAVCR2','CXCL13',   # CD8 exhaustion
                  'FASLG','IFNG','CCL5',  # CD8 efector
                  'GZMA','GZMB','GZMK','GZMH','NKG7','GNLY','IL2', # Cytotoxicity
                  'CXCR5','BCL6','TOX','TOX2','BTLA','ICOS','CD40LG','ICA1',  # CD4 Tfh
                  'FOXP3','IL32','ILRA','CCR8','LAYN','IKZF2', # CD4 Treg
                  'CD69','ITGAE','CTLA4' # Resident
                  )
# B cell and plasma markers 
new_b_marker_list.1 <- c('IFIT1','IFI44L','HSPA6','ACP5','RHOB',
                         'ZNF331','NR4A2','AREG','IL4R','FCER2',
                         'TCL1A')
new_b_marker_list.2 <- c('ISG15','IFI27','CXCR3','TBX21','ITGAX','CR2','FGR',
                         'TFEC','FCRL2','FCRL5','FCRL3','CCR7','IL6','CXCR4',
                         'CD5','CD9','AICDA','TLR7','PRDM1')
other_b_marker <- c('CD40LG','CD80','CD86','FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCER1A',
                    'FCER2','CCR7','CXCR5','CXCR4','CR2','CD1A','CD1B','CD1C')
#b_marker_list <- Reduce( union, list(b_marker_list.traditional,new_b_marker_list.1,
#                                    new_b_marker_list.2, other_b_marker))

# Monocyte, Macrophage and DC markers
macro_main_marker <- c('CD163','C1QC')
macro_sub_marker <- c('BAG3','G0S2','AREG','THBS1','IL1B','EREG','RNASE1','SEPP1','FOLR2','MRC1', # tissue resident
                      'APOE','APOC1','C1QC','HLA-DRB5','RGS1','SGK1','FCGR2B','CXCL10','CXCL9',
                      'GBP5','GBP1','SPP1','TREM2','ALOX5AP','FN1','GPNMB','FTL'
)
mono_dc_main_marker <-c('CD1C','CLEC9A','CD163','C1QC','CSF1R','TPSAB1')

# Platelet markers


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
                            file_path = getwd(), keep_row_name = T, threads = 12){
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
do_harmony <- function(seu_obj,harmony_slot='orig.ident',from_begin=F,max.iter=20,
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
    source('/data/TLB/scripts/function_R/utils.R')
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

# Detach All current loading R packages
detachAllPackages <- function() {
  basic.packages.blank <- c(    
    "stats",    
    "graphics",    
    "grDevices",    
    "utils",   
    "datasets",  
    "methods",    
    "base"    
  )    
  basic.packages <- paste("package:", basic.packages.blank, sep = "")   
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]   
  package.list <- setdiff(package.list, basic.packages)   
  if (length(package.list) > 0) {   
    for (package in package.list) {   
      detach(package, character.only = TRUE)   
    }   
  }    
}

# Plot sankey plot
# @reference : a colunm which is the source of the sankey
# @cluster   : a colunm which is the target of the sankey
# >>> sankey plot
getSankey <- function(reference, clusters, plot_width = 400, plot_height = 600, colors = NULL) {
  library(googleVis)
  Var1 <- value <- NULL
  res.all <- NULL
  for (j in names(table(reference))) {
    res <- NULL
    for (i in names(table(clusters))) {
      tmp <- length(intersect(which(clusters == i), which(reference == j)))
      res <- c(res, tmp)
    }
    res.all <- rbind(res.all, res)
  }
  colnames(res.all) <- names(table(clusters))
  rownames(res.all) <- names(table(reference))
  
  if (ncol(res.all) > 1) {
    res.all <- res.all[order(as.numeric(table(reference)), decreasing = TRUE), 
                       order(as.numeric(table(clusters)), decreasing = TRUE), drop = FALSE]
  }
  
  res <- reshape2::melt(res.all)
  res <- res[res$value != 0, ]
  
  if (ncol(res.all) > 1) {
    maxs <- res %>% dplyr::group_by(Var1) %>% dplyr::summarise(max = max(value))
    
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    maxs <- maxs[order(maxs$value, decreasing = TRUE), ]
    res <- res[res$value != res$max, ]
    res <- rbind(maxs, res)
    res <- res[, 1:3]
  }
  
  # remove cycles from the data
  res[, 1] <- paste0(res[, 1], " ")
  res[, 2] <- paste0(" ", res[, 2])
  
  colnames(res) <- c("From", "To", "# of cells")
  
  if (!is.null(colors)) {
    colors <- paste(colors, collapse = "', '")
    colors <- paste0("['", colors, "']")
  }
  
  Sankey <- gvisSankey(res, from = "From", to = "To", weight = "# of cells", 
                       options = list(width = plot_width, 
                       height = plot_height, sankey = paste0("{
                node:{
                    label:{
                        fontName:'Arial',
                        fontSize:11,color:
                        '#000000',
                        bold:true,
                        italic:false
                    },
                    colors:'#FFFFFF',
                    nodePadding:12
                },", 
        if (!is.null(colors)) {
        paste0("link:{
                    colorMode: 'source',
                    colors: ", 
                    colors, "
                },")
                              }, "iterations:0
            }")))
  
  return(Sankey)
}

# Plot sankey plot
# @reference : a colunm which is the source of the sankey
# @cluster   : a colunm which is the target of the sankey
# >>> sankey plot
getmy.sankey=function (reference, clusters,reorder_label=F,plot_width = 400, plot_height = 600,fontsize=12) 
{
  Var1 <- value <- NULL
  stopifnot(class(reference)=="factor")
  stopifnot(class(clusters)=="factor")
  res.all <- NULL
  for (j in names(table(reference))) {
    res <- NULL
    for (i in names(table(clusters))) {
      tmp <- length(intersect(which(clusters == i), which(reference == 
                                                            j)))
      res <- c(res, tmp)
    }
    res.all <- rbind(res.all, res)
  }
  colnames(res.all) <- names(table(clusters))
  rownames(res.all) <- names(table(reference))
  
  if (ncol(res.all) > 1) {
    res.all <- res.all[order(as.numeric(table(reference)), 
                             decreasing = TRUE), order(as.numeric(table(clusters)), 
                                                       decreasing = TRUE), drop = FALSE]
  }
  res <- reshape2::melt(res.all)
  res <- res[res$value != 0, ]
  if (ncol(res.all) > 1) {
    maxs <- res %>% dplyr::group_by(Var1) %>% dplyr::summarise(max = max(value))
    res <- merge(res, maxs)
    maxs <- res[res$value == res$max, ]
    maxs <- maxs[order(maxs$value, decreasing = TRUE), ]
    res <- res[res$value != res$max, ]
    if(reorder_label){
      maxs=maxs[match(rev(levels(reference)),maxs$Var1),]
    }
    res <- rbind(maxs, res)
    res <- res[, 1:3]
  }
  res[, 1] <- paste0(res[, 1], " ")
  res[, 2] <- paste0(" ", res[, 2])
  colnames(res) <- c("From", "To", "# of cells")
  colors_link <- colorRampPalette(brewer.pal(5, "Set2"))(length(unique(reference)))
  colors_link_array <- paste0("[", paste0("'", colors_link,"'", collapse = ','), "]")
  colors_node <- colorRampPalette(brewer.pal(5, "Dark2"))(length(unique(reference))+length(unique(clusters)))
  colors_node_array <- paste0("[", paste0("'", colors_node,"'", collapse = ','), "]")
  opts <- paste0("{
                 link: { colorMode: 'source',
                 colors: ", colors_link_array ," },
                 node: { colors: ", colors_node_array ," ,
                 label:{fontSize:",fontsize,",fontName:'Arial'},
                 fontName: 'Arial',
                 fontSize: ", fontsize, " },
                 iterations:0
}" )
  Sankey=gvisSankey(res, from = "From", to = "To", weight = "# of cells",options=list(width = plot_width, height = plot_height, sankey=opts))
  return(Sankey)
}

# Help the function seu_plot_heatmap() to choose the DEGs
get_top_genes_new <- function (dataset, markers, n) 
{
  library(Scillus)
  int_features <- rownames(dataset@assays$RNA@scale.data)
  df <- markers %>% filter(.data$gene %in% int_features) %>% 
    arrange(desc(.data$avg_log2FC)) %>% group_by(.data$cluster) %>% 
    filter(row_number() <= n) %>% arrange(.data$cluster)
  return(df$gene)
}
are_colors_new <- function (x) 
{
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)
  })
}
set_colors_new <- function (pal, n) 
{
  if (all(pal %in% rownames(brewer.pal.info))) {
    num <- c()
    for (i in seq(length(pal))) {
      num[i] <- brewer.pal.info[pal[i], ][[1]]
    }
    full_pal <- do.call(c, map2(.x = num, .y = pal, .f = brewer.pal))
  }
  else if (all(are_colors(pal))) {
    full_pal <- pal
  }
  else {
    stop("Incorrect palette setup. Please input valid RColorBrewer palette names or color names.")
  }
  if (n <= length(full_pal)) {
    return(full_pal[1:n])
  }
  else {
    warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.", 
            immediate. = TRUE)
    return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
  }
}


# A modify method for DoHeatmap which can plot heatmap with multi panel
# @ dataset     : a seurat object whose data has been scaled
# @ markers     : marker gene table of all variable(results of FindAllMarkers())
# @ sort_var    : which ident used to group(must inculed in para 'anno_var')
# @ n           : genes choose to plot for every variable
# @ anno_var    : which idents used to anno on the top of heatmap
# @ anno_colors : colors of anno var, must have same length
# @ hm_limit    : plot para
# @ hm_colors   : plot para
# @ row_font_size: plot para
seu_plot_heatmap <- function (dataset, markers, sort_var = c("seurat_clusters"), 
                              n = 8, anno_var, anno_colors, hm_limit = c(-2, 0, 2), hm_colors = c("#4575b4", 
                                                                                                  "white", "#d73027"), row_font_size = 12) 
{
  library(Scillus)
  library(ComplexHeatmap)
  library(circlize)
  mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), 
                      slot = "scale.data")
  if (is.data.frame(markers)) {
    genes <- get_top_genes_new(dataset, markers, n)
  }
  else if (is.character(markers)) {
    genes <- markers
  }
  else {
    stop("Incorrect input of markers")
  }
  mat <- mat[match(genes, rownames(mat)), ]
  anno <- dataset@meta.data %>% rownames_to_column(var = "barcode") %>% 
    arrange(!!!syms(sort_var))
  mat <- t(mat)
  mat <- mat[match(anno$barcode, rownames(mat)), ]
  mat <- t(mat)
  annos <- list()
  for (i in seq_along(1:length(anno_var))) {
    err_msg <- paste("Incorrect specification for annotation colors for", 
                     anno_var[i])
    value <- anno[[anno_var[i]]]
    if (is.numeric(value)) {
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 
                                                              "qual"])) {
        n <- brewer.pal.info[anno_colors[[i]], ]["maxcolors"][[1]]
        pal <- brewer.pal(n = n, name = anno_colors[[i]])
        col_fun <- colorRamp2(c(min(value), stats::median(value), 
                                max(value)), c(pal[2], pal[(n + 1)/2], pal[n - 
                                                                             1]))
      }
      else if (length(anno_colors[[i]]) == 3 & all(are_colors_new(anno_colors[[i]]))) {
        col_fun <- colorRamp2(c(min(value), stats::median(value), 
                                max(value)), anno_colors[[i]])
      }
      else {
        stop(err_msg)
      }
      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]], 
                              col = list(a = col_fun), border = TRUE, annotation_label = anno_var[i])
    }
    else {
      l <- levels(factor(anno[[anno_var[i]]]))
      if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
        col <- set_colors_new(anno_colors[[i]], length(l))
      }
      else if (length(anno_colors[[i]]) >= length(l) & 
               all(are_colors_new(anno_colors[[i]]))) {
        col <- anno_colors[[i]]
      }
      else {
        stop(err_msg)
      }
      names(col) <- l
      col <- col[!is.na(names(col))]
      col <- list(a = col)
      ha <- HeatmapAnnotation(a = anno[[anno_var[i]]], 
                              col = col, border = TRUE, annotation_label = anno_var[i])
    }
    names(ha) <- anno_var[i]
    annos[[i]] <- ha
  }
  annos <- do.call(c, annos)
  annos@gap <- rep(unit(1, "mm"), length(annos))
  ht <- Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, 
                heatmap_legend_param = list(direction = "horizontal", 
                                            legend_width = unit(6, "cm"), title = "Expression"), 
                col = colorRamp2(hm_limit, hm_colors), show_column_names = FALSE, 
                row_names_side = "left", row_names_gp = gpar(fontsize = row_font_size), 
                top_annotation = annos)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
}

# Assign a color to each cell based on some value
# @cell_vars Vector indicating the value of a variable associated with cells.
# @pal_fun   Palette function that returns a vector of hex colors, whose
#            argument is the length of such a vector.
# @ ...      Extra arguments for pal_fun.
# >>>return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

