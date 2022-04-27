# -------------------------------- Initialize ----------------------------------
setwd('/data/sle')
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
.libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.0")
# .libPaths("/home/liny/R/x86_64-pc-linux-gnu-library/4.0")

packages <- c("Seurat", "tidyverse","harmony","cowplot")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

# source other scripts
source('/data/sle/scripts/function_R/do_seurat.R')
source('/data/sle/scripts/function_R/markers.R')


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

t_gd_marker <- c('TRDV1','TRDV2','TRGV9')
t_MAIT_marker <- c('TRAV1-2','TRAJ33','TRAJ12','KLRB1','IL18RAP')
t_invarient_NKT_marker <- c('TRAV10','TRAJ18','CD1D')
t_rare_marker <- c(t_gd_marker,t_MAIT_marker,t_invarient_NKT_marker)
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
b_marker_list <- Reduce( union, list(b_marker_list.traditional,new_b_marker_list.1,
                                     new_b_marker_list.2, other_b_marker))

# Monocyte, Macrophage and DC markers
macro_main_marker <- c('CD163','C1QC')
macro_sub_marker <- c('BAG3','G0S2','AREG','THBS1','IL1B','EREG','RNASE1','SEPP1','FOLR2','MRC1', # tissue resident
                      'APOE','APOC1','C1QC','HLA-DRB5','RGS1','SGK1','FCGR2B','CXCL10','CXCL9',
                      'GBP5','GBP1','SPP1','TREM2','ALOX5AP','FN1','GPNMB','FTL'
)
dc_marker <- c('IRF4','ITGAM')
mono_dc_main_marker <-c('CD1C','CLEC9A','CD163','C1QC','CSF1R','TPSAB1')

# Platelet markers

#Add scores (from https://github.com/scCOVID-19/COVIDPBMC/blob/69d1047eb8fd1f775acabb9d8453452853076dde/myeloid/myeloid_signature_scoring_publish.ipynb)
#https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html
TNF_genes = c("ABCA1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2",
              "BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44",
              "CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11",
              "CXCL2","CXCL3","CXCL6","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1",
              "EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2",
              "GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5",
              "IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1",
              "IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF",
              "MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2",
              "NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PER1","PFKFB3","PHLDA1",
              "PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1",
              "REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6",
              "SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2",
              "TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10",
              "TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")

#https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_IL6_JAK_STAT3_SIGNALING
JAK_genes = c("A2M","ACVR1B","ACVRL1","BAK1","CBL","CCL7","CCR1","CD14","CD36","CD38","CD44","CD9","CNTFR","CRLF2",
             "CSF1","CSF2","CSF2RA","CSF2RB","CSF3R","CXCL1","CXCL10","CXCL11","CXCL13","CXCL3","CXCL9","DNTT","EBI3",
             "FAS","GRB2","HAX1","HMOX1","IFNAR1","IFNGR1","IFNGR2","IL10RB","IL12RB1","IL13RA1","IL15RA","IL17RA",
             "IL17RB","IL18R1","IL1B","IL1R1","IL1R2","IL2RA","IL2RG","IL3RA","IL4R","IL6","IL6ST","IL7","IL9R","INHBE",
             "IRF1","IRF9","ITGA4","ITGB3","JUN","LEPR","LTB","LTBR","MAP3K8","MYD88","OSMR","PDGFC","PF4","PIK3R5","PIM1",
             "PLA2G2A","PTPN1","PTPN11","PTPN2","REG1A","SOCS1","SOCS3","STAM2","STAT1","STAT2","STAT3","TGFB1","TLR2","TNF",
             "TNFRSF12A","TNFRSF1A","TNFRSF1B","TNFRSF21","TYK2")

#https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_RESPONSE_TO_TYPE_I_INTERFERON
IFN_genes = c("ABCE1","ADAR","BST2","CACTIN","CDC37","CNOT7","DCST1","EGR1","FADD","GBP2","HLA-A","HLA-B","HLA-C",
             "HLA-E","HLA-F","HLA-G","HLA-H","HSP90AB1","IFI27","IFI35","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2",
             "IFITM3","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7",
             "IFNA8","IFNAR1","IFNAR2","IFNB1","IKBKE","IP6K2","IRAK1","IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8",
             "IRF9","ISG15","ISG20","JAK1","LSM14A","MAVS","METTL3","MIR21","MMP12","MUL1","MX1","MX2","MYD88","NLRC5",
             "OAS1","OAS2","OAS3","OASL","PSMB8","PTPN1","PTPN11","PTPN2","PTPN6","RNASEL","RSAD2","SAMHD1","SETD2","SHFL",
             "SHMT2","SP100","STAT1","STAT2","TBK1","TREX1","TRIM56","TRIM6","TTLL12","TYK2","UBE2K","USP18","WNT5A","XAF1",
             "YTHDF2","YTHDF3","ZBP1")



# -------------------------------- Functions -----------------------------------
# calculate the overlap before and after the treatment (just simplify the scatterClonotype())
# @ <all>    : please refer scatterClonotype()
# @ filter_NA: filter the cell without heavy chain
treat_overlap <- function (df, cloneCall = "gene+nt", x.axis = NULL, dot.size = 'total',
                           y.axis = NULL, chain = "both", split.by = NULL ,filter_NA = T){
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  axes <- which(names(df) %in% c(x.axis, y.axis, dot.size))
  if (chain != "both") {
    for (i in axes) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  x.df <- as.data.frame(table(df[[x.axis]][, cloneCall]))
  colnames(x.df)[2] <- x.axis
  y.df <- as.data.frame(table(df[[y.axis]][, cloneCall]))
  colnames(y.df)[2] <- y.axis
  combined.df <- merge(x.df, y.df, by = "Var1", all = TRUE)
  combined.df[is.na(combined.df)] <- 0
  combined.df <- combined.df[!grepl('NA_NA_',combined.df$Var1),]
  return(combined.df)
}
environment(treat_overlap) <- asNamespace('scRepertoire')


# Beautify the vizGenes()
# @combined_xcr : a seurat object which
# @gene         : see vizGenes()
# @chain        : see vizGenes()
# @title        : title of the plot
# @limit        : set limit of max value(value greater than this vill be white !!!)
# >>> Nothing (plotting)
vizGenes2 <- function(combined_xcr, gene = 'V', chain='IGH',y.axis = 'J',title= 'Need Title',limit = c(0,0.1)){
  vizGenes(combined_xcr, gene = gene, chain = chain, y.axis = y.axis , 
           plot = "heatmap", scale = TRUE, order = "gene") +  theme_void()+
    ggtitle(title) + theme(axis.text.x = element_text(angle = 90,vjust = 1,size = 10,hjust = 1)) + 
    theme(axis.text.x = element_text(size = 8),axis.text.y =element_text(size = 8) ) +
    scale_fill_viridis( na.value = "grey",limit = limit,space = "Lab",name = "Preference")
}


# Beautify the clonalOverlap()
# @seu_with_xcr : a seurat object which
# @cloneCall    : see clonalOverlap()
# @method       : see clonalOverlap()
# @title        : title of the plot
# @limit        : set limit of max value(value greater than this vill be white !!!)
# >>> Nothing (plotting)
clonalOverlap2 <- function(seu_with_xcr, cloneCall = 'gene+nt',method='jaccard',
                           title= 'Need Title',limit = c(0,0.5) ){
  clonalOverlap(seu_with_xcr, cloneCall=cloneCall, method=method,exportTable = T) %>% reshape2::melt('names') %>% 
    ggplot(aes(x = names,y = variable,fill = value))+
    geom_tile()+theme_bw()+
    theme_minimal()+ # 设置主题为无边框
    scale_fill_viridis( na.value = "white",limit =limit,space = "Lab",name = paste0(str_to_title(method) ,' index')) +
    labs(title = title)+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,size = 10,hjust = 1),
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
          legend.title = element_text(size = 12),
          axis.title = element_blank(), # 去除横纵坐标轴标题标签
          panel.grid.major = element_blank(), # 去除背景网格线
          panel.background = element_blank(), # 去除背景颜色
          legend.justification = c(1,0),
          legend.position = c(1.3,0.3), # 图例位置
          legend.direction = "horizontal")+ # vertical(垂直)/horizontal(水平)
    coord_fixed()+ # 确保x轴一个单位与y轴一个单位长度相同
    guides(fill = guide_colorbar(barwidth = 9,barheight = 1.5, # 图例长宽
                                 title.position = "top",title.hjust = 0.1))
}


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
                            file_path = '/data/sle/tmp/', keep_row_name = T, threads = 12){
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

# Retrun head and tail simultaneously
ht <- function(d, m=5, n=m){
  # print the head and tail together
  list(HEAD = head(d,m), TAIL = tail(d,n))
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
                    species ="Homo sapiens",category = 'H'){
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
  fgsea_results<- fgsea(fgsea_sets, stats = ranks, nperm = 10000)
  
  fgsea_results %<>% arrange (padj,desc(NES)) %>% select (pathway, padj, NES) %>%
    {rbind(head(., 8), tail(., 8))} 

  print(dim(fgsea_results))
  fgsea_results$pathway <- str_split_fixed(fgsea_results$pathway , pattern = '_',n=2)[2] %>% str_replace('_',' ') %>% str_to_title()
  View(fgsea_results)
  
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
do_harmony <- function(seu_obj = seu_obj,harmony_slot='orig.ident',theta = 2,from_begin=F,max.iter=20,
                       scale_all = F, res = 0.8, qc = F,
                      feature_exclude = '^HLA-|^IG[HJKL]|^RNA|^MT-|^RP[SL]|^TR[ABDY][VJ]',
                      reg = F,
                      vars.reg = c('nCount_RNA','pct_counts_mt'),save_tmp = F, save_path = './pca_before_harmony.rdata'){
  
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
    if(save_tmp){
      save(seu_obj,file = save_path)
    }
    
    # seu_obj <- FindNeighbors(object = seu_obj)
    # seu_obj <- FindClusters(object = seu_obj, resolution = res)
    # seu_obj <- RunUMAP(object = seu_obj,dims = 1:30)
  }
# pdf()
  seu_obj <- seu_obj %>% 
    RunHarmony(group.by.vars = harmony_slot,theta = theta, plot_convergence = F, max.iter.harmony = max.iter)%>% 
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
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab('') +
      theme_bw() +theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "grey", size=1)) + xlab('') + 
      theme(axis.text.x = element_text(angle = 90))
  )
}

# Beautify the seurat function DoHeatmap()
# @para seu_obj:     seurat object
# @para marker_list: marker genes list
DoHeatmap2 <- function(seu_obj,marker_list,group.by=NULL){
  if(length(Cells(seu_obj)) > 15000){
    print('please downsample')
  }else{
    print(
      DoHeatmap(seu_obj, marker_list,group.by = group.by)+guides(color = FALSE) +
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
