do_seurat <- function(seu_obj, scale_all = F, res = 0.8, qc = F,
                      var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT-|^RP[SL]|^TR[ABDY][VJ]',
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
  
  var_regex = var_regex # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
  hvg = grep(var_regex, hvg, invert=T, value=T)
  
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
