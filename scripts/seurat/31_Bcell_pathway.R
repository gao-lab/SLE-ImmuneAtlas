
# load( file = './final/seurat/b_cell/03-b_cell_anno_filter_harm.rdata') # bcell_filter
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)

# ---------------------------- Prepare data ------------------------------------
Idents(bcell_filter) <- 'subtype'
# we should do gsea in subtype level, first I test the IFN-response subtype
B.ifn<- subset(bcell_filter, idents = 'B.IFN-response')
B.ifn_gesa <- wilcoxauc(B.ifn , 'treatment')
# head(B.ifn_gesa)


# ---------------------------- Do enrichment -----------------------------------
# prepare fgsea data
msigdbr_collections()
msig_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- msig_df %>% split(x = .$gene_symbol, f = .$gs_name)

# rank the genes
ranks<- B.ifn_gesa %>%
    dplyr::filter(group == "untreated") %>%
    arrange(desc(logFC)) %>% 
    dplyr::select(feature, logFC) %>% deframe()

# do fgsea
fgseaRes<- fgsea(fgsea_sets, stats = ranks,  nperm = 10000, minSize = 0, maxSize =1000)
fgseaResTidy <- fgseaRes %>% as_tibble() %>%   arrange(desc(NES))

ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= NES < 7.5)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal()

# -------------------------- Copy from github  ---------------------------------
# https://github.com/scCOVID-19/COVIDPBMC/blob/master/B_cell/Analysis_and_plotting/GSEA/021921_R_NewcastleSangerCambridge_GSEA.ipynb

celltypes = c("B.transition","B.naive","B.IFN-response","B.mem.IGHM+","B.mem","B.mem.CXCR3+","B.mem.CD27-")
# celltypes = c('B.IFN-response')
groups = c('untreated', 'treated')
Idents(bcell_filter) <- 'subtype'
for (celltype in celltypes) {
    for (group in groups) {
        tmp<- subset(bcell_filter, idents = celltype)
        Idents(tmp) <- 'treatment'
        tmp <- subset(tmp, idents = c('HC',group))
        tmp_gesa <- wilcoxauc(tmp , 'treatment')
        write.table(tmp_gesa,file = paste0('/data/sle/final/gsea/bcell/',celltype,'_',group,'_gsea.csv'))
    }
}

# function
makeGeneList <- function(filename){
    gl <- read.table(filename)
    # gl <- readr::read_csv(filename)
    y <- grepl('^RPS|^RPL|^MRPL|^MRPS|^MT-', gl$feature)
    gl <- gl[!y, ]
    gl <- gl  %>% filter(!group == 'HC') %>% # not use HC
        dplyr::select(feature, logFC) %>% arrange(desc(logFC)) %>% deframe()
    return(gl)
}

plotGSEA_Hallmark <- function(gsea, group_ref = NULL, cols = NULL, newlabels = NULL, keep_significant_only = FALSE) {
    require(ggplot2)
    if(!is.null(cols)){
        gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
        }
        cols = gg_color_hue(dplyr::n_distinct(gsea$group, na.rm = TRUE))
    } else {
        cols = c('#696969','#F5DEB3','#800000')
    }   
    
    gsea$NES[which(is.na(gsea$NES))] <- 0
    gsea$pval[which(is.na(gsea$pval))] <- 1
    gsea$padj[which(is.na(gsea$padj))] <- 1
    gsea$ranking[which(is.na(gsea$ranking))] <- 0
    gsea <- gsea[order(gsea$ranking),]      
    # gsea_spl <- split(gsea, gsea$group)
    # if(!is.null(group_ref)){
    #     gsea_spl[[group_ref]] <- gsea_spl[[group_ref]][order(gsea_spl[[group_ref]]$ranking),]
    #     gsea_spl[[group_ref]]$ranking <- gsea_spl[[group_ref]]$ranking*999
    # } else {
    #     gsea_spl[[2]] <- gsea_spl[[2]]$ranking*999
    # }
    # names(gsea_spl) <- NULL
    # 
    # gsea <- do.call(rbind, gsea_spl)
    
    if (keep_significant_only){
        gseax <- split(gsea, gsea$pathway)
        for (i in 1:length(gseax)){
            if (all(gseax[[i]]$pval >= 0.05)|all(gseax[[i]]$padj >=0.25)){
                gseax[[i]] <- NA        
            }
        }
        gseax <- gseax[!is.na(gseax)]
        gsea <- do.call(rbind, gseax)
        cols = c('#F5DEB3','#800000')
    }
    # gsea <- gsea[order(gsea$ranking), ]
    gsea %<>% arrange(desc(group),-NES)
    gsea$pathway <- gsub("HALLMARK_|", "", gsea$pathway)
    gsea$group[which(gsea$pval >= 0.05 & gsea$padj >= 0.25)] <- 'NotSig'
    gsea$group[which(gsea$pval < 0.05 & gsea$padj >= 0.25)] <- 'NotSig'
    gsea$group[which(gsea$pval >= 0.05)] <- 'NotSig'
    gsea$group <- factor(gsea$group, levels = c('NotSig', 'untreated', 'treated'))
    
    x_lim_min <- abs(ceiling(min(-log10(gsea$padj))))
    x_lim_max <- abs(ceiling(max(-log10(gsea$padj))))
    
    if(x_lim_min > x_lim_max){
        xval1 <- x_lim_min * -1
        xval2 <- x_lim_min
    } else {
        xval1 <- x_lim_max * -1
        xval2 <- x_lim_max
    }
    
    g <- ggplot(gsea, aes(x = -log10(padj)*sign(NES), y = reorder(pathway, ranking), fill  = group, size = abs(NES))) + 
        scale_fill_manual(values = cols) +
        geom_point(shape = 21, colour = "black",alpha = 0.7,position = position_jitter(width = 0,height = 0)) +
        labs(x = expression(paste("Signed", " -log" ["10"], "adjusted pval")), y = "Hallmarks") +
        theme_bw() +
        geom_vline(xintercept = 0) + geom_vline(xintercept = -log10(0.25),linetype="dashed",colour = 'gray') +
        geom_vline(xintercept = -log10(0.25)*-1,linetype="dashed",colour = 'gray') + xlim(xval1, xval2) +
        scale_size_area( max_size = 7, limits = c(0,3)) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line = element_blank(), 
              axis.ticks = element_blank())
    
    g$data <- g$data[order(g$data$group, na.last = TRUE), ]
    return(g)
}

comparisons = list()
res = list()
h <- as.list(kelvinny::parse_gmt("/data/sle/source/gsea/h.all.v7.2.symbols.gmt"))
for (c in celltypes){
    for (g in groups){
        comparisons[[c]][[g]] = makeGeneList(paste0('/data/sle/final/gsea/bcell/',celltype,'_',group,'gsea.csv'))
        res[[c]][[g]] <- fgsea(pathways = h, stats = comparisons[[c]][[g]], nperm = 10000, minSize = 0, maxSize =1000) %>%
            as.data.frame()
    }
}


result = res
# names(result) <- gsub('.csv', "",files)
for(i in 1:length(comparisons)){
    result[[i]] <- lapply(result[[i]], function(x){
        x$ranking <- -log10(x$pval)*sign(x$NES) 
        x <- x[order(x$pathway), ]
        return(x)
    })
}

result <- lapply(result ,function(x){
    x[['untreated']]$group = "untreated"
    x[['treated']]$group = "treated"
    # x[['Moderate']]$group = "Moderate"
    # x[['Mild']]$group = "Mild"
    # x[['Asymptomatic']]$group = "Asymptomatic"
    return(x)
})
result2 <- lapply(result, function(x) {
    y <- do.call(rbind, x)
    return(y)
})

plotGSEA_Hallmark(result2$`B.IFN-response`,keep_significant_only = T, group_ref = "untreated") 

