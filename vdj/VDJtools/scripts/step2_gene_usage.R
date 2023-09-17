library(immunarch)

imm_10x <- repLoad('data/10x_bcr/' ,mode = "paired")
# need add a meta file in the dir but we can modify the meta slot
meta

imm_10x$meta$name  <- str_split_fixed(imm_10x$meta$Sample,pattern = '_',n=2)[,1]
imm_10x <- repFilter(imm_10x, .method = "by.meta", .query = list(name = exclude("WH1","WH2","LLS")))

imm_10x$meta %<>% left_join(meta %>% rownames_to_column(), by = c('name' = 'rowname'))

imm_ov1 <- repOverlap(imm_10x$data, .method = "public", .verbose = F)
repOverlap(imm_10x$data, .method = "public", .verbose = F) %>% vis()
repOverlap(imm_10x$data, .method = "public", .verbose = F) %>% vis('heatmap2')
repOverlapAnalysis(imm_ov1, "mds")
repOverlapAnalysis(imm_ov1, "mds+kmeans") %>% vis()

imm_gu <- geneUsage(imm_10x$data, "hs.ighd",.norm = T)
vis(imm_gu, .by = "group", .meta = imm_10x$meta)
