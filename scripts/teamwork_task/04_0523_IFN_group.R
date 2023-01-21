# Need seurat and pbmc_all
library(Seurat)
library(reshape2)
library(assert)

SLEDAI <- c(6,2,4,1,11,4,7,4,3,2,2,4,7,8,5,8,2,2,0,0,0,0)
names(SLEDAI) <- c('XH','XYY','XYY2','GW','WYF','WYF2','HXR','HXR2','ZPP','ZPP2',
                   'LGY','HXX','GZR','MXY','SQ','LL','LL2','WYY','QJY','ZMY1','ZH','ZS')
SLEDAI <- as.data.frame(SLEDAI)
meta <- unique(pbmc_all@meta.data[,c('treatment','orig.ident','pair')]) %>% remove_rownames() %>%column_to_rownames('orig.ident')
meta <-left_join(rownames_to_column(meta),rownames_to_column(SLEDAI),
                 by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
################################################################################
#
#  Use total ratio of PBMC to cluster
#
################################################################################
tmp <- table(pbmc_all$orig.ident, pbmc_all$subtype) %>% data.frame()
tmp <- reshape2::dcast(tmp,formula = Var1~ Var2, value.var = 'Freq' ) %>% 
    column_to_rownames('Var1') %>% mutate(all = rowSums(.)) 
for(i in c(1:dim(tmp)[1])){
    for(j in c(1:(dim(tmp)[2]-1))){
        # print(tmp[i,j])
        tmp[i,j] <- tmp[i,j]/tmp[i,dim(tmp)[2]]
    }
}
tmp <- select(tmp, -all)
rowSums(tmp)
tmp2 <-left_join(rownames_to_column(tmp),rownames_to_column(SLEAI),
                         by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
# cor(tmp2) %>% View()

b_sum <- table(bcell_filter$orig.ident) %>% as.data.frame()
cd4_sum <- table(cd4_filter$orig.ident) %>% as.data.frame()
cd8_sum <- table(cd8_filter$orig.ident) %>% as.data.frame()
b_sum$Var1 == cd4_sum$Var1
cd4_sum$Var1 == cd8_sum$Var1
tmp <- tmp[match(b_sum$Var1 %>% as.vector(), rownames(tmp)),]
rownames(tmp) == cd4_sum$Var1

tmp$`B.IFN-response` <- tmp$`B.IFN-response`/b_sum$Freq
tmp$`T.CD4.IFN-response` <- tmp$`T.CD4.IFN-response`/cd4_sum$Freq
tmp$`T.CD8.IFN-response` <- tmp$`T.CD8.IFN-response`/cd8_sum$Freq

# pheatmap(tmp, scale = 'row')
# pheatmap(tmp, scale = 'column')
# pheatmap(tmp, scale = 'none')

ifn_ratio_df <- tmp[,colnames(tmp) %>% str_detect( 'IFN')] * 100
head(ifn_ratio_df)

# pheatmap(ifn_ratio_df,scale = 'column')
# pheatmap(ifn_ratio_df, scale = 'column', cutree_rows =2)
# # pheatmap(ifn_ratio_df, scale = 'row') # not reasonable
# pheatmap(ifn_ratio_df, scale = 'none', cutree_rows =2)

ifn_ratio_df <-left_join(rownames_to_column(ifn_ratio_df),rownames_to_column(SLEAI),
                             by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')

cor.test(ifn_ratio_df[,1],ifn_ratio_df[,4])
cor.test(ifn_ratio_df[,2],ifn_ratio_df[,4])
cor.test(ifn_ratio_df[,3],ifn_ratio_df[,4]) # 0.03791

meta <- unique(pbmc_all@meta.data[,c('treatment','orig.ident')]) %>% remove_rownames() %>%column_to_rownames('orig.ident')
ifn_ratio_df <-left_join(rownames_to_column(ifn_ratio_df),rownames_to_column(meta),
                         by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')

scatterplotMatrix(ifn_ratio_df[4:1],groups =ifn_ratio_df$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='IFN-response and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = F)

fit<-lm(sledai~`B.IFN-response` + `T.CD8.IFN-response` + `T.CD4.IFN-response` ,data=ifn_ratio_df)
summary(fit)
plot(fit)

write.csv(ifn_ratio_df, file = './tmp/ifn_ratio_df.csv')

################################################################################
#
#  + and - correlation 
#
################################################################################
tmp2 <-left_join(rownames_to_column(tmp2),rownames_to_column(meta),
                 by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
# positive cor: T.prolife Treg plasmablast

scatterplotMatrix(tmp2[,c('T.prolife','plasmablast','SLEAI.x')],groups =tmp2$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='IFN-response and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )

# nega cor: memory
scatterplotMatrix(tmp2[,c('B.mem','SLEAI')],groups =tmp2$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='IFN-response and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = list(coords= "topleft", cex = 0.01,  text.width =0.01) )


################################################################################
#
#  Use cell type specific ratio of PBMC to cluster
#
################################################################################
# we should get the IFN-response ratio in B cell, CD4 and CD8 respectively

ifn_b_ratio <-  bcell_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(bcell_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='B.IFN-response') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname')) 
ifn_cd4_ratio <-  cd4_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd4_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='T.CD4.IFN-response') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname')) 
ifn_cd8_ratio <-  cd8_filter@meta.data  %>% 
    group_by(orig.ident,subtype) %>% summarise(sub_num = n()) %>% 
    mutate(sample_num = sum(sub_num)) %>% mutate(Ratio = sub_num/sample_num*100) %>%
    left_join(cd8_filter@meta.data[,c(1,4)]  %>%  distinct()) %>%
    filter(subtype =='T.CD8.IFN-response') %>% left_join(meta %>% rownames_to_column(), by = c('orig.ident' = 'rowname')) 
stopifnot(ifn_b_ratio$orig.ident == ifn_cd4_ratio$orig.ident && ifn_b_ratio$orig.ident== ifn_cd8_ratio$orig.ident)
ifn_sub_ratio <- data.frame(row.names = ifn_b_ratio$orig.ident,b_ifn = ifn_b_ratio$Ratio,
                            cd4_ifn = ifn_cd4_ratio$Ratio, cd8_ifn = ifn_cd8_ratio$Ratio)

# pheatmap(ifn_sub_ratio[,1:3], scale = 'row')
# pheatmap(ifn_sub_ratio[,1:3], scale = 'column')
# pheatmap(ifn_sub_ratio[,1:3], scale = 'none')

ifn_sub_ratio <- left_join(rownames_to_column(ifn_sub_ratio),rownames_to_column(meta),
          by = c('rowname' = 'rowname')) %>% column_to_rownames('rowname')
cor(ifn_sub_ratio[c(1,2,3,6)])
library(car)
# scatterplotMatrix(ifn_sub_ratio)
scatterplotMatrix(ifn_sub_ratio[c(1,2,3,6)],groups =ifn_sub_ratio$treatment ,
                  smooth = list(spread = T,  lty.smooth=2, lwd.smooth=3, lty.spread=3, lwd.spread=2),
                  spread=FALSE,main='IFN-response and SLEDAI', by.groups = F,
                  pch = c(15,16,17),col = c("#88a16f", "#9DB0D3", "#DA9494"),  diagonal=list(method ="boxplot"), 
                  legend = F )
cor.test(ifn_sub_ratio[,1],ifn_sub_ratio[,4])
cor.test(ifn_sub_ratio[,2],ifn_sub_ratio[,4])
cor.test(ifn_sub_ratio[,3],ifn_sub_ratio[,4]) # 0.03791

fit<-lm(sledai~b_ifn +cd4_ifn +cd8_ifn ,data=ifn_sub_ratio)
summary(fit)


################################################################################
#
#  Use Random Forrest(RF) to choose the better feature about SLEDAI
#
################################################################################
# we should build a table, which every row is a sample and every col is a feature
# features include cell type ratio and IFN-related gene expression score





