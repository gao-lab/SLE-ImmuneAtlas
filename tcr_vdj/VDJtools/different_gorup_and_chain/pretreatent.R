library(data.table)
library(tidyverse)
library(magrittr)
library(ggpubr)

# ----------------------- organize data -----------------------------
meta <- read.csv('../all_sample_meta.csv', header = T)
files <- list.files('../tcr_data/',pattern = 'csv$')
HC_list <- data.frame(matrix(NA, ncol = 18, nrow = 0))
SLE_list <- data.frame(matrix(NA, ncol = 18 , nrow = 0))
for (file in files){
    name = str_split(file,'_')[[1]][1]
    # assign(name, read.csv(paste0('./tcr_data/',file),header = T))
    tmp <-  read.csv(paste0('../tcr_data/',file),header = T)
    # get(name) %<>% mutate(sample = name)
    # tmp %<>% mutate(sample = name)
    if(meta$group[which(meta$name == name)] == 'HC' ){
        HC_list <- rbind(HC_list,tmp)
    }else if(meta$group[which(meta$name == name)] == 'SLE'){
        SLE_list <- rbind(SLE_list, tmp)
    }
}

write.csv(HC_list,'./HC_tcr_list.csv',quote =F,row.names = F )
write.csv(SLE_list,'./SLE_tcr_list.csv',quote =F,row.names = F )


#-------------------------- 10x to immunarch -----------------------------------
library(immunarch)
HC_10x <-repLoad('./HC_tcr_list.csv', .mode = 'single')
SLE_10x <-repLoad('./SLE_tcr_list.csv', .mode = 'single')
repSave(HC_10x,.format = 'vdjtools',.path = './HC_vdjtools_format',.compress = F)
repSave(SLE_10x,.format = 'vdjtools',.path = './SLE_vdjtools_format',.compress = F)


#------------------------ immunarch to VDJtools --------------------------------
# SLE
SLE_TRA <- read.csv(list.files(path= './SLE_vdjtools_format/',pattern='TRA',full.names=T),sep = '\t')
SLE_TRB <- read.csv(list.files(path='./SLE_vdjtools_format/' ,pattern='TRB',full.names=T),sep = '\t')

SLE_tcr <- rbindlist(list(SLE_TRA,SLE_TRB))
SLE_tcr <- dplyr::rename(SLE_tcr,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                     CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(SLE_tcr,file = './SLE_tcr_VDJtools.txt', sep = '\t', row.names = F,quote = F)

# HC
HC_TRA <- read.csv(list.files(path= './HC_vdjtools_format/',pattern='TRA',full.names=T),sep = '\t')
HC_TRB <- read.csv(list.files(path='./HC_vdjtools_format/' ,pattern='TRB',full.names=T),sep = '\t')

HC_tcr <- rbindlist(list(HC_TRA,HC_TRB))
HC_tcr <- dplyr::rename(HC_tcr,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                         CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(HC_tcr,file = './HC_tcr_VDJtools.txt', sep = '\t', row.names = F,quote = F)


#------------------------Different chains to VDJtools---------------------------
SLE_TRA  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                         CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
SLE_TRB  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
HC_TRA  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
HC_TRB  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)

write.table(SLE_TRA,file = './SLE_TRA_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(SLE_TRB,file = './SLE_TRB_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(HC_TRA,file = './HC_TRA_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(HC_TRB,file = './HC_TRB_VDJtools.txt', sep = '\t', row.names = F,quote = F)


#-------------------------- different TR class ratio ---------------------------
# HC_IgClass <- data.frame(matrix(NA, ncol = 18 + 1, nrow = 0))
# SLE_IgClass <- data.frame(matrix(NA, ncol = 18 + 1, nrow = 0))
# for (file in files){
#     name = str_split(file,'_')[[1]][1]
#     # assign(name, read.csv(paste0('./tcr_data/',file),header = T))
#     tmp <-  read.csv(paste0('./tcr_data/',file),header = T)
#     # get(name) %<>% mutate(sample = name)
#     tmp %<>% mutate(sample = name)
#     if(meta$group[which(meta$name == name)] == 'HC' ){
#         HC_IgClass <- rbind(HC_IgClass,tmp)
#     }else if(meta$group[which(meta$name == name)] == 'SLE'){
#         SLE_IgClass <- rbind(SLE_IgClass, tmp)
#     }
# }

# HC_ig.df <- HC_IgClass %>% group_by(sample) %>% count(c_gene, sample) %>% 
#     filter(str_detect(c_gene, "^IGH")) %>% group_by(sample) %>% mutate(ratio = n/sum(n) * 100) 
# ggplot(HC_ig.df, aes(x = sample,y=n, fill = c_gene))+
#     geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
#     # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
#     # scale_fill_manual("legend", values = ) +  # change the color
#     labs(fill="IG class")

# SLE_ig.df <- SLE_IgClass %>% group_by(sample) %>% count(c_gene, sample) %>% 
#     filter(str_detect(c_gene, "^IGH")) %>% group_by(sample) %>% mutate(ratio = n/sum(n) * 100)

# ggplot(SLE_ig.df,aes(x = sample,y=n, fill = c_gene))+
#     geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
#     # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
#     # scale_fill_manual("legend", values = ) +  # change the color
#     labs(fill="IG class")

# rbind(HC_ig.df,SLE_ig.df) %>% ggplot(aes(x = sample,y=n, fill = c_gene))+
#     geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
#     # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
#     # scale_fill_manual("legend", values = ) +  # change the color
#     labs(fill="IG class")

# ggplot(rbind(HC_ig.df,SLE_ig.df), aes(x=c_gene,y=sample)) + geom_tile(aes(fill=ratio)) + 
#     scale_fill_gradient(low = "white", high = "red")+ xlab("samples") + 
#     theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank())

# all_ig.mat <- rbind(HC_ig.df,SLE_ig.df) %>% select(-n) %>% spread(c_gene,ratio,fill =0 ) %>% 
#     column_to_rownames('sample')
#     # left_join(meta)
# heatmap(all_ig.mat %>% as.matrix(),scale = 'column')


# all_ig.df <-rbind(HC_ig.df,SLE_ig.df) %>% left_join(meta,by=c('sample'='name')) %>%
#     mutate(tmp=paste0(group,'_', treatment)) %>% filter(c_gene != 'IGHE')
# my_comparisons = list(c('HC_HC','SLE_none'),c('SLE_none','SLE_treatment'),c('HC_HC','SLE_treatment'))
# ggpubr::ggboxplot(all_ig.df, x='tmp',y='ratio', color = 'tmp',
#                   palette =c("#FC4E07", "#00AFBB","#FFC0CB"))+ 
#     facet_wrap(~c_gene,scales = "free",nrow = 2) +
#     stat_compare_means(method = 't.test',comparisons=my_comparisons)
