library(data.table)
library(tidyverse)
library(magrittr)
library(ggpubr)

meta <- read.csv('./vdj/VDJtools/all_sample_meta.csv', header = T)
files <- list.files('./vdj/VDJtools/bcr_data/',pattern = 'csv$')
HC_list <- data.frame(matrix(NA, ncol = 18, nrow = 0))
SLE_list <- data.frame(matrix(NA, ncol = 18 , nrow = 0))
for (file in files){
    name = str_split(file,'_')[[1]][1]
    # assign(name, read.csv(paste0('./vdj/VDJtools/bcr_data/',file),header = T))
    tmp <-  read.csv(paste0('./vdj/VDJtools/bcr_data/',file),header = T)
    # get(name) %<>% mutate(sample = name)
    # tmp %<>% mutate(sample = name)
    if(meta$group[which(meta$name == name)] == 'HC' ){
        HC_list <- rbind(HC_list,tmp)
    }else if(meta$group[which(meta$name == name)] == 'SLE'){
        SLE_list <- rbind(SLE_list, tmp)
    }
}

write.csv(HC_list,'./vdj/VDJtools/different_Ig_class/HC_bcr_list.csv',quote =F,row.names = F )
write.csv(SLE_list,'./vdj/VDJtools/different_Ig_class/SLE_bcr_list.csv',quote =F,row.names = F )


#-------------------------- 10x to immunarch -----------------------------------
library(immunarch)
HC_10x <-repLoad('./vdj/VDJtools/different_Ig_class/HC_bcr_list.csv', .mode = 'single')
SLE_10x <-repLoad('./vdj/VDJtools/different_Ig_class/SLE_bcr_list.csv', .mode = 'single')
repSave(HC_10x,.format = 'vdjtools',.path = './vdj/VDJtools/different_Ig_class/HC_vdjtools_format',.compress = F)
repSave(SLE_10x,.format = 'vdjtools',.path = './vdj/VDJtools/different_Ig_class/SLE_vdjtools_format',.compress = F)


#------------------------ immunarch to VDJtools --------------------------------
# SLE
SLE_IGH <- read.csv(list.files(path= './vdj/VDJtools/different_Ig_class/SLE_vdjtools_format/',pattern='IGH',full.names=T),sep = '\t')
SLE_IGL <- read.csv(list.files(path='./vdj/VDJtools/different_Ig_class/SLE_vdjtools_format/' ,pattern='IGL',full.names=T),sep = '\t')
SLE_IGK <- read.csv(list.files(path='./vdj/VDJtools/different_Ig_class/SLE_vdjtools_format/' ,pattern='IGK',full.names=T),sep = '\t')

SLE_BCR <- rbindlist(list(SLE_IGH,SLE_IGL,SLE_IGK))
SLE_BCR <- dplyr::rename(SLE_BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                     CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(SLE_BCR,file = './vdj/VDJtools/different_Ig_class/SLE_BCR_VDJtools.txt', sep = '\t', row.names = F,quote = F)

# HC
HC_IGH <- read.csv(list.files(path= './vdj/VDJtools/different_Ig_class/HC_vdjtools_format/',pattern='IGH',full.names=T),sep = '\t')
HC_IGL <- read.csv(list.files(path='./vdj/VDJtools/different_Ig_class/HC_vdjtools_format/' ,pattern='IGL',full.names=T),sep = '\t')
HC_IGK <- read.csv(list.files(path='./vdj/VDJtools/different_Ig_class/HC_vdjtools_format/' ,pattern='IGK',full.names=T),sep = '\t')

HC_BCR <- rbindlist(list(HC_IGH,HC_IGL,HC_IGK))
HC_BCR <- dplyr::rename(HC_BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                         CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(HC_BCR,file = './vdj/VDJtools/different_Ig_class/HC_BCR_VDJtools.txt', sep = '\t', row.names = F,quote = F)


#------------------------Different chains to VDJtools---------------------------
SLE_IGH  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                         CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
SLE_IGL  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
SLE_IGK  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
HC_IGH  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
HC_IGL  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
HC_IGK  %<>%  dplyr::rename(count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                             CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(SLE_IGH,file = './vdj/VDJtools/different_Ig_class/SLE_IGH_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(SLE_IGL,file = './vdj/VDJtools/different_Ig_class/SLE_IGL_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(SLE_IGK,file = './vdj/VDJtools/different_Ig_class/SLE_IGK_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(HC_IGH,file = './vdj/VDJtools/different_Ig_class/HC_IGH_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(HC_IGL,file = './vdj/VDJtools/different_Ig_class/HC_IGL_VDJtools.txt', sep = '\t', row.names = F,quote = F)
write.table(HC_IGK,file = './vdj/VDJtools/different_Ig_class/HC_IGK_VDJtools.txt', sep = '\t', row.names = F,quote = F)


#-------------------------- different Ig class ratio ---------------------------
HC_IgClass <- data.frame(matrix(NA, ncol = 18 + 1, nrow = 0))
SLE_IgClass <- data.frame(matrix(NA, ncol = 18 + 1, nrow = 0))
for (file in files){
    name = str_split(file,'_')[[1]][1]
    # assign(name, read.csv(paste0('./vdj/VDJtools/bcr_data/',file),header = T))
    tmp <-  read.csv(paste0('./vdj/VDJtools/bcr_data/',file),header = T)
    # get(name) %<>% mutate(sample = name)
    tmp %<>% mutate(sample = name)
    if(meta$group[which(meta$name == name)] == 'HC' ){
        HC_IgClass <- rbind(HC_IgClass,tmp)
    }else if(meta$group[which(meta$name == name)] == 'SLE'){
        SLE_IgClass <- rbind(SLE_IgClass, tmp)
    }
}

HC_ig.df <- HC_IgClass %>% group_by(sample) %>% count(c_gene, sample) %>% 
    filter(str_detect(c_gene, "^IGH")) %>% group_by(sample) %>% mutate(ratio = n/sum(n) * 100) 
ggplot(HC_ig.df, aes(x = sample,y=n, fill = c_gene))+
    geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
    # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
    # scale_fill_manual("legend", values = ) +  # change the color
    labs(fill="IG class")

SLE_ig.df <- SLE_IgClass %>% group_by(sample) %>% count(c_gene, sample) %>% 
    filter(str_detect(c_gene, "^IGH")) %>% group_by(sample) %>% mutate(ratio = n/sum(n) * 100)

ggplot(SLE_ig.df,aes(x = sample,y=n, fill = c_gene))+
    geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
    # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
    # scale_fill_manual("legend", values = ) +  # change the color
    labs(fill="IG class")

rbind(HC_ig.df,SLE_ig.df) %>% ggplot(aes(x = sample,y=n, fill = c_gene))+
    geom_bar(stat = 'identity',position = 'fill') + labs(y = 'proportions', x = "") + 
    # scale_fill_discrete(names(table(Ig_ratio$c_gene))) + 
    # scale_fill_manual("legend", values = ) +  # change the color
    labs(fill="IG class")

ggplot(rbind(HC_ig.df,SLE_ig.df), aes(x=c_gene,y=sample)) + geom_tile(aes(fill=ratio)) + 
    scale_fill_gradient(low = "white", high = "red")+ xlab("samples") + 
    theme_bw() + theme(panel.grid.major = element_blank()) + theme(legend.key=element_blank())

all_ig.mat <- rbind(HC_ig.df,SLE_ig.df) %>% select(-n) %>% spread(c_gene,ratio,fill =0 ) %>% 
    column_to_rownames('sample')
    # left_join(meta)
heatmap(all_ig.mat %>% as.matrix(),scale = 'column')


all_ig.df <-rbind(HC_ig.df,SLE_ig.df) %>% left_join(meta,by=c('sample'='name')) %>%
    mutate(tmp=paste0(group,'_', treatment)) %>% filter(c_gene != 'IGHE')
my_comparisons = list(c('HC_HC','SLE_none'),c('SLE_none','SLE_treatment'),c('HC_HC','SLE_treatment'))
ggpubr::ggboxplot(all_ig.df, x='tmp',y='ratio', color = 'tmp',
                  palette =c("#FC4E07", "#00AFBB","#FFC0CB"))+ 
    facet_wrap(~c_gene,scales = "free",nrow = 2) +
    stat_compare_means(method = 't.test',comparisons=my_comparisons)
