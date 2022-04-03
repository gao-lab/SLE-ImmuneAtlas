
# since vdjtools not support the 10x data, I try use immunarch to do transition 

setwd('/data/sle')
library(immunarch)
imm_10x <- repLoad('./vdj/VDJview/test_data/GW_10X_bcr/filtered_contig_annotations.csv', .mode = 'single')

repSave(imm_10x,.format = 'vdjtools',.path = './vdj/VDJtools/GW_bcr_vdj',.compress = F)


GW_IGH <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGH.tsv',sep = '\t')
GW_IGL <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGL.tsv',sep = '\t')
GW_IGK <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGK.tsv',sep = '\t')

GW_BCR <- rbindlist(list(GW_IGH,GW_IGK,GW_IGL))
GW_BCR <- dplyr::rename(GW_BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
              CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(GW_BCR,file = './vdj/VDJtools/GW_bcr_vdj/gw_bcr.tsv', sep = '\t', row.names = F,quote = F)


# test LL which file in snakemake

library(immunarch)
imm_10x <- repLoad('./bcr_data/LL_bcr.csv', .mode = 'single')

repSave(imm_10x,.format = 'vdjtools',.path = './vdj/VDJtools/GW_bcr_vdj',.compress = F)


GW_IGH <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGH.tsv',sep = '\t')
GW_IGL <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGL.tsv',sep = '\t')
GW_IGK <- read.csv('./vdj/VDJtools/GW_bcr_vdj/filtered_contig_annotations_IGK.tsv',sep = '\t')

GW_BCR <- rbindlist(list(GW_IGH,GW_IGK,GW_IGL))
GW_BCR <- dplyr::rename(GW_BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
              CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(GW_BCR,file = './vdj/VDJtools/GW_bcr_vdj/gw_bcr.tsv', sep = '\t', row.names = F,quote = F)
