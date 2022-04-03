

library(immunarch)
imm_10x <-repLoad(snakemake@input[[1]], .mode = 'single')

repSave(imm_10x,.format = 'vdjtools',.path = snakemake@output[[1]],.compress = F)

print(list.files(path=snakemake@output[[1]]),full.names=T)

print('-----')
IGH <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGH',full.names=T),sep = '\t')
IGL <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGL',full.names=T),sep = '\t')
IGK <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGK',full.names=T),sep = '\t')

BCR <- rbindlist(list(IGH,IGL,IGK))
BCR <- dplyr::rename(BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
              CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(BCR,file = snakemake@output[[2]], sep = '\t', row.names = F,quote = F)