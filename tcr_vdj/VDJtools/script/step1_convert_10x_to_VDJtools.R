

library(immunarch)
imm_10x <-repLoad(snakemake@input[[1]], .mode = 'single')

repSave(imm_10x,.format = 'vdjtools',.path = snakemake@output[[1]],.compress = F)

print(list.files(path=snakemake@output[[1]]),full.names=T)

print('-----')
TRA <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='TRA',full.names=T),sep = '\t')
TRB<- read.csv(list.files(path=snakemake@output[[1]] ,pattern='TRB',full.names=T),sep = '\t')


TCR <- rbindlist(list(TRA,TRB))
TCR <- dplyr::rename(TCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
              CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
write.table(TCR,file = snakemake@output[[2]], sep = '\t', row.names = F,quote = F)