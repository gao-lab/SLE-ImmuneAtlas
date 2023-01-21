library(immunarch)
library(stringr)

sample = snakemake@params[[1]]
imm_10x <-repLoad(snakemake@input[[1]], .mode = 'single')
repSave(imm_10x,.format = 'vdjtools',.path = snakemake@output[[1]],.compress = F)

# single chain mode
if(str_detect(snakemake@input[[1]],"TR")) {
    # print(list.files(path=snakemake@output[[1]]),full.names=T)
    print('done')
}else{
    # double chain mode
    print('double chain mode')
    TRA <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='TRA',full.names=T),sep = '\t')
    TRB<- read.csv(list.files(path=snakemake@output[[1]] ,pattern='TRB',full.names=T),sep = '\t')

    TCR <- rbindlist(list(TRA,TRB))
    TCR <- dplyr::rename(TCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
    out_path <-  paste0(snakemake@output[[1]],'/',sample,'_tcr.tsv')
    print(out_path)
    print(head(TCR))
    write.table(TCR,file = out_path, sep = '\t', row.names = F,quote = F)
}