library(immunarch)
library(stringr)

sample = snakemake@params[[1]]
imm_10x <-repLoad(snakemake@input[[1]], .mode = 'single')
repSave(imm_10x,.format = 'vdjtools',.path = snakemake@output[[1]],.compress = F)

# single chain mode
if(str_detect(snakemake@input[[1]],"IG")) {
    # print(list.files(path=snakemake@output[[1]]),full.names=T)
    print('done')
}else{
    # double chain mode
    print('double chain mode')
    IGH <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGH',full.names=T),sep = '\t')
    IGL <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGL',full.names=T),sep = '\t')
    IGK <- read.csv(list.files(path=snakemake@output[[1]] ,pattern='IGK',full.names=T),sep = '\t')
    BCR <- rbindlist(list(IGH,IGL,IGK))
    BCR <- dplyr::rename(BCR,count=X.Seq..Count,frequency=Percent,CDR3nt=N.Sequence,
                CDR3aa=AA.Sequence,V=V.Segments,D=D.Segment,J=J.Segments)
    out_path <-  paste0(snakemake@output[[1]],'/',sample,'_bcr.tsv')
    print(out_path)
    print(head(BCR))
    write.table(BCR,file = out_path, sep = '\t', row.names = F,quote = F)
}