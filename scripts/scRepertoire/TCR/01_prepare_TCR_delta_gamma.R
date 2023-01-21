library(scRepertoire)
library(stringr)
# library(immunarch)
library(Seurat)

################################################################################
#
# This script had been desert becasue we do not find delta gamma TCR
#
################################################################################

stop('you are running desert script')

setwd('/data/sle/')
# source()
output_path <- './final/scRepertotre/TCR/'

#------------------------------ Read File --------------------------------------
tcr_sample <- list.files('./data/10x_tcr/',,pattern = 'csv$')
tcr_sample_name <- str_split_fixed(tcr_sample,pattern = '_',2)[,1]

tcr_list <- list()
for(i in tcr_sample){
    assign( i,as.data.frame(read.csv(paste0('./data/10x_tcr/',i)))) 
}


#------------------------------ Add Meta ---------------------------------------
tcr_id = c("SLE", "SLE", "SLE","SLE", "SLE", "SLE","SLE","SLE","SLE","HC",
           "SLE", "SLE","SLE","SLE","SLE","SLE","SLE","HC",'HC','SLE',
           "SLE","HC")
tcr_group = c("none", "none","before","after", "none", 
              "none", "before", "after","none", "HC",
              "none", "before","after","none", "before",
              "after", "after", "HC", "HC","before", 
              "after", "HC")

contig_list_tcr <- list(GW_tcr.csv,GZR_tcr.csv,HXR_tcr.csv,HXR2_tcr.csv,HXX_tcr.csv,
                        LGY_tcr.csv,LL_tcr.csv,LL2_tcr.csv,MXY_tcr.csv,QJY_tcr.csv,
                        SQ_tcr.csv,WYF_tcr.csv, WYF2_tcr.csv,WYY_tcr.csv,XH_tcr.csv,
                        XYY_tcr.csv,XYY2_tcr.csv,ZH_tcr.csv,ZMY1_tcr.csv,ZPP_tcr.csv,
                        ZPP2_tcr.csv,ZS_tcr.csv)

# combined_tcr <- combineTCR(contig_list_tcr,
#                            samples = tcr_sample_name, 
#                            ID = tcr_id)
back_run(combineTCR, out_name = 'combined_tcr',job_name = 'combined_tcr',
         contig_list_tcr,samples = tcr_sample_name, ID = tcr_id)

# rm(list(tcr_sample,tcr_sample_name,tcr_id,tcr_group,contig_list))


#------------------------------ Save File --------------------------------------
save(combined_tcr ,file = 'final/scRepertoire/TCR/combined_tcr_delta_gamma.rdata')

