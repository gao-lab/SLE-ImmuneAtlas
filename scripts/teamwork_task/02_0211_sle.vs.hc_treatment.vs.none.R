
setwd('/data/sle')
source('./scripts/function_R/utils.R')
library(future)

load('./output_file/seurat/all_pbmc/modify_subtype_final_concat_pbmc_without_pSS.rdata')

pbmc_final$group %>% table()

pbmc_final$compare_group <- paste0(pbmc_final$group,'_',pbmc_final$treatment)


plan("multiprocess", workers = 16)
options(future.globals.maxSize= 200 * 1000 * 1024 ^2 )
marker_pbmc_sle.noneVShc <- FindMarkers(object = pbmc_final, ident.1 = 'SLE_pah_untreated',
                                        ident.2 = 'HC_HC', group.by = 'compare_group',only.pos = T)

back_run(func = FindMarkers,out_name ='marker_pbmc_slepah.noneVShc',job_name = 'marker_pbmc_slepah.noneVShc',
         object = pbmc_final, ident.1 = 'SLE_pah_untreated',
         ident.2 = 'HC_HC', group.by = 'compare_group',only.pos = T)

back_run(func = FindMarkers,out_name ='marker_pbmc_slepah.noneVSsle.none',job_name = 'marker_pbmc_slepah.noneVSsle.none',
         object = pbmc_final, ident.1 = 'SLE_pah_untreated',
         ident.2 = 'SLE_untreated', group.by = 'compare_group',only.pos = T)

back_run(func = FindMarkers,out_name ='marker_pbmc_slepah.noneVSslepah.treat',job_name = 'marker_pbmc_slepah.noneVSslepah.treat',
         object = pbmc_final, ident.1 = 'SLE_pah_untreated',
         ident.2 = 'SLE_pah_treated', group.by = 'compare_group',only.pos = T)


write.csv(marker_pbmc_slepah.noneVShc,file = './output_file/teamwork_task/marker_pbmc_slepah.noneVShc.csv')
write.csv(marker_pbmc_slepah.noneVSsle.none,file = './output_file/teamwork_task/marker_pbmc_slepah.noneVSsle.none.csv')
write.csv(marker_pbmc_slepah.noneVSslepah.treat,file = './output_file/teamwork_task/marker_pbmc_slepah.noneVSslepah.treat.csv')
