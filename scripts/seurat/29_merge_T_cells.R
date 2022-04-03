library(Seurat)

t_cell_merge <- merge(cd4_filter, y=c(cd8_filter,other_T_filter,nk_filter,prolife_T_fliter))
back_run(do_harmony,out_name = 't_cell_merge', job_name = 't_cell_merge',
         seu_obj = t_cell_merge, harmony_slot = 'orig.ident',max.iter = 30,res =c(1.0,0.8),from_begin = T)


DimPlot(t_cell_merge, group.by = 'subtype',label = T)

save(t_cell_merge, file = './final/seurat/t_cell/04-t_cell_merge.rdata')

