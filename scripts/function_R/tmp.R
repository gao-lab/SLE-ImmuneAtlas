

bcell_harm$RNA_snn_res.0.8
bcell_harm$RNA_snn_res.1

library(ggplot2)
library(ggalluvial)
library(data.table)
library(magrittr)


temp.dt <- rbindlist(list(
    bcell_harm$RNA_snn_res.0.8 %>% {data.table(res="res_0.8", cell=names(.), cluster=.)},
    bcell_harm$RNA_snn_res.1 %>% {data.table(res="res_1", cell=names(.), cluster=.)}
))

ggplot(temp.dt, aes(x=res, stratum=cluster, alluvium=cell, fill=cluster, label=cluster)) +
    geom_flow() +
    geom_stratum()

## plasma

temp.dt <- lapply(
    paste(sep="", "RNA_snn_res.", c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)),
    FUN=function(temp.name){
        plasma_filter_harmony[[temp.name]] %>% {data.table(res=temp.name, cell=rownames(.), cluster=.[, 1])}
    }
) %>% rbindlist

date()
ggplot(temp.dt, aes(x=res, stratum=cluster, alluvium=cell, fill=cluster, label=cluster)) +
    geom_flow() +
    geom_stratum()
date()

