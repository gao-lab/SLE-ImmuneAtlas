# Figures

This notebook keeps the main figures of paper

```{r}
library(Seurat)
library(tidyverse)
library(ggpubr)
```

```{r, setup, include=FALSE}
knitr::opts_knit$set(root.dir = '/data/sle')
```

```{r}
source('./scripts/function_R/utils.R')
```

## Figure1

```{r}
load('./final/seurat/pbmc/04-pbmc_all_anno_modify_meta.rdata')
```

```{r}
ggplot(data = pbmc_all@meta.data, aes(x = pbmc_all$sample_name, 
                                      fill =pbmc_all$main_type))+
    geom_bar(stat = 'count',position = 'fill') + labs(y = 'proportions', x = "") + 
    scale_fill_discrete(labels= names(table(pbmc_all$main_type))) + xlab('')+ labs(fill="") + 
    scale_fill_manual(values=get_color(14 ,set = 'Paired',set_len = 9)) +
    # scale_color_npg() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                        text = element_text(size = 15,colour = 'black',face="bold"),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

```{r}

```
