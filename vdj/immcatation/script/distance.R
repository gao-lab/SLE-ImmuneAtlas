library(shazam)
library(tidyverse)

# ----------------------- read data ---------------------------
sample <- read.csv(snakemake@input[[1]], header=T, sep='\t') %>% as.tibble()


# --------------------- remove duplicate -----------------------
sample_index <- sample %>% group_by(cell_id) %>% summarise(umi_count=max(umi_count))
sample_filter <- sample[which(paste0(sample$cell_id,sample$umi_count) %in% paste0(sample_index$cell_id,sample_index$umi_count) ),]
# sample_filter %>% dim()

dup_index <- (table(sample_filter$cell_id) >1) %>% as.data.frame() %>% filter(.== TRUE) %>% rownames()
sample_single <-sample_filter[ !sample_filter$cell_id %in% dup_index, ]


# -------------------- calculate distance -----------------------
dist_sc <- distToNearest(sample_single, cellIdColumn="cell_id", locusColumn="locus", VJthenLen=FALSE, onlyHeavy=FALSE, nproc =4)

output2 <- findThreshold(dist_sc$dist_nearest, method="gmm", model="gamma-gamma")

pdf(snakemake@output[[2]])
plot(output2, binwidth=0.02, title="GMM Method: gamma-gamma")
dev.off()


# -------------------- calculate distance -----------------------
distance <- output2@threshold
write(distance, snakemake@output[[1]])