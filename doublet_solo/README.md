# 用solo定义doublet

solo是一个利用半监督学习来寻找single cell RNA 数据中doublet的算法，在多个benchmark中都优于常用的DoubletFinfer和Scrublet算法，所以想用solo来定义数据集中的doublet，并寻找这些可能的doublet和TLB之间的关系

### work on smartseq2 data

基于smartseq2的实验原理，在正常情况下不容易出现doublet，为了测试算法的性能，我们在6个包含/不包含tlb的smartseq2数据集上运行了solo，使用默认的参数，运行结果如下

| dataset   | cell_number | doublet_number |
| --------- | ----------- | -------------- |
| GSE103322 | 5902        | 0              |
| GSE115978 | 7186        | 0              |
| GSE120575 | 16289       | 2              |
| GSE146771 | 10468       | 0              |
| GSE72056  | 4645        | 0              |
| GSE140228 | 7074        | 0              |

说明smart-seq2中的TLB不是由于实验造成的假象

### work on 10X data

我们从10x官网上下载了全部的PBMC的单细胞测序数据，涵盖了10x公司不同的化学试剂、测序平台、分选仪器的数据，首先用smart-seq2使用同样的流程进行doublet分析

但是考虑到10x工作原理中容易产生一定比率的doublet，按照官方的说法doublet的比率在5%以下，所以我们将solo中的`-e EXPECTED_NUMBER_OF_DOUBLETS`参数的值设置为数据集大小的5%（默认为None），其他参数不变，重新执行了分析

两次分析的结果都如下表所示

| dataset                           | cell_number | doublet_number_without_expect | doublet_ratio_without_expect | doublet_number_with_expect | doublet_ratio_with_expect |
| --------------------------------- | ----------- | ----------------------------- | ---------------------------- | -------------------------- | ------------------------- |
| PBMC_10k_dual_indexed             | 9972        | 459                           | 4.60                         | 498                        | 4.99                      |
| PBMC_10k_single_indexed           | 10008       | 385                           | 3.85                         | 500                        | 5.00                      |
| PBMC_10k_surface_protein          | 7865        | 408                           | 5.19                         | 393                        | 5.00                      |
| PBMC_10k_v3                       | 10959       | 616                           | 5.62                         | 547                        | 4.99                      |
| PBMC_33k                          | 33052       | 1412                          | 4.27                         | 1652                       | 5.00                      |
| PBMC_3k                           | 2668        | 10                            | 0.37                         | 133                        | 4.99                      |
| PBMC_4k                           | 4244        | 17                            | 0.40                         | 212                        | 5.00                      |
| PBMC_50_50_mixture                | 8040        | 80                            | 1.00                         | 402                        | 5.00                      |
| PBMC_5k_chromium_connect          | 3642        | 19                            | 0.52                         | 182                        | 5.00                      |
| PBMC_5k_manual                    | 4027        | 12                            | 0.30                         | 201                        | 4.99                      |
| PBMC_5k_Next_GEM                  | 4563        | 44                            | 0.96                         | 228                        | 5.00                      |
| PBMC_5k_Next_GEM_2                | 4580        | 17                            | 0.37                         | 229                        | 5.00                      |
| PBMC_5k_Next_GEM_surface_protein  | 5247        | 84                            | 1.60                         | 262                        | 4.99                      |
| PBMC_5k_Next_GEM2_surface_protein | 5527        | 58                            | 1.05                         | 276                        | 4.99                      |
| PBMC_6k                           | 5352        | 7                             | 0.13                         | 267                        | 4.99                      |
| PBMC_8k                           | 6027        | 0                             | 0.00                         | 301                        | 4.99                      |
| PBMC_90_10_mixture                | 6945        | 123                           | 1.77                         | 347                        | 5.00                      |
| PBMC_99_1_mixture                 | 8185        | 140                           | 1.71                         | 409                        | 5.00                      |
| PBMC_donor_A_68k                  | 68389       | 2714                          | 3.97                         | 3419                       | 5.00                      |
| PBMC_donor_B                      | 7647        | 58                            | 0.76                         | 382                        | 5.00                      |
| PBMC_donor_C                      | 9302        | 145                           | 1.56                         | 465                        | 5.00                      |
| PBMC_whole_transcriptome          | 9365        | 385                           | 4.11                         | 468                        | 5.00                      |

### 10x中定义的doublet与TLB关系

ongoing...