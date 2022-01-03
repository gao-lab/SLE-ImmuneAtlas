### GSE134355

 郭国骥老师的 human cell landscape (HCL)项目，采用的是Microwell-seq，在141个检测的组织和器官，包含了成年人和胎儿的大部分重要器官，也有少量小鼠的数据（几乎所有器官都有重复实验）中，在8个数据集中发现了tlb分别如下

| 数据集                              | 重复实验次数    | 占B细胞比例(%) |
| ----------------------------------- | --------------- | -------------- |
| Adult-Epityphlon1                   | 1(1个中发现tlb) | 7.25           |
| Adult-Spleen1-2                     | 2(1个中发现tlb) | 2.1            |
| Adult-Rectum1                       | 1(1个中发现tlb) | 10.1           |
| Adult-Transverse-Colon2-1           | 3(2个中发现tlb) | 7.3            |
| Adult-Transverse-Colon2-2           | 3(2个中发现tlb) | 1.6            |
| Fetal-Pancreas1                     | 3(1个中发现tlb) | 53.8           |
| Adult-Liver4-2                      | 5(1个中发现tlb) | 27.2           |
| Adult-Ileum2                        | 1(1个中发现tlb) | 20.7           |
| 注意，其他没有发现TLB的组织没有列出 |                 |                |

与之前的结果相比，

1. 相同：在肠道、脾脏、肝脏发现了tlb
2. 不同：在骨髓、外周血没有发现tlb；此外，在胎儿的组织中发现的tlb明显少，只在Pancreas的一个重复中发现，并且tlb的比例特别高

### GSE135779

发表在**nat immunity**上的sle文章，之前已经对该文章进行了详细的分析，没有在不同group之间找到显著的差异

可能有一些混杂因素需要去掉、回归

- 测序深度的不同
- 不同group样本数量相差较大

### GSE139324

HNSCC数据集，一共有63个样本

| 分组       | 重复数量 | 有TLB的样本的数量 |
| ---------- | -------- | ----------------- |
| HNSCC_PBMC | 26       | 0                 |
| HNSCC_TIL  | 26       | 4                 |
| HD_PBMC    | 6        | 0                 |
| HD_Tonsil  | 5        | 4                 |
| 总计63     | 63       | 8                 |

这个数据集的分布同样不均匀，但是发现对单个样本的 PBMC进行的分析均没有发现TLB。可能与TLB在外周血中的比例过低有关。                                                                                                                                                                                        

### GSE139555

来自 14 名癌症患者的预处理样本的单细胞 RNA-seq 和单细胞 TCR-seq，取样自肿瘤、正常邻近组织和外周血，一共32个样本

| 分组               | 重复数量 | 有TLB的样本的数量 |
| ------------------ | -------- | ----------------- |
| Lung tumor         | 6        | 2                 |
| Lung  normal       | 6        | 1                 |
| Lung PBMC          | 1        | 0                 |
| Endometrial tumor  | 3        | 0                 |
| Endometrial normal | 3        | 0                 |
| Colorectal tumor   | 2        | 1                 |
| Colorectal normal  | 2        | 0                 |
| Renal tumor        | 3        | 0                 |
| Renal normal       | 3        | 1                 |
| Renal  PBMC        | 3        | 0                 |
| 总计               | 32       | 5                 |

同样的，没有发现tlb分布的明显偏好性。且在pbmc中没有发现

### GSE154109

健康人、儿童 B 细胞急性淋巴细胞白血病患者(ALL)、急性髓系白血病患者(AML)的骨髓

| 分组               | 重复数量 | 有TLB的样本的数量 |
| ------------------ | -------- | ----------------- |
| health bone marrow | 4        | 4                 |
| AML bone marrow    | 8        | 0                 |
| ALL bone marrow    | 7        | 3                 |
| 总计               | 19       | 7                 |

值得注意的是在该数据集所有的健康人骨髓中都发现了tlb

### GSE161529

69 个 scRNA-seq 样本，有 52 名患者：包括 4 个 TNBC、4 个 BRCA1 TNBC、6 个 HER+ 肿瘤、19 个 ER+ 肿瘤和 6 个 ER+ 肿瘤的淋巴结转移的概况。它还包括来自 13 名没有乳腺癌的正常患者的总乳腺细胞、来自 11 名正常患者的上皮乳腺细胞以及来自 4 名具有 BRCA1 突变的癌前患者的总乳腺细胞的概况。

| 分组                                     | 重复数量 | 有TLB的样本的数量 |
| ---------------------------------------- | -------- | ----------------- |
| Normal Total cells                       | 13       | 0                 |
| Normal Epithelial cells                  | 11       | 0                 |
| BRCA1 pre-neoplastic Total cells         | 4        | 0                 |
| Triple negative tumour Total cells       | 4        | 0                 |
| Triple negative BRCA1 tumour Total cells | 4        | 1                 |
| HER2+ tumour Total cells                 | 6        | 4                 |
| ER+ tumour Total cells                   | 19       | 1                 |
| ER+ tumour Lymph-node cells              | 7        | 4                 |
| PR+ tumour Total cells                   | 1        | 0                 |
| 总计                                     | 69       | 10                |

分布的情况和数据集本身是否包含免疫细胞关系也挺大的

### GSE163314

研究克罗恩病和脊椎关节炎这两种疾病的患者的结肠和外周血中的免疫学特征的相关性

| 分组                    | 重复数量 | 有TLB的样本的数量 |
| ----------------------- | -------- | ----------------- |
| Crohn_colon             | 2        | 1                 |
| Crohn_blood             | 2        | 1                 |
| Spondyloarhtritis_colon | 2        | 0                 |
| Spondyloarhtritis_blood | 2        | 0                 |
| Co-morbid_colon         | 2        | 2                 |
| Co-morbid_blood         | 2        | 0                 |
| Health_colon            | 2        | 1                 |
| Health_blood            | 2        | 1                 |
| 总计                    | 16       | 6                 |

### GSE165087

因为所有样本都是是针对白血病的研究，所以省略

### GSE171555

COVID-19中人PBMC的scRNA-seq，有四个队列一共四十八个样本，TLB分布情况如下

| 分组              | 重复数量 | 有TLB的样本的数量 |
| ----------------- | -------- | ----------------- |
| PBMC_hospitalized | 6        | 1                 |
| PBMC_infected     | 29       | 21                |
| PBMC_exposed      | 8        | 6                 |
| PBMC_healthy      | 5        | 5                 |
| 总计              | 48       | 34                |

