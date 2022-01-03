# smart-seq2数据处理和分析

### 原始数据

数据源自丁阳师兄的大规模筛选和tisch数据库中记录的smart-seq2实验（二者大部分是重合的）

| gse       | gsm.O1.pasted | final.single.file                                           | info | raw.all_cell | raw.TLB | raw.TLB.relaxed | raw.TLB.more.relaxed | organism_ch1.pasted | source_name_ch1.pasted            | smartseq | matrix    | all_cell | B cell | TLB  | clustering | 之前的结果 | note                                                         |
| --------- | ------------- | ----------------------------------------------------------- | ---- | ------------ | ------- | --------------- | -------------------- | ------------------- | --------------------------------- | -------- | --------- | -------- | ------ | ---- | ---------- | ---------- | ------------------------------------------------------------ |
| GSE146771 | Homo sapiens  | GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz               | gsxc | 10468        | 174     | 617             | 617                  | Homo sapiens        | CD45+ and CD45- cells;CD45+ cells | TRUE     | TPM       | 10468    | 978    | 19   | √          | none       | TLB是plamsa表型，只考虑B cell聚类中的TLB                     |
| GSE115978 | Homo sapiens  | GSE115978_tpm.csv.gz                                        | gsxc | 7186         | 164     | 256             | 286                  | Homo sapiens        | melanoma tumor                    | TRUE     | TPM       | 7186     | 943    | 144  | √          | 84         |                                                              |
| GSE120575 | Homo sapiens  | GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz | gsxc | 16292        | 104     | 413             | 413                  | Homo sapiens        | Melanoma single cell              | TRUE     | TPM       | 16292    | 1782   | 0    | ×          | 21         |                                                              |
| GSE72056  | Homo sapiens  | GSE72056_melanoma_single_cell_revised_v2.txt.gz             | gsxc | 4645         | 66      | 190             | 190                  | Homo sapiens        | human melanoma;melanoma tumor     | TRUE     | TPM       | 4645     | 602    | 83   | √          | 83         |                                                              |
| GSE140228 | Homo sapiens  | GSE140228_read_counts_Smartseq2.csv.gz                      | gsxc | 7074         | 109     | 151             | 560                  | Homo sapiens        | CD45+ cells                       | TRUE     | raw_count | 7074     | 1132   | 0    | ×          | none       | how to convert to TPM->tisch include raw data(finally use TISCH pre-processed data) |
| GSE103322 | Homo sapiens  | GSE103322_HNSCC_all_data.txt.gz                             | gsxc |              |         |                 |                      | Homo sapiens        | HNSCC tumor                       | TRUE     | TPM       | 5092     | 205    | 32   | ？         | 没看       | 多次聚类后才出现（不表达CD79A,但是表达MS4A1），singleR注释更偏向CLP |

数据与丁阳师兄筛查的结果的出入原因

1. 比丁阳的多：基于聚类，能够部分克服一些细胞drop out事件的影响
2. 比丁阳的少：只对B cell进行分析，其实在属于T cell的cluster中也有部分的细胞同时表达CD79A 和CD3D，这类细胞没有进入我的分析流程

### CCA alignment



<img src="https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726101939147.png" alt="image-20210726101939147" style="zoom:50%;" /><img src="https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726101958925.png" alt="image-20210726101958925" style="zoom:50%;" />

一些marker gene的表达情况

<img src="https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726102113112.png" alt="image-20210726102113112" style="zoom:67%;" />

![image-20210726102211984](https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726102211984.png)









使用singleR进行注释，大部分注释的结果是在B cell，少数是CLP和Monocytes，注释到T cell中的细胞非常的少

![image-20210726102525332](https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726102525332.png)

![image-20210726102655858](https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726102655858.png)

tlb细胞的注释

![image-20210726103924135](https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726103924135.png)

### check doublets via ‘solo’

solo是一个基于半监督深度学习的的doublet检测工具，在基准测试中胜出了之前的工具DoubletFinder和Scrublet，但是只有python版本，在我们全部的6个smart-seq2的数据集的所有细胞中运行了solo，只在GSE120575中检测到2个doublet（一共有16289个细胞），之后会用10X平台的数据作为对照

![image-20210726191453174](https://gitee.com/huhansan666666/picture/raw/master/img/image-20210726191453174.png)

