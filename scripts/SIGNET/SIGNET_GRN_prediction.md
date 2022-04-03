# SIGNET: GRN Prediction 

*Dec 31, 2021 with SIGNET version 1.0.0*

**Abstract**

SIGNET (Single-cell RNA-seq-based Gene Regulatory Network Prediction using Multiple-Layer Perceptrons Bagging) is a deep learning -based framework for genes regulatory network (GRN) prediction relying only on single-cell RNA-seq data. This tutorial offers a brief description on the whole workflow, preparation and prediction of GRN using SIGNET.

**Package**

* Python(3.6.13):  PyTorch, Scanpy
* version of packages:
      scanpy         --1.7.2
      scikit-learn   --0.23.2
      torch          --1.5.0+cpu
      numpy          --1.19.5              
      pandas         --1.1.5
      torchvision    --0.6.0+cpu

* R:  RcisTarget, AUCell



## Contents

[toc]



## 1 Introduction to SIGNET

SIGNET is a new method for gene regulatory network (GRN) reconstruction and prediction using only single-cell RNA-seq data. It  is mainly based on a deep learning framework, multiple-layers perception (MLP), for relationship prediction between transcription factors (TF) and the non-transcription factors (NTF), combining RcisTarget for motif analysis and AUCell for GRN activity scoring. The predicted network graph and AUC score matrix  obtained by SIGNET can be used further for multiple downstream analysis. For more information on SIGNET, we recommend the user to check the following article:

> SIGNET: single-cell RNA-seq-based gene regulatory network prediction using multiple-layer perceptron bagging (https://doi.org/10.1093/bib/bbab547)

Please cite this article if you use SIGNET in your research. 

## 2 Requirements

### 2.1 TF list

Since SIGNET uses the expression of  TFs to predict the expression of NTFs, a gene list for TFs is necessary. There are multiply sources for the TF list and we recommend the TF list provided on the GitHub of pySCENIC (https://github.com/aertslab/pySCENIC/tree/master/resources).  There are 1,789 and 1,721 transcription factors recorded for human (hs_hgnc_tfs.txt) and mouse (mm_mgi_tfs.txt) respectively in each of the txt file. 

> Note: The TF list is just used for determining the argument for the input of the MLP model. If the users would like to predict the relationships between TF and NTF in other species, they can appoint any other credited gene lists for TFs. 

### 2.2 Packages installation 

In order to speed up the whole method, SIGNET should be conducted on both Python and R platform during different steps. Here we present the installation of the core packages used in SIGNET.

For Python platform, the user should install `PyTorch` framework for MLP training and `Scanpy` for data analysis. 

For R platform, the core dependent packages are `RcisTarget` and `AUCell`. 

``` R
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()

# Install RcisTarget & AUCell
BiocManager::install(c("AUCell", "RcisTarget"))
```

### 2.3 Input: expression matrix

The input of SIGNET is the **expression matrix** of scRNA-seq:

* Each column represents a cell sample and each row represents a gene. 
* The row name of the matrix should be the gene-symbol of gene ID.
* Expression units: The preferred expression values are raw values. Since we will use the binarized data matrix for MLP training, the expression units (raw, TPM, FPKM/RPKM) have almost no effect on the binarization matrix.

## 3 SIGNET workflow

This tutorial goes through the steps in the SIGNET workflow by following: 

1. Data preparation: 
   * Screen for expression feature genes;
   * Data transformation by binarization;
2. Predicting the relationship between TFs and NTFs using MLP model:
   * Split the whole dataset into training set and validation set with a ratio of 7:3;
   * Bootstrap to adjust the proportion of positive cases not less than 1/6;
   * MLP training;
   * Use case-deletion method to obtain accuracy curve for each target NTF gene;
   * Use generalized extreme Studentized deviate test (gESD test) to screen for outliers;
   * Obtain predicted co-expressed TF-NTF pairs (i.e. a gene list).

3. Pruning using RcisTarget:
   * Screen for NTFs with high motif enrichment scores for each TF;
   * Obtain final TF-NTF regulatory pairs (i.e. regulons).
4. GRN activity scoring using AUCell:
   * Score regulons in each cell sample and earning AUC score matrix.
5. Downstream analysis demo: network analysis using regulons list,  dimension reduction and clustering using AUC score matrix, heatmap...

### 3.1 Data preparation & TF-NTF relationships prediction

In this part, SIGNET uses 

* **scRNA-seq expression matrix**, 

* **TF list** 

as input and exports the

* **predicted TF-NTF co-pairs list**. 

#### 3.1.1 MLP training and accuracy curve

First, conduct the `SIGNET.py` file. The user should give the pathway where the scRNA-seq data matrix stored and the pathway that the output files will be stored. Assume a fake pathway: ./home/test and all files needed are stored in the `test` file folder (including your TF files).  

``` 
## On Linux server 
$ python -u SIGNET.py --data_file /home/test/data.csv --output_path /home/test
```

`SIGNET.py` uses scRNA-seq expression matrix and TF list as input and exports 6 files:

* `co_tf_fc.txt` and `co_fc.txt`: The accuracy curve, which is readay for gESD test using R.
* `gene_tf.csv` and `gene_ntf.csv`: Two gene lists contain the tf and ntf used in MLP training.
* `data_tf_binary.txt` and `data_ntf_binary.txt`: The binarized scRNA-seq data matrix after qc, feature genes selection.

>Note: As for why SIGNET output the accuracy curve in two files, that's because when a tf gene treated as the prediction target, it is also a variable as the input of the neural network and this will cause overfitting and in order to get over this difficulty, whenever a tf gene treated as the prediction target, the input variable will not include itself (but to substitute  the true binarized expression with a zero vector). Thus `co_tf_fc.txt` stores the accuracy curve targeting tf genes and `co_fc.txt` stores the accuracy curve targeting ntf genes.

There are several parameters that could be modified personally  listed as following: 

* `data_file`: The input scRNA-seq gene expression file.
* `tf_list_file`: The input transcription factor genes file used for prediction.
* `species`: The species for SIGNET. Default is mouse.  

* `n_genes`: Number of feature genes for training SIGNET. Default is 5,000. We recommend that it can change from 5,000 to 7,000. Too much feature genes may slow down the training process significantly unless the user is accessible to GPU accelerating  or parallel structure.
* `batch_size`: The batch size used in the training process.  Default is 64.
* `n_epochs`: Number of Epochs for training SIGNET. Default is 30.
* `lr`: Learning rate used for SGD. Default is 10^-3. 
* `output_path`: The output file path.

>**Advance parameters:** 
>
>For these parameters, the users cannot assign the specific value unless they modify the primary `SIGNET.py` file. Since these parameters may have an influence on the model training and  performance, under some circumstances, users may change these values to earn better results. However, users should be extremely careful when modifying these parameters.
>
>* `times`:  Number of initialization times for single target gene training. Default is 10. This is to avoid the the accuracy for predicting target gene be lower than 50%, which means the training is failed. Since this is also a hyperparameter, users should turn it carefully or avoid turning it. 

#### 3.1.2 gESD test and co-expressed pairs

For the afterward procedure, we will turn to `R` platform. First load packages needed and extra functions in `SIGNET.R`. 

```R
# Load packages and annotation data
library(reshape2)
library(umap)
library(pheatmap)
library(igraph)
library(GGally) # ggplot2 >= 3.3.0
library(ggplot2)
windowsFonts(TNM = windowsFont("Times New Roman"))
windowsFonts(HNR = windowsFont("Helvetica Neue CE 55 Roman")) # Use specific font on Windows platform
library(RcisTarget)
library(AUCell)
# Load extra functions
source("SIGNET.R")
```

In order to obtain the final predicted co-expressed TF-NTF pairs, users need copy the `co_fc.txt`, `co_tf_fc.txt`, `gene_tf.csv` and `gene_ntf.csv` to the work catalogue. 

```R
# Load RcisTarget package if not
library(RcisTarget)
# Load gene sets to analyze, i.e the 'genesets' obtained in last step
# Select motif database to use (annotations)
data(motifAnnotations_mgi)
# Import the motif databases for RcisTarget
motifRankings <- importRankings("mm9-tss-centered-10kb-7species.mc9nr.feather")
gene <- colnames(motifRankings)
# Load data
ac_ntf <- read.table("co_fc.txt")
ac_tf <- read.table("co_tf_fc.txt")
gene_ntf <- read.csv("gene_ntf.csv")
gene_tf <- read.csv("gene_tf.csv")
# Merge data
coexpressed <- Merge(ac_ntf,ac_tf,gene_ntf,gene_tf,gene)
# gESD test and outliers screening
O <- Screen(coexpressed)
# Transform to gene list
copaired <- Copaired(O)
genesets <- Genesets(copaired)
genesets_list <- levels(as.factor(copaired$V1))
```

Thus we obtain the co-expressed TF-NTF pairs.

### 3.2 Pruning using RcisTarget

For RcisTarget, users need to provide the gene list and two databases: 

1. Gene-motif rankings: the rankings of all the genes for each motif;
2. The annotation of motifs to transcription factors.

And then the RcisTarget can narrow the target genes down for each transcription factor and thus pruning. 

#### 3.2.1 Databases download

##### 3.2.1.1 Gene-motif rankings

For gene-motif rankings, users can download from https://resources.aertslab.org/cistarget/ or use the codes below:

```R
featherURL <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather" 
download.file(featherURL, destfile = basename(featherURL))
```

>**Note:** the codes above just show a human example. Users may use different files in the databases and choose various motif-rankings according to  version of database(mc9nr, etc.), species, search space and number of orthologous species taken into account to score the motifs. For further information, see https://resources.aertslab.org/cistarget/  for help.

##### 3.2.1.2 Annotation of motifs to transcription factors

In order to connect the motifs with potentially combined transcription factors, users should provide the annotations. For database of `mc9nr` version, users can simply load the file from the `RcisTarget` package:

```R
# mouse:
data(motifAnnotations_mgi)
# human:
data(motifAnnotations_hgnc)
```

For other annotations, users can use `importAnnotations` function and import the files from the source of their need. 

For simplicity, we only show a mouse example below. For other species, the codes are similar. 

#### 3.2.2 Using RcisTarget

```R
# Motif enrichment analysis
Regulons <- Prune(O, genesets, genesets_list, motifRankings, motifAnnotations_mgi)
copaired2 <- Copaired2(Regulons)
genesets2 <- Genesets(copaired2)
genesets_list2 <- levels(as.factor(copaired2$V1))
# Save the predicted regulons
write.csv(genesets2, "genesets2.csv")
```

Thus we obtain the final predicted regulons in the `genesets2` file. Pay attention to that `genesets` is the gene list predicted by the MLP model and `genesets2` is what we want. 

### 3.3 GRN activity scoring with AUCell

In order to evaluate the activity of regulons in each cell sample, we need the primary scRNA-seq data matrix along with the annotations of cell types. 

```R
# Load AUCell package if not
library(AUCell)
# Load scRNA-seq data
data <- read.csv("data.csv", header = T, row.names = 1) # Users may use any method to load the data matrix whose entries are exactly data counts
meta <- read.csv("meta.csv", header = T, row.names = 1) # A vector whose length equals to the number of samples in data
# Obtain AUC score matrix
set.seed(123)
par(mfrow=c(1,1))
data_auc <- AUCellMatrix(data, genesets2, genesets_list2)
```

Thus we obtain the AUC score data matrix `data_auc`. 

### 3.4 Downstream analysis demo

We mainly show three downstream analysis in this tutorial, i.e. dimension reduction and scatter plot, heatmap and network visualization. Dimension reduction and heatmap are both based on the AUCell score matrix and the network visualization is based on the `copaired2` file. In order to draw similar figure in the paper, users should load another R file which contains the auxiliary functions used for visualization.

```R
source("SIGNET_plot.R")
```

For any one of the plots listed below, we just give the simple version of using, and if the users are not satisfied with the characteristics with the figures, users can change them in the way whatever they want.

#### 3.4.1 Dimension reduction and scatter plot

``` R
# Dimension reduction
dataumap <- umap::umap(as.data.frame(t(data)))
dataplot <- data.frame(dataumap$layout[,1:2], CellType = meta)
colnames(dataplot) <- c("UMAP_1", "UMAP_2","CellType")
# Scatter plot
plot.scatter(dataplot)
```

#### 3.4.2 Heatmap

```R
plot.heatmap(data_auc,meta,color)
```

`color` is a color vector, with the same length as the levels of meta.

#### 3.4.3 Network visualization

```R
# load data: copaired2
# Get network object: mode = "crosshighlight" - under this mode the ntf nodes with degrees over one will be highlight
start = c("FOS","FOSL2","IRF4","MYC","JUN","KLF4")
edgecol = c("#0073c2","#ff7f0e","#2ca02c","#d62728","#9467bd","#efc000")
net = getNetElement(copaired2, center = start, color = edgecol, mode = "crosshighlight")
# personlize node characteristics
mnodecol = c("#868686", "#003c67", "#a73030") # length of mnodecol equals to net$cross_max, in this case is 3
nodecol = c(edgecol, mnodecol) 
names(nodecol) = c(start, "NTF", paste0("MNTF", 2:net$cross_max))
nodesize = c(rep(30, length(start)), 10, rep(20, net$cross_max - 1))
names(nodesize) = c(start, "NTF", paste0("MNTF", 2:net$cross_max))

# network visualization
set.seed(102400)
plot.net(net$net)
```



## Reference

[1] Luo, Qinhuan ,  Y. Yu , and  X. Lan . "SIGNET: single-cell RNA-seq-based gene regulatory network prediction using multiple-layer perceptron bagging." Briefings in Bioinformatics (2021).

[2] https://resources.aertslab.org/cistarget/

[3] https://bioconductor.org/packages/devel/bioc/vignettes/RcisTarget/inst/doc/RcisTarget.html
