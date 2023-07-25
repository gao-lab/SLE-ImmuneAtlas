# Code of single cell analysis of SLE samples

This repo will keep the reproducible analysis pipeline of the manuscript "Global immune repertoire profiling suggests multifaceted unresolved dysregulation upon systemic lupus erythematosus immunosuppressant therapy".


## Setup environment

### R env
We use `R-4.1.3`, and record the all R packages and their version in `renv.lock` file. To reproduce the environment, just run

```R
install.packages("renv")
renv::restore()
```

### Python env
We use `conda` to manage Python environment. To reproduce the environment, run

```bash
conda env create --name <ENV_NAME> --file=environments.yml
```

## Code structure

The code is organized as follows:

```bash
.
│── Figure/            # dir keeps the scripts for reproducing figures in the manuscript
│── other_sc_data/     # dir keeps the scripts for building the single cell SLE atlas
│── output_file/       # dir keeps the key outputs of the analysis pipeline
├── scripts/           # dir keeps the scripts for RNA analysis
│── vdj/               # dir keeps the scripts for BCR analysis
│── tcr_vdj/           # dir keeps the scripts for TCR analysis
│── renv.lock          # file keeps R packages and their version
├── environments.yml   # file keeps conda environment
```

## Data

Single cell RNA matrix and VDJ data will be available upon publication. 

## Reproduce the analysis

### RNA 
We mainly use R package `Seurat` to analyze single cell RNA data, and analysis pipeline is in directory  `scripts/seurat`. Other specific analysis is also in `scripts/`, such as `scripts/monocle2`, `scripts/CellChat` and `scripts/scvelo`.

### TCR and BCR
We use R package `immcatation`, `VDJtools` and `scRepertoire` to analyze the single cell TCR/BCR data. The analysis pipeline is in `tcr_vdj/` and `vdj` respectively.


### Reproduce plots in the manuscript

Run `Figure/Figures_plot.R` to reproduce main figures in manuscript. Input data needs to be requested from the authors. (we can not put it online due to the privacy of the patients)


## Publication

Our manuscript is under review.

## Contact

If you have any questions, please contact Chen-Rui Xia (xiachenrui@mail.cbi.pku.edu.cn) or open an issue on GitHub.

