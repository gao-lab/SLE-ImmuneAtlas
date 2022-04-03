import pandas as pd
from snakemake.utils import R
import scanpy as sc
import os
import numpy as np
import loompy as lp
# NOTE: for 20k cell, this flow will take about 8 hours 

SAMPLES = ['bcell']

rule all:
    input:
        loom = expand('{sample}/pySCENIC_out.loom', sample = SAMPLES)
        
rule seurat_to_h5ad:
    input:
        seu_obj = '{sample}/seurat.rdata'
    output:
        h5ad = '{sample}/convert.h5ad'
    script:
        'scripts/seurat_to_h5ad.R'

        
rule h5ad_to_loom:
    input:
        h5ad = '{sample}/convert.h5ad'
    output:
        loom = '{sample}/input.loom'
    run:
        adata = sc.read_h5ad(input.h5ad)
        # create basic row and column attributes for the loom file:
        row_attrs = {
            "Gene": np.array(adata.var_names) ,
        }
        col_attrs = {
            "CellID": np.array(adata.obs_names) ,
            "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
            "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
            'subtype': np.array(adata.obs['subtype'])
        }
        lp.create( output.loom, adata.X.transpose().todense(), row_attrs, col_attrs)

        
rule pySCENIC:
    input:
        loom = '{sample}/input.loom'
    output:
        loom = '{sample}/pySCENIC_out.loom'
    params: 
        # f_loom_path_scenic="/rd2/user/xiacr/sle/output_file/seurat/mono_dc/mono_dc_filter_harmony_sceasy_scanpy.loom"
        f_tfs="/rd2/user/xiacr/sle/source/SCENIC/resource/hs_hgnc_curated_tfs.txt",
        f_db_names="/rd2/user/xiacr/sle/source/SCENIC/resource/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
        f_motif_path="/rd2/user/xiacr/sle/source/SCENIC/resource/motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
        # f_pyscenic_output="/rd2/user/xiacr/sle/output_file/seurat/mono_dc/pyscenic_out_mono_dc_filter_harmony.loom"
        cores  = 36
    shell:
        """
        eval "$(conda shell.bash hook)"
        conda activate scenic_protocol
        
        # echo  --------------------------------step1 -------------------------------
        pyscenic grn -o adj.csv --num_workers {params.cores} --sparse {input.loom} {params.f_tfs} 

        # echo --------------------------------step2 -------------------------------
        pyscenic ctx adj.csv \
            {params.f_db_names} \
            --annotations_fname {params.f_motif_path} \
            --expression_mtx_fname {input.loom} \
            --output reg.csv \
            --mask_dropouts \
            --num_workers {params.cores}

        # echo -------------------------------step3 ----------------------------- 
        pyscenic aucell \
            {input.loom} \
            reg.csv \
            --output {output.loom} \
            --num_workers {params.cores}
    
        """
        