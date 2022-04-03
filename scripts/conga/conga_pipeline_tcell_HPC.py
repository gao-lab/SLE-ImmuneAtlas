import scanpy as sc
import pandas as pd
import numpy as np
import snakemake

# This pipeline can not run automaticly ,we should do these manually
# choose the input dir
# prepare the tcell_all_meta.csv

meta = pd.read_csv("./all_sample_meta.csv")
SAMPLES = meta["name"]


rule all:
    input:
        final ='/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step3/all_tcell_results_summary.html'
    threads: 1

        
# STEP1: Set up clone file(can not para)
rule set_up:
    input:
        contig = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/source/10x_tcr/{sample}_tcr.csv',
    output:
        clone_file = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step1/{sample}_out.tsv',
    params:
        organism = 'human',
    threads: 1,
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga_new_env

        python /gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/conga/scripts/setup_10x_for_conga.py \
        --filtered_contig_annotations_csvfile {input.contig} \
        --output_clones_file {output.clone_file} \
        --organism {params.organism} \
        # --condense_clonotypes_by_tcrdist --tcrdist_threshold_for_condensing 50 \
        --no_kpca 
        '''
        
# STEP2: Create meta(use a trick to copy file)
rule creat_meta:
    input:
        clone_file = expand('/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step1/{sample}_out.tsv',sample = SAMPLES),
        h5 = expand('/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/source/10x_rna/{sample}_rna.h5',sample = SAMPLES),
    output:
        meta = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/tcell_all_meta.tsv',
    params:
        file = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/tcell_all_meta.tsv'
    threads: 1,
    shell:
        '''
        cp  -f {params.file} {output.meta}
        '''
        
    
# STEP3: Merge samples        
rule merge:
    input:
        meta = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/tcell_all_meta.tsv',
    output:
        clone_merged = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step2/all_tcell_merged_clones.tsv',
        h5ad_merged = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step2/all_tcell_merged_gex.h5ad',    
    params:
        organism = 'human'
    threads: 32,
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga_new_env

        echo starting step2
        python /gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/conga/scripts/merge_samples.py \
        --samples {input.meta} \
        --output_clones_file  {output.clone_merged} \
        --output_gex_data {output.h5ad_merged} \
        --organism {params.organism} \
        --condense_clonotypes_by_tcrdist --tcrdist_threshold_for_condensing 50 
        
        '''

# STEP4: conga analysis
rule conga_analysis:
    input:
        clone_merged = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step2/all_tcell_merged_clones.tsv',
        h5ad_merged = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step2/all_tcell_merged_gex.h5ad',  
    output:
        report = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step3/all_tcell_results_summary.html'
    params:
        organism = 'human',
        dir = '/gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/TCR/step3/all_tcell'
    threads: 32,   
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga_new_env
        
        # STEP3: Analysis 
        python /gpfs2/gaog_pkuhpc/users/xiacr/xiehe/conga/conga/scripts/run_conga.py \
        --gex_data {input.h5ad_merged} \
        --gex_data_type h5ad \
        --clones_file {input.clone_merged} \
        --organism {params.organism} \
        --all \
        --outfile_prefix {params.dir}
        '''