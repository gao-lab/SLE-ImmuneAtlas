import scanpy as sc
import pandas as pd
import numpy as np
import snakemake

# This pipeline can not run automaticly ,we should do these manually
# choose the input dir
# prepare the bcell_all_meta.csv

meta = pd.read_csv("./all_sample_meta.csv")
SAMPLES = meta["name"]




rule all:
    input:
        final ='/rd2/user/xiacr/sle/scripts/conga/BCR/all/step3/all_bcell_results_summary.html'

        
# STEP1: Set up clone file(can not para)
rule set_up:
    input:
        contig = '/rd2/user/xiacr/sle/data/10x_bcr/{sample}_bcr.csv',
    output:
        clone_file = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step1/{sample}_out.tsv',
    params:
        organism = 'human_ig',
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga

        python /rd2/user/xiacr/sle/source/conga/scripts/setup_10x_for_conga.py \
        --filtered_contig_annotations_csvfile {input.contig} \
        --output_clones_file {output.clone_file} \
        --organism {params.organism} \
        --condense_clonotypes_by_tcrdist --tcrdist_threshold_for_condensing 50 \
        --no_kpca 
        '''
        
# STEP2: Create meta(use a trick to copy file)
rule creat_meta:
    input:
        clone_file = expand('/rd2/user/xiacr/sle/scripts/conga/BCR/all/step1/{sample}_out.tsv',sample = SAMPLES),
        h5 = expand('/rd2/user/xiacr/sle/data/10x_rna/{sample}_rna.h5',sample = SAMPLES),
    output:
        meta = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/bcell_all_meta.tsv',
    params:
        file = '/rd2/user/xiacr/sle/scripts/conga/bcell_all_meta.tsv'
    shell:
        '''
        cp  {params.file} {output.meta}
        '''
        
    
# STEP3: Merge samples        
rule merge:
    input:
        meta = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/bcell_all_meta.tsv',
    output:
        clone_merged = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step2/all_bcell_merged_clones.tsv',
        h5ad_merged = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step2/all_bcell_merged_gex.h5ad',    
    params:
        organism = 'human_ig'
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga

        echo starting step2
        python /rd2/user/xiacr/sle/source/conga/scripts/merge_samples.py \
        --samples {input.meta} \
        --output_clones_file  {output.clone_merged} \
        --output_gex_data {output.h5ad_merged} \
        --organism {params.organism} \
        --condense_clonotypes_by_tcrdist --tcrdist_threshold_for_condensing 50 
        
        '''

# STEP4: conga analysis
rule conga_analysis:
    input:
        clone_merged = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step2/all_bcell_merged_clones.tsv',
        h5ad_merged = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step2/all_bcell_merged_gex.h5ad',  
    output:
        report = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step3/all_bcell_results_summary.html'
    params:
        organism = 'human_ig',
        dir = '/rd2/user/xiacr/sle/scripts/conga/BCR/all/step3/all_bcell'
        
    shell:
        '''
        eval "$(conda shell.bash hook)"
        conda activate conga
        
        # STEP3: Analysis 
        python /rd2/user/xiacr/sle/source/conga/scripts/run_conga.py \
        --gex_data {input.h5ad_merged} \
        --gex_data_type h5ad \
        --clones_file {input.clone_merged} \
        --organism {params.organism} \
        --all \
        --outfile_prefix {params.dir}
        '''