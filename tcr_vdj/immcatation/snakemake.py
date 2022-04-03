import snakemake
import os
import pandas as pd

# README for TCR 
# To run this pipeline: you should prepare
# 1. (software): igblast, R package in the script, add change-o python scripts into path
# 2. filtered_contig fasta and csv file of cell ranger output in a dir named by sample such as
#    ---GW
#      \___filtered_contig.fasta
#      \___filtered_contig.csv
# output: the output will collect in the sample dir 

meta = pd.read_csv("./all_sample_meta.csv")
# SAMPLES = meta["name"]
SAMPLES = ['GW'] # test only

rule all:
    input:
        file = expand('data/{sample}/10X_clone-pass_germ-pass.tsv', sample = SAMPLES),
    
rule cellrangerToAIRR:
    input:
        fq = 'data/{sample}/filtered_contig.fasta',
        contig = 'data/{sample}/filtered_contig_annotations.csv',
    output:
        file = 'data/{sample}/filtered_contig_igblast_db-pass.tsv',
        tmp = temp('data/{sample}/filtered_contig_igblast.fmt7'),
    params:
        igblast_path = '/rd2/user/xiacr/sle/vdj/immcatation/source/igblast',
        organism = 'human',
        loci = 'tr',
        IMGT_Human_source = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj',
    # threads = min(5, cores)
    shell:
        '''
        AssignGenes.py igblast -s {input.fq} -b {params.igblast_path} \
            --organism {params.organism} --loci {params.loci} --format blast --outdir data/{wildcards.sample}
        MakeDb.py igblast -i {output.tmp} -s {input.fq} \
            -r {params.IMGT_Human_source} --10x {input.contig} --extended --outdir data/{wildcards.sample}
        '''

# TRA and TRB file
rule split_heavy_light:
    input:
        passed = 'data/{sample}/filtered_contig_igblast_db-pass.tsv',
    output:
        light = 'data/{sample}/light_parse-select.tsv',
        heavy = 'data/{sample}/heavy_parse-select.tsv',
    shell:
        '''
        ParseDb.py select -d {input.passed} -f locus -u "TRA" \
            --logic all --regex --outname heavy --outdir data/{wildcards.sample}
        ParseDb.py select -d {input.passed} -f locus -u "TRB" \
            --logic all --regex --outname light  --outdir data/{wildcards.sample}
        '''
        
# get the clone distance of the repo
rule heavy_chain_clone_cluster:
    input:
        heavy = 'data/{sample}/heavy_parse-select.tsv' ,
    output:
        distance = 'data/{sample}/shazam_clone_distance.txt',
        distribution_pdf = 'data/{sample}/shazam_clone_distance.pdf',
    script:
        './script/distance.R'
        
# NOTE: use TRA or TRB as heavy chain??
rule heavy_chain_define_clone:
    input:
        heavy = 'data/{sample}/heavy_parse-select.tsv',
        distance = 'data/{sample}/shazam_clone_distance.txt',
    output:
        heavy_cluster = 'data/{sample}/heavy_parse-select_clone-pass.tsv',
    shell:
        '''
        DefineClones.py -d {input.heavy} --act set --model ham \
            --norm len --dist `cat {input.distance}`
        '''
        
rule light_chain_clone_assign:
    input:
        light = 'data/{sample}/light_parse-select.tsv',
        heavy_cluster = 'data/{sample}/heavy_parse-select_clone-pass.tsv',
    output:
        passed_clone = 'data/{sample}/10X_clone-pass.tsv',
    shell:
        '''
        light_cluster.py -d {input.heavy_cluster} -e {input.light} -o {output.passed_clone}
        '''
    

rule  germline_sequences_to_database:
    input:
        passed_clone = 'data/{sample}/10X_clone-pass.tsv',
    output:
        germ = 'data/{sample}/10X_clone-pass_germ-pass.tsv',
    params:
        TRAV = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj/imgt_human_TRAV.fasta',
        TRBV = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj/imgt_human_TRBV.fasta',
        TRAJ = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj/imgt_human_TRAJ.fasta',
        TRBJ = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj/imgt_human_TRBJ.fasta',
        TRBD = '/rd2/user/xiacr/sle/vdj/immcatation/source/germlines/imgt/human/vdj/imgt_human_TRBD.fasta',
    shell:
        '''V
        CreateGermlines.py -d {input.passed_clone} -g dmask --cloned \
            -r {params.TRAV} {params.TRBV} {params.TRAJ} {params.TRBJ} {params.TRBD}
        '''


