import snakemake
import pandas as pd

# README 
# VDJ tools for TCR in groups
# NOTE: we not need step1 and step2 in this pipeline 
# NOTE: please run `Rscript pretreatment.R` in cmd before

# meta = pd.read_csv("./all_sample_meta.csv")
# SAMPLES = meta["name"]
SAMPLES = ["HC_TRA","HC_TRB",'SLE_TRA','SLE_TRB','HC_tcr','SLE_tcr']
# BEFORE = ['XH','HXR','WYF','LL','ZPP']
# AFTER = ['XYY','HXR2','WYF2','LL2','ZPP2']
OVERLAP_METHOD =['R', 'D', 'F' ,'F2', 'vjJSD' ,'vj2JSD','sJSD', 'Jaccard' ,'MorisitaHorn']
# low2up = {i:i.upper() for i in OVERLAP_METHOD}
# up2low = {i.upper():i for i in OVERLAP_METHOD}
# OVERLAP_METHOD =[ 'D' ,'F2', 'Jaccard']

# meta["#file.name"] = "./out/vdj_format/"+meta["name"]+"/tcr_vdjtools.tsv"
# meta["sample.id"] = meta["name"]
# meta.insert(0, "sample.id", meta.pop("sample.id"))
# meta.insert(0, "#file.name", meta.pop("#file.name"))
# meta.to_csv("./vdjtools_meta.tsv", index=False,sep = "\t")

rule all:
    input:
        # expand("./out/vdj_format/{sample}/tcr_vdjtools.tsv", sample = SAMPLES), # STEP1
        # "tcr_vdjtools.basicstats.txt",  # STEP2
        expand("{sample}_tcr.fancyvj.wt.pdf", sample = SAMPLES), # STEP3
        # expand("{before}_{after}.aa.paired.scatter.pdf", zip, before=BEFORE, after=AFTER), # STEP4
        # "pairwise_distances.intersect.batch.aa.pdf", # STEP5
        # "pairwise_distances.intersect.batch.aa.txt", # STEP5
        # expand("cluster_sample.hc.aa.{method}.pdf", method= OVERLAP_METHOD),# STEP6
        # expand("cluster_sample.mds.aa.{method}.pdf", method= OVERLAP_METHOD), # STEP6
        # expand("cluster_sample.perms.aa.{method}.pdf", method= OVERLAP_METHOD)
         
        
# STEP1: convert 10x cellranger format to vdjtools format(single chain mode) √
rule convert_10x_to_VDJtools:
    input:
        raw_file = "./tcr_data/{sample}_tcr.csv"
    output:
        tmp_dir = directory("./out/vdj_format/{sample}"),
        vdj_file = "./out/vdj_format/{sample}/tcr_vdjtools.tsv"
    script:
        "script/step1_convert_10x_to_VDJtools.R"
        
        
# STEP2:(need STEP1) vdj tools analysis basic: whole
rule vdjtools_basic_whole:
    input:
        # rules.convert_10x_to_VDJtools.output.vdj_file # not useful
        meta_path = "./vdjtools_meta.tsv",
        # vdj_file = "./out/vdj_format/{sample}/tcr_vdjtools.tsv"
    output:
        stat_file = "tcr_vdjtools.basicstats.txt",  #can not use params
        # plot_heatmap = "",
        # plot_circle = "{sample}.fancyvj.wt.pdf"
    params:
        prefix = "tcr_vdjtools"
    shell:
        # java -jar ../vdjtools/vdjtools-1.2.1.jar PlotFancyVJUsage  {input.vdj_file}  {wildcards.sample}_tcr
        """
        java -jar ../vdjtools/vdjtools-1.2.1.jar CalcBasicStats -m {input.meta_path} {params.prefix} 
        java -jar ../vdjtools/vdjtools-1.2.1.jar CalcSegmentUsage -p -m {input.meta_path} -f group -l name {params.prefix}
        """
        
# STEP3:(need STEP1) vdj tools analysis basic: single sample 
# LL2 MXY LAY failed for some reaeson
rule vdjtools_basic_single:
    input:
        vdj_file =  "{sample}_VDJtools.txt"
        # rules.convert_10x_to_VDJtools.output.vdj_file # not useful
        # meta_path = "./vdjtools_meta.tsv",
    output:
        # stat_file = "tcr_vdjtools.basicstats.txt",  #can not use params
        # plot_heatmap = "",
        plot_circle = "{sample}_tcr.fancyvj.wt.pdf"
    # params:
    #     prefix = "tcr_vdjtools"
    shell:
        # java -jar ../vdjtools/vdjtools-1.2.1.jar CalcBasicStats -m {input.meta_path} {params.prefix} 
        # java -jar ../vdjtools/vdjtools-1.2.1.jar CalcSegmentUsage -p -m {input.meta_path} -f group -l name {params.prefix} 
        """
        java -jar ../vdjtools/vdjtools-1.2.1.jar PlotFancyVJUsage  {input.vdj_file}  {wildcards.sample}_tcr
        """

# STEP4:(need STEP1) vdj tools compare paird samples 
rule compare_share:
    input:
        vdj_file_before =  "./out/vdj_format/{before}/tcr_vdjtools.tsv",
        vdj_file_after =  "./out/vdj_format/{after}/tcr_vdjtools.tsv",
    output:
        plot_pair = '{before}_{after}.aa.paired.scatter.pdf'
    shell:
        # before and after can not have same name(even in different dir)
        # rm {input.vdj_file_before}_tmp
        '''
         cp {input.vdj_file_before} ./out/vdj_format/{wildcards.before}/{wildcards.before}_tcr_vdjtools.tsv
         java -jar ../vdjtools/vdjtools-1.2.1.jar OverlapPair -p -i aa  ./out/vdj_format/{wildcards.before}/{wildcards.before}_tcr_vdjtools.tsv {input.vdj_file_after} {wildcards.before}_{wildcards.after}

         rm ./out/vdj_format/{wildcards.before}/{wildcards.before}_tcr_vdjtools.tsv
        '''
        
# STEP5:(need STEP1) vdj tools pairwise distances
rule pairwise_distances:
    input:
        meta_path = "./vdjtools_meta.tsv"
    output:
        pic = "pairwise_distances.intersect.batch.aa.pdf",
        file = "pairwise_distances.intersect.batch.aa.txt"
    params:
        prefix = "pairwise_distances"
    shell:
        '''
        java -jar ../vdjtools/vdjtools-1.2.1.jar CalcPairwiseDistances -p -m  {input.meta_path} -i aa {params.prefix} 
        '''

        
# STEP6:(need STEP5) vdj tools compare cluster 
rule cluster_sample:
    input:
        "pairwise_distances.intersect.batch.aa.txt"
    output:
        file1 = "cluster_sample.mds.aa.{method}.txt",
        # file2 = "cluster_sample.hc.aa.{method}.txt",
        pic1 = "cluster_sample.hc.aa.{method}.pdf",
        pic2 = "cluster_sample.mds.aa.{method}.pdf"
    params:
        prefix = "cluster_sample",
        # method = lambda wildcards: up2low[wildcards.method]
        # measure = 'F' # can be one of {R D F F2 vjJSD vj2JSD sJSD Jaccard MorisitaHorn} 
    shell:
        # for file in *$(typeset -l {method})*; do mv $file ${file//ABC/XYZ} ; done
        # for file in `ls|grep `echo {wildcards.method}| tr a-z A-Z` `; do mv $file ${file//`echo {wildcards.method}| tr a-z A-Z`/{wildcards.method}} ; done
        # bash rename.sh
        '''
        cp {input} cluster_sample.intersect.batch.aa.txt
        java -jar ../vdjtools/vdjtools-1.2.1.jar ClusterSamples -e {wildcards.method} -i aa -p  -f treatment -l name   {params.prefix} 
        upper=`echo {wildcards.method}| tr a-z A-Z`
        if [[ $upper !=  {wildcards.method} ]] ;then
            mv cluster_sample.mds.aa.$upper.txt cluster_sample.mds.aa.{wildcards.method}.txt
            mv cluster_sample.mds.aa.$upper.pdf cluster_sample.mds.aa.{wildcards.method}.pdf
            mv cluster_sample.hc.aa.$upper.pdf cluster_sample.hc.aa.{wildcards.method}.pdf
        fi
        '''
        
# STEP7:(need STEP6) vdj tools testclusters cluster 
rule test_cluster:
    input:
        # lambda wildcards: "cluster_sample.mds.aa.{}.txt".format(low2up[wildcards.method]),
        file = 'cluster_sample.mds.aa.{method}.txt'
        # "cluster_sample.hc.aa.{method}.txt"
    output:
        pic = "cluster_sample.perms.aa.{method}.pdf"
    params:
        prefix = "cluster_sample",
        # measure = 'F' # can be one of {R D F F2 vjJSD vj2JSD sJSD Jaccard MorisitaHorn} 
    shell:
        '''
        upper=`echo {wildcards.method}| tr a-z A-Z`
        if [[ $upper !=  {wildcards.method} ]] ;then
            mv cluster_sample.mds.aa.{wildcards.method}.txt cluster_sample.mds.aa.$upper.txt
            mv cluster_sample.mds.aa.{wildcards.method}.pdf cluster_sample.mds.aa.$upper.pdf
            mv cluster_sample.hc.aa.{wildcards.method}.pdf cluster_sample.hc.aa.$upper.pdf
        fi
        java -jar ../vdjtools/vdjtools-1.2.1.jar TestClusters  -e {wildcards.method} -i aa  {params.prefix}      
        if [[ $upper !=  {wildcards.method} ]] ;then
            mv cluster_sample.perms.aa.$upper.pdf cluster_sample.perms.aa.{wildcards.method}.pdf
        fi
        '''