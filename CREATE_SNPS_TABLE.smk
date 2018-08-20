#####
### Ariane's pipeline
### Remember to include all the requird scripts
### Problems:
### 1. where and how to include the required scripts?
### 2. python2
### 3. requiring flatcounts

rule covert_snps_table:
    input:
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab'),
        syn_output= os.path.join(TMP_D, '/syn_SNPs.tab')
    output: 
        config['nonsyn_snps_table'],
        config['syn_snps_table']
    params:
        strains= STRAINS
    script: 'lib/collect_snps_data.py'

rule syn_snps_table:
    input:
        all_snps_tab=os.path.join(TMP_D, '/all_SNPs.tsv'),
        ref_gbk=config['reference_annotation']
#    output: config['syn_snps_table']
    output:
        syn_output= os.path.join(TMP_D, '/syn_SNPs.tab')
        
    params:
        script_f='lib/Snp2Amino.v2.py'
    shell:
        """
        source activate py27
        python {params[script_f]} -f {input[all_snps_tab]} -g {input[ref_gbk]} \
        -n all -o {output}
        source deactivate
        """

rule nonsyn_snps_table:
    input:
        all_snps_tab=os.path.join(TMP_D, '/all_SNPs.tsv'),
        ref_gbk=config['reference_annotation']
#    output: config['nonsyn_snps_table']
    output: 
        non_syn_output= os.path.join(TMP_D, '/non-syn_SNPs.tab')
    params:
        script_f='lib/Snp2Amino.v2.py'
    shell:
        """
        source activate py27
        python {params[script_f]} -f {input[all_snps_tab]} -g {input[ref_gbk]} \
        -n non-syn -o {output} 
        source deactivate
        """

rule all_snps_table:
    input:
        flt_files=expand('{TMP_D}/{strain}/{mapper}/st_variant.flatcount', 
            TMP_D=TMP_D, strain=STRAINS, mapper= SOFTWARE['mapper']),
        snp_vcf_files=expand("{TMP_D}/{strain}/{mapper}/st_variant.snp-vcf", 
            TMP_D=TMP_D, strain=STRAINS, mapper= SOFTWARE['mapper']),
        dict_f='strain_list',
        anno_f='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
    output: 
        all_snps_tab=os.path.join(TMP_D, '/all_SNPs.tsv'),
    params:
        script_f='lib/mutation_table.v5.py'
    shell:
        """
        source activate py27
        python {params[script_f]}  -f {input[dict_f]} \
        -a {input[anno_f]} -o {output[all_snps_tab]}
        source deactivate
        """

rule compute_flatcounts:
    input:
        cov_f='{TMP_D}/{strain}/{mapper}/st_mapping.coverage'
    output:
        flt_f='{TMP_D}/{strain}/{mapper}/st_variant.flatcount'
    params:
        script_f="lib/bam_cov2flatcount.py" 
    shell:
       'python {params[script_f]} -in_f {input[cov_f]} -out_f {output[flt_f]}'

rule count_cov:
    input:
        bam_f='{TMP_D}/{strain}/{mapper}/st_sorted.bam'
    output:
        cov_f='{TMP_D}/{strain}/{mapper}/st_mapping.coverage'
    shell:
        'samtools depth -a {input[bam_f]} > {output[cov_f]};'

rule filter_vcf:
    input:
        vcf_gz="{TMP_D}/{strain}/{mapper}/st_vcf.gz",
        vcf_gz_index= "{TMP_D}/{strain}/{mapper}/st_vcf.gz.tbi"
    output:
        snp_vcf_f="{TMP_D}/{strain}/{mapper}/st_variant.snp-vcf"
    shell:
        """
        vcftools --gzvcf {input[vcf_gz]} --remove-indels --recode \
        --recode-INFO-all --stdout > {output[snp_vcf_f]}
        """
