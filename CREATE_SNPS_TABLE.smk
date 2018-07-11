#####
### Ariane's pipeline
### Remember to include all the requird scripts
### Problems:
### 1. where and how to include the required scripts?
### 2. python2
### 3. requiring flatcounts

rule all_snps_table:
    input:
        flt_files=expand("{TMP_D}/{strains}/{mapper}.vcf.gz", 
            TMP_D=TMP_D, strains=STRAINS, mapper= ['bwa']),
        snp_vcf_files=expand("{TMP_D}/{strains}/{mapper}.snp.vcf", 
            TMP_D=TMP_D, strains=STRAINS, mapper= ['bwa']),
        dict_f='strain_list',
        anno_f='anno.tab'
    output: 'all_SNPs.v4.tab'
    shell:
        'source activate py27;'
        'python ./mutation_table.v4.py  -f {input[dict_f]} '
        '-a {input[anno_f]} -o {output};'
        'source deactivate'

rule syn_snps_table:
    input:
        all_snps_tab='all_SNPs.v4.tab',
        ref_gbk='reference.gbk'
    output: 'syn_SNPs_final.v4.tab' 
    shell:
        'source activate py27;'
        'python ./Snp2Amino.py -f {input[all_snps_tab]} -g {input[ref_gbk]} '
        '-n all -o {output} ;'
        'source deactivate'

rule nonsyn_snps_table:
    input:
        all_snps_tab='all_SNPs.v4.tab',
        ref_gbk='reference.gbk'
    output: 'non-syn_SNPs_final.v4.tab' 
    shell:
        'source activate py27;'
        'python ./Snp2Amino.py -f {input[all_snps_tab]} -g {input[ref_gbk]} '
        '-n non-syn -o {output} ;'
        'source deactivate'

rule filter_vcf:
    input:
        vcf_gz_f="{TMP_D}/{strains}/{mapper}.vcf.gz"
    output:
        snp_vcf_f="{TMP_D}/{strains}/{mapper}.snp.vcf"
    shell:
        'vcftools --gzvcf {input[vcf_gz_f]} --remove-indels --recode '
        '--recode-INFO-all --stdout > {output[snp_vcf_f]}'
        

rule compute_flatcounts:
    input:
       cov_f='{TMP_D}/{strain}/{mapper}.coverage'
    output:
       flt_f='{TMP_D}/{strain}/{mapper}.flatcount'
    shell:
       'python bam_cov2flatcount.py -in_f {input[cov_f]} -out_f {output[flt_f]}'
    
rule count_cov:
    input:
        bam_f='{TMP_D}/{strain}/{mapper}.sorted.bam'
    output:
        cov_f='{TMP_D}/{strain}/{mapper}.coverage'
    shell:
        'source activate Ariane_dna;'
        'samtools depth -a {input[bam_f]} > {output[cov_f]};'
        'source deactivate'
