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
        syn_output= os.path.join(TMP_D, 'syn_SNPs.tab')
    output: 
        config['nonsyn_snps_table'],
        config['syn_snps_table']
    params:
        strains= STRAINS
    script: 'lib/collect_snps_data.py'

rule syn_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output:
        syn_output= os.path.join(TMP_D, 'syn_SNPs.tab')
    params:
        script_f='snps_original_methods/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -n syn \
        -o {output.syn_output} 
        source deactivate
        """

rule nonsyn_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output: 
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab'),
    params:
        script_f='snps_original_methods/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -n non-syn \
        -o {output.non_syn_output} 
        source deactivate
        """

rule all_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output:
        all_snps_tab=os.path.join(TMP_D, 'all_SNPs.tsv')
    params:
        script_f='snps_original_methods/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -o {output.all_snps_tab} 
        source deactivate
        """

rule for_tab_copy_files:
    ## duplicate everyone to create 'STRAIN.flt.vcf' and 'STRAIN.flatcount'
    input: 
        flatcount_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, SOFTWARE['mapper'], 'st_variant.fltcnt'),
        snp_vcf_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, SOFTWARE['mapper'], 'st_variant.snp-vcf')
    output:
        flatcount=temp("{strain}.flatcount"),
        flt_vcf=temp("{strain}.flt.vcf")
    shell:
        """
        ln -s {input.snp_vcf_file} {output.flt_vcf}
        ln -s {input.flatcount_file} {output.flatcount}
        """

rule create_dict_file: 
    input: 
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=STRAINS),
        flatcount_files=expand("{strain}.flatcount", 
            strain=STRAINS)
    output:
        dict_file= 'dict.txt'
    params:
        strains=STRAINS
    shell:
        """
        echo {input.flt_vcf_files}| \
        sed 's/\.flt\.vcf\W*/\\n/g'| \
        grep '\w' > {output.dict_file}
        """

rule all_snps_list:
    input:
        ## the exact files should also be, so they cannot be deleted before 
        flatcount_file=expand(os.path.join(TMP_D, '{strain}', SOFTWARE['mapper'], 'st_variant.fltcnt'), strain=STRAINS),
        snp_vcf_file=expand(os.path.join(TMP_D, '{strain}', SOFTWARE['mapper'], 'st_variant.snp-vcf'), strain=STRAINS),
        flt_files=expand("{strain}.flatcount", 
            strain=STRAINS),
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=STRAINS),
        dict_f='dict.txt',
        anno_f='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
    output: 
#        all_snps_tab=os.path.join(TMP_D, '/all_SNPs.tsv'),
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab'))
    params:
        #script_f='lib/mutation_table.v5.py'
        script_f='snps_original_methods/mutation_table.py'
    shell:
        """
        source activate py27
        python {params.script_f}  -f {input.dict_f} \
        -a {input.anno_f} \
        -o {output.snps_list}
        source deactivate
        """


rule compute_flatcounts:
    input:
        STAMPY_SAM= '{TMP_D}/{strain}/{mapper}/st_paired.sam'
    output:
        flt_f='{TMP_D}/{strain}/{mapper}/st_variant.fltcnt',
        sin_f='{TMP_D}/{strain}/{mapper}/st_variant.sin',
        art_f='{TMP_D}/{strain}/{mapper}/st_variant.art'
    params:
        script_f="lib/sam2art.pl"
    shell:
        """
        {params.script_f} -s 2 -p -4 {input.STAMPY_SAM} > {output.art_f} 
        {params.script_f} -s 2 -l -p -4 {input.STAMPY_SAM} > {output.sin_f}
        {params.script_f} -f -s 2 -p {input.STAMPY_SAM} > {output.flt_f}
        """

rule count_cov:
    input:
        bam_f='{TMP_D}/{strain}/{mapper}/st_sorted.bam'
    output:
        cov_f='{TMP_D}/{strain}/{mapper}/st_mapping.coverage'
    shell:
        'samtools depth -a {input[bam_f]} > {output[cov_f]};'
'''
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
'''
