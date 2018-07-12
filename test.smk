import pandas as pd
import os
import re
configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
TMP_D=(config['tmp_d'] if re.search('\w', config['tmp_d']) else '.')
CORES=config['cores']
STAMPY_EXE=config['stampy_exe']
RAXML_EXE=config['raxml_exe']
RESULT_D=config['result_d']

rule rna_mapping:
    input:
        FQ=lambda wildcards: SAMPLES_DF.loc[wildcards.sample, 'rna_reads'],
        REF=REF_FA,
        REF_STAMPY_INDEX1=REF_FA+".stidx",
        REF_STAMPY_INDEX2=REF_FA+".sthash",
        REF_BWA_INDEX1=REF_FA+".bwt",
        REF_BWA_INDEX2=REF_FA+".fai"
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY=STAMPY_EXE,
    output:
        RMDUP_SAM=temp(TMP_D+'/{sample}/rna_sorted.rmdup.sam'),
        BAM=temp(TMP_D+'/{sample}/rna_sorted.bam'),
        RMDUP_BAM=temp(TMP_D+'/{sample}/rna_sorted.rmdup.bam')
    #threads: 16
    shell:
       # 'samtools view -bS {input} |samtools sort | samtools rmdup -s |samtools view -h > {output} '
        #"source activate Ariane_dna; source activate py27;" ## needs to be dealt with wrappers
        "source activate py27;" ## needs to be dealt with wrappers
        '{params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" '
        '-t{params[CORES]} -g {params[REF_PREFIX]} -h {params[REF_PREFIX]} '
        '-M {input[FQ]} | samtools view -bS -@ {params[CORES]} |'
        'samtools sort > {output[BAM]};'
        'samtools rmdup -s {output[BAM]} {output[RMDUP_BAM]};'
        'samtools view -h {output[RMDUP_BAM]} > {output[RMDUP_SAM]} ;'
        #"samtools index {output};"
        "source deactivate;"

rule convert_to_art:
    input: 
        '{TMP_D}/{sample}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{sample}/rna_sorted.rmdup.art')
    params:
        SAM2ART=config['sam2art_exe']
    shell:
        '{params[SAM2ART]} -s 2 {input} > {output} '

rule convert_to_sin:
    input: 
        '{TMP_D}/{sample}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{sample}/rna_sorted.rmdup.sin')
    params:
        SAM2ART=config['sam2art_exe']
    shell:
        '{params[SAM2ART]} -s 2 -l {input} > {output} '

rule convert_to_flatcount:
    input: 
        '{TMP_D}/{sample}/rna_sorted.rmdup.sam'
    output:
        temp('{TMP_D}/{sample}/rna_sorted.rmdup.flatcount')
    params:
        SAM2ART=config['sam2art_exe']
    shell:
        '{params[SAM2ART]} -f -s 2 {input} > {output} '
'''
rule make_annot_tab:
    input: 
        REF_GBK
    output:
        ANNO_F='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
'''
rule gene_counts:
    input:
        SIN_F='{TMP_D}/{sample}/rna_sorted.rmdup.sin'
    output:
        RPG_F='{TMP_D}/{sample}/rna_sorted.rmdup.rpg'
    params:
        ART2GENECOUNT=config['art2genecount_exe'], 
        ANNO_F='Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab'
    shell:
        '{params[ART2GENECOUNT]} -b -a {input[SIN_F]} -t tab '
        '-r {params[ANNO_F]} > {output[RPG_F]}'

rule write_dict_file:
    input:
        RPG_FILES=expand('{TMP_D}/{sample}/rna_sorted.rmdup.rpg', sample= STRAINS,
            TMP_D= TMP_D)
    output:
        #dict_f= TMP_D+'/rpg_dict'
        dict_f= 'tmp/rpg_dict'
    run:
        out_f= output.dict_f
        out_fh= open(out_f, 'w')
        out_fh.write('name\tpath\n')
        for n in range(len(input.RPG_FILES)):
            out_fh.write('{}\t{}\n'.format(STRAINS[n], input.RPG_FILES[n]))
        out_fh.close()

rule create_and_make_expr_table:
    input:
        dict_f= 'tmp/rpg_dict'
    output:
        expr_table='results/expr.tab'
    params:
        ANNO_F='Pseudomonas_aeruginosa_PA14_12genes_R_annotation'
    script: 'collect_rpg_data.R'
