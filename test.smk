import pandas as pd
import os
import re
configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
#TMP_D=(config['tmp_d'] if re.search('\w', config['tmp_d']) else '.')
TMP_D='seq2geno_temp'
CORES=config['cores']
STAMPY_EXE=config['stampy_exe']
RAXML_EXE=config['raxml_exe']
#RESULT_D=config['result_d']
RESULT_D='seq2geno'
#SOFTWARE={'mapper': 'bwa'}
SOFTWARE= config['software']
SOFTWARE['annotator']= 'prokka'
SOFTWARE['gene_sorter']= 'roary'
print(SAMPLES_DF)
print(SOFTWARE)

rule all: 
    input:
        aln= dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln'))

rule sort_gene_family :
    # compute family-wise alignments
    input:
        roary_gpa= lambda wildcards: os.path.join(TMP_D,
            SOFTWARE['gene_sorter'],'gene_presence_absence.csv'),
        gene_dna_seqs= lambda wildcards: [os.path.join(TMP_D, strain,
            SOFTWARE['assembler'], SOFTWARE['annotator'], 
            'de_novo.ffn') for strain in STRAINS]
    output:
        aln= dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln'))
    params:
        CORES= CORES,
        MIN_SAMPLES= int(len(STRAINS)/2),
        STRAINS= STRAINS,
        TMP_D= TMP_D+'/extracted_proteins_nt'
    script: 'lib/makeGroupAln.py'

rule aln_to_vcf_per_fam:
    input:
        FAM_ALN=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.aln')
    output:
        FAM_VCF=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.vcf') 
    params:
        cores= CORES,
        msa2vcf_script='lib/jvarkit/dist/msa2vcf.jar',
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3'
    shell:
        """
        java -jar {params.msa2vcf_script} \
        < {input[FAM_ALN]} > {output[FAM_VCF]} 
        """

rule vcf_to_indels_per_fam:
    input:
        FAM_VCF=os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}.vcf')       
    output:
        FAM_INDELS_TXT=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indels.txt'),
        FAM_INDELS_GFF=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indels.gff'),
        FAM_INDELS_STATS=os.path.join(TMP_D, 'extracted_proteins_nt',
'{fam}_indel_stats.txt'),
    params:
        cores= CORES,
        working_dir='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v3',
        strains_perc_cutoff= 0.5,
        len_cutoff= 8
    script: 'lib/indel_detection/vcf2indel.py'

rule expand_by_family:
    input:
        FAM_INDELS_FILES=dynamic(os.path.join(TMP_D, 'extracted_proteins_nt', '{fam}_indels.txt'))        
    output:
        indel_list= os.path.join(TMP_D, 'indel.list')
    run:
        out_fh=open(output.indel_list, 'w')
        out_fh.write("\n".join(input.FAM_INDELS_FILES))
        out_fh.close()

rule create_indel_table:
    input:
        indel_list= os.path.join(TMP_D, 'indel.list')
    output:
        annot_f=config['indel_table'],
        indel_stat_f='indels_stats.txt'
    script: 'lib/indel_detection/generate_indel_features.py'
