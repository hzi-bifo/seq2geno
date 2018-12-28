import pandas as pd
import pprint

TMP_D='/net/metagenomics/data/from_moni/old.tzuhao/\
seq2geno/dev_versions/v9/seq2geno_temp'
#DNA_READS_F='/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/data/samples.3.dna_subset.tsv'
DNA_READS= {'ESP088':['/net/metagenomics/data/from_moni/old.tzuhao/\
seq2geno/data/strains/dna/subsampled/ESP088_1.fq.gz',
'/net/metagenomics/data/from_moni/old.tzuhao/\
seq2geno/data/strains/dna/subsampled/ESP088_2.fq.gz']}

#for l in open(DNA_READS_F):
#    strain,files= l.strip().split('\t')
#    files_d= files.split(',')
#    DNA_READS[strain]= files_d

rule all:
    input:
        spades_out='{}/{}/spades/scaffolds.fasta'.format(TMP_D, 'ESP088')
        
rule spades_create_assembly:
    input:
        READS=lambda wildcards: DNA_READS[wildcards.strain]
    output:
        spades_out="{TMP_D}/{strain}/spades/scaffolds.fasta"
    params:
        spades_outdir=lambda wildcards: os.path.join(wildcards.TMP_D, 
            wildcards.strain, 'spades'),
        SPADES_OPT='--careful',
        SPADES_BIN='spades.py'
    threads: 20
#    threads: 1
#    conda: 
    run:
        import subprocess
        import sys
        out_dir= params.spades_outdir
        cmd= []
        try:
            if len(input.READS) == 2:
                cmd= [params.SPADES_BIN, params.SPADES_OPT, '--threads',
                    str(threads), '-o', out_dir, 
                    '-1', input.READS[0], 
                    '-2', input.READS[1]]
            elif len(input.READS) == 1:
                cmd= [params.SPADES_BIN, params.SPADES_OPT, '--threads',
                    str(threads), '-o', out_dir, 
                    '-s', input.READS[0]]
            else:
                raise Exception('The reads are neither single nor paired')
            print(cmd)
            subprocess.run(cmd)
        except Exception as e:
            print(e)
            sys.exit()
