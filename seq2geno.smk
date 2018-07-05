import pandas as pd
import os

configfile: "config.yaml"
SAMPLES_DF=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
STRAINS=SAMPLES_DF['strain'].tolist()
REF_FA=config['reference_sequence']
REF_GBK=config['reference_annotation']
TMP_D=config['tmp_d']
CORES=config['cores']

rule all:
    input:
        #TMP_D+'/whole_genomes.aln' 
        #expand(TMP_D+"/{strain}/target_regions.fa", strain= STRAINS)
        #TMP_D+'/OneBig.aln' 
        #TMP_D+'/roary/gene_presence_absence.csv'
        'gpa.mat'
   
rule prepare_reference:
    input:
        REF_FA
    output:
        REF_FA+".fai",
        REF_FA+".bwt"

    shell:
        "samtools faidx {input};" 
        "bwa index {input}"

rule prepare_target_regions:
    input:
        REF_GBK
    output:
        TMP_D+"/reference.target_regions.fa"
    run:
        from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Alphabet import IUPAC
        import re
	# read the gbk
        gbk_f= input[0]
	rec= SeqIO.read(gbk_f, 'gb')
	chromosome= rec.id
	# filter target regions types (gene, CDS...etc)
	features= [fea for fea in rec.features if fea.type == 'CDS']

	seq_records=[]
	for fea in features:
	    #print(fea.strand)
            name=fea.qualifiers['locus_tag'][0]
	    start= int(re.sub('\W','',str(fea.location.start)))
	    end= int(re.sub('\W','', str(fea.location.end)))
	    seq= rec.seq[start:end]
	    fasta_header='{}:{}-{}'.format(chromosome, str(start+1), str(end)) # notice that the coordinate should be gff style as the bcftools uses it
	    seq_rec= SeqRecord(seq, id= fasta_header, name= name, description= '')
	    seq_records.append(seq_rec)
        with open(output[0], 'w') as out_fh:
            SeqIO.write(seq_records, out_fh, 'fasta')

rule map2reference:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".bwt"
    output:
        TMP_D+"/{strain}/bwa.sorted.bam",
        TMP_D+"/{strain}/bwa.sorted.bam.bai"
    #threads: CORES 
    shell:
        "bwa mem -v 2 -M -t {CORES} {input[REF]} {input[READS1]} {input[READS2]}| samtools view -b -@ {CORES} | bamtools sort > {output[0]} ;"
        "samtools index {output[0]}"

rule call_var:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM=TMP_D+"/{strains}/bwa.sorted.bam",
        BAM_INDEX=TMP_D+"/{strains}/bwa.sorted.bam.bai"
    output:
        TMP_D+"/{strains}/bwa.vcf"
    #threads: CORES
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input[REF_FA_INDEX]} 100000) {CORES} -f {input[REF]} {input[BAM]} > {output};"

rule compress_vcf:
    input:
        TMP_D+"/{strain}/bwa.vcf"
    output:
        TMP_D+"/{strain}/bwa.vcf.gz",
        TMP_D+"/{strain}/bwa.vcf.gz.tbi"
    shell:
        "bgzip {input}; "
        "tabix -p vcf {output[0]}"

rule process_vcf:
    input:
        TMP_D+"/{strain}/bwa.vcf.gz"
    output:
        TMP_D+"/{strain}/bwa.vcf.non-indel.gz",
        TMP_D+"/{strain}/bwa.vcf.non-indel.gz.tbi"
    shell:
        "vcftools --remove-indels --gzvcf {input} --recode --stdout| bgzip -c > {output[0]}; "
        "tabix -p vcf {output[0]}"

rule consensus_whole_genome:
    input:
        REF=REF_FA,
        VCF=TMP_D+"/{strains}/bwa.vcf.non-indel.gz",
    output:
        TMP_D+"/{strains}/consensus_genome.fa"
    shell: 
        "bcftools consensus -f {input[REF]} {input[VCF]} > {output}"

rule whole_genome_alignment:
    input:
        expand(TMP_D+"/{strains}/consensus_genome.fa", strains= STRAINS)
    output:
        TMP_D+'/whole_genomes.aln' 
    shell:
        'cat {input} > {output}'

rule consensus_coding_regions:
    input:
        REGIONS= TMP_D+"/reference.target_regions.fa",
        VCF= TMP_D+"/{strain}/bwa.vcf.gz",
        VCF_INDEX= TMP_D+"/{strain}/bwa.vcf.gz.tbi"
    output:
        TMP_D+"/{strain}/target_regions.fa"
    shell:
        "bcftools consensus -f {input[REGIONS]} {input[VCF]} > {output}"


rule coding_regions_aln:
    input:
        reg_files= expand(TMP_D+"/{strain}/target_regions.fa", strain=STRAINS)
    output:
        TMP_D+'/OneBig.aln' 
    run:
        from Bio import SeqIO
        import textwrap
        from Bio.Align.Applications import MafftCommandline
        from Bio import AlignIO

        ## all the files
        #files={'CH2500': 'tmp/CH2500/target_regions.fa', 'CH2522': 'tmp/CH2522/target_regions.fa', 'CH2502': 'tmp/CH2502/target_regions.fa', 'F1659': 'tmp/F1659/target_regions.fa'}
        files= {STRAINS[n]: input['reg_files'][n] for n in range(len(STRAINS))}
        seq_format= 'fasta'
        #cores= 20
        cores= CORES
        #one_big_aln_f= 'OneBig.aln'
        one_big_aln_f= output[0]


        ## concatenate and then split
        seq_dict= {}## the values are arrays, in which the odd ones are strain names and the even ones are sequences
        for s in files:
            records_dict= SeqIO.to_dict(SeqIO.parse(files[s], seq_format))
            for key in records_dict :
                if not (key in seq_dict):
                    seq_dict[key]= []
                seq_dict[key].append('>'+s) # formated in fasta
                seq_dict[key].append(textwrap.fill(str(records_dict[key].seq), width= 60)) # formated in fasta

        alignments= []
        for k in seq_dict:
            out_fasta= k+'.fa'
            out_aln= k+'.aln'
            alignments.append(out_aln)
            with open(out_fasta, 'w') as out_fh:
                out_fh.write('\n'.join(seq_dict[k]))
            aln=  MafftCommandline(quiet=True, retree= 1, thread= cores, nuc= True, globalpair= True, input= out_fasta)
        #        print(aln())
            with open(out_aln, 'w') as out_fh:
                out_fh.write('\n'.join(aln()))

        one_big_aln= AlignIO.read(alignments[0], 'fasta')
        one_big_aln.sort()
        for f in alignments:
            aln= AlignIO.read(f, 'fasta')
            aln.sort()
            one_big_aln= one_big_aln+aln

        with open(one_big_aln_f, 'w') as out_fh:
            AlignIO.write(one_big_aln, out_fh, 'fasta')

rule convert_gpa_table:
    input:
        TMP_D+'/roary/gene_presence_absence.csv'
    output:
        'gpa.mat'
    script: 'roary_gpa2bin.R'

rule run_roary:
    input:
        expand(TMP_D+"/annotations/{strain}/{strain}.gff", strain= STRAINS)
    output:
        TMP_D+'/roary/gene_presence_absence.csv'
    shell:
        "rm -r {TMP_D}/roary"
        "roary -f {TMP_D}/roary -e -n -v -r -p CORES -g 100000 -z {input}"

rule run_prokka:
    input:
        scaf_f=TMP_D+'/assemblies/{strain}/scaffolds.fasta'
    output:
        TMP_D+"/annotations/{strain}/{strain}.gff"
    params:
        OUT_DIR= TMP_D+"/annotations/{strain}"
    shell:
        "prokka --cpus {CORES} --force --prefix {wildcards.strain} --outdir {params.OUT_DIR} {input}"

rule make_assembly:
    input:
        lambda wildcards: [SAMPLES_DF.loc[wildcards.strain, 'reads1'], SAMPLES_DF.loc[wildcards.strain, 'reads2']]
    output:
        scaf_f=TMP_D+'/assemblies/{strain}/scaffolds.fasta'
    shell:
        "spades.py --careful --threads {CORES} -o {TMP_D}/assemblies/{wildcards.strain} -1 {input[0]} -2 {input[1]}"

