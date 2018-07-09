
'''
rule all:
    input:
        tree
'''
rule find_best_tree:
    input:
        one_big_var_aln=TMP_D+'/OneBig.var.aln'
    output:
        tree=TREE
    params:
        RAXML= RAXML_EXE,
        PREFIX= 'OneBig.var',
        RESULT_D= RESULT_D,
        RAXML_OUTDIR= RESULT_D,
       	CORES= CORES
    shell:
        '{params[RAXML]} -T {params[CORES]} -w {params[RAXML_OUTDIR]} '
        '-m GTRGAMMA -s {input} -n {params[PREFIX]} -p 1 -N 1;' 
        'cp -s {params[RAXML_OUTDIR]}/RAxML_bestTree.{params[PREFIX]} {output}'

rule postprocess_alignment:
    input:
        one_big_aln='{TMP_D}/OneBig.aln'
    output:
        one_big_var_aln='{TMP_D}/OneBig.var.aln'
    shell:
        "trimal -st 1 -gt 1 -complementary -in {input} -out {output}"

rule create_coding_regions_aln:
### concatenate 
    input:
        cons_coding_seqs_every_strain=expand(
            "{TMP_D}/{strains}/{mapper}.cons.fa", 
            TMP_D= TMP_D, strains= STRAINS, mapper= 'bwa')
    output:
        one_big_aln='{TMP_D}/OneBig.aln'

    params:
        CORES=CORES, 
        TMP_D=TMP_D+"/families", 
        STRAINS=STRAINS 
    script: 'makeAlignment.py'
'''            
rule load_vcf:
    output:
        vcf_gz=TMP_D+"/{strains}/{mapper}.vcf.gz"

rule load_reference
    input:
        REF_FA
    output:
        REF_FA+".fai",
        REF_FA+".bwt"
'''

'''
rule prepare_reference:
    input:
        REF_FA
    output:
        REF_FA+".fai",
        REF_FA+".bwt"

    shell:
        "samtools faidx {input};" 
        "bwa index {input}"

rule prepare_reference_Ariane:
    input:
        REF_FA
    output:
        REF_FA+".stidx",
        REF_FA+".sthash"
    params:
        STAMPY=STAMPY_EXE
    shell:
        "source activate py27;"
        "{params[STAMPY]} -G {input} {input};"
        "{params[STAMPY]} -g {input} -H {input};"
        "source deactivate"

rule prepare_target_regions:
    input:
        REF_GBK
    output:
        #TMP_D+"/reference.target_regions.fa"
        "{tmp_dir}/reference.target_regions.fa"
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

rule map2reference_Ariane:
    input:
        REF=REF_FA,
        READS1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        READS2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF_INDEX=REF_FA+".stidx"
    output:
        "{tmp_dir}/{strain}/stampy.sorted.bam",
        "{tmp_dir}/{strain}/stampy.sorted.bam.bai"
    params:
        STAMPY=STAMPY_EXE,
    	CORES=CORES
    shell:
        # stampy runs with python2.6
        "source activate py27;"
        "{params[STAMPY]} --bwaoptions=\"-q10 {input[REF]}\" -g {input[REF]} -h  {input[REF]} -M  {input[READS1]} {input[READS2]} | samtools view -bS -@ {params[CORES]} |  samtools sort>  {output[0]};"
        "samtools index {output[0]};"
        "source deactivate;"

rule call_var:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{tmp_dir}/{strains}/stampy.sorted.bam",
        BAM_INDEX="{tmp_dir}/{strains}/stampy.sorted.bam.bai"
    output:
        "{tmp_dir}/{strains}/"+PREFIX+".vcf"
    params: 
        CORES=CORES
    shell:
        "freebayes-parallel <(fasta_generate_regions.py {input[REF_FA_INDEX]} 100000) {params[CORES]} -f {input[REF]} {input[BAM]} > {output};"

rule compress_vcf:
    input:
        "{tmp_dir}/{strain}/{PREFIX}.vcf"
    output:
        #TMP_D+"/{strain}/bwa.vcf.gz",
        #TMP_D+"/{strain}/bwa.vcf.gz.tbi"
        "{tmp_dir}/{strain}/{PREFIX}.vcf.gz",
        "{tmp_dir}/{strain}/{PREFIX}.vcf.gz.tbi"
    shell:
        "bgzip {input}; "
        "tabix -p vcf {output[0]}"

rule process_vcf:
    input:
        TMP_D+"/{strain}/{PREFIX}.vcf.gz"
    output:
        TMP_D+"/{strain}/{PREFIX}.vcf.non-indel.gz",
        TMP_D+"/{strain}/{PREFIX}.vcf.non-indel.gz.tbi"
    shell:
        "vcftools --remove-indels --gzvcf {input} --recode --stdout| bgzip -c > {output[0]}; "
        "tabix -p vcf {output[0]}"

rule consensus_coding_regions:
    input:
        REGIONS= TMP_D+"/reference.target_regions.fa",
        VCF= TMP_D+"/{strain}/"+PREFIX+".vcf.gz",
        VCF_INDEX= TMP_D+"/{strain}/"+PREFIX+".vcf.gz.tbi"
    output:
        #TMP_D+"/{strain}/target_regions.fa"
        TMP_D+"/{strain}/"+PREFIX+"_target_regions.fa"
    shell:
        "bcftools consensus -f {input[REGIONS]} {input[VCF]} > {output}"

rule coding_regions_aln:
    input:
        reg_files= expand("{tmp_dir}/{strain}/{PREFIX}_target_regions.fa", strain=STRAINS, tmp_dir= TMP_D, PREFIX= PREFIX)
    output:
        TMP_D+'/'+PREFIX+'_OneBig.aln' 
    params:
        PREFIX=PREFIX,
	CORES=CORES, 
	TMP_D=TMP_D, 
	STRAINS=STRAINS 
    script: 'makeAlignment.py'

rule process_aln:
    input:
        TMP_D+'/'+PREFIX+'_OneBig.aln' 
    output:
        TMP_D+'/'+PREFIX+'_OneBig.var.aln'  
    shell:
        "trimal -st 1 -gt 1 -complementary -in {input} -out {output} "

rule find_besttree:
    input:
        TMP_D+'/'+PREFIX+'_OneBig.var.aln' 
    output:
        'RAxML_bestTree.'+PREFIX+'_OneBig.var'
    params:
        RAXML= RAXML_EXE,
        PREFIX= PREFIX+'_OneBig.var',
        RAXML_OUTDIR= RESULT_D,
    	CORES= CORES
    shell:
        '{params[RAXML]} -T {params[CORES]} -w {params[RAXML_OUTDIR]} -m GTRGAMMA -s {input} -n {params[PREFIX]} -p 1 -# 1' 
'''
