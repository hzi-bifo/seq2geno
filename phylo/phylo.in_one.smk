import os
import pandas as pd

list_f=config['list_f']
dna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        dna_reads[d[0]]= d[1].split(',')

REF_FA=config['REF_FA']
REF_GFF=config['REF_GFF']
results_dir=config['results_dir']
seq_list=config['seq_list']
families_seq_dir=config['families_seq_dir']
aln_f=config['aln_f']
tree_f=config['tree_f']
adaptor_f= config['adaptor']
new_reads_dir= config['new_reads_dir']

rule all:
    input:
        tree_f

rule tree:
    input:
        alignment=aln_f
    output:
        tree= tree_f
    threads:20
    shell:
        '''        
        which FastTreeMP
        export OMP_NUM_THREADS={threads}
        FastTreeMP -nt -gtr -gamma \
-log {output.tree}.log -out {output.tree} {input.alignment}
        '''

rule process_aln:
    input:
        out_aln='OneLarge.aln'
    output: 
        out_var_aln=aln_f
    params:
        tmp_aln1= 'OneLarge.gapReplaced.aln',
        tmp_aln2= 'OneLarge.gapReplaced.var2.aln'
    conda: 'concatenate_seq_env.yml'    
    shell:
        """
        cat {input} | sed '/^>/!s/[^ATCGatcg]/-/g' > {params.tmp_aln1}
        removeInvariant.py --in {params.tmp_aln1} --out {params.tmp_aln2} --cn 2
        trimal -gt 0.90 -in {params.tmp_aln2} -out {output.out_var_aln}
        """
        
rule concatenate:
    input:
        fam_list='aln_to_concatenate'
    output:
        out_aln='OneLarge.aln'
    conda: 'concatenate_seq_env.yml'    
    params:
        conc_script='concatenateAln.py'
    shell:
        '''
        {params.conc_script} --l {input.fam_list} --o {output.out_aln} 
        '''

rule list_families:
    input:
        out_seq=dynamic(os.path.join(families_seq_dir, '{fam}.aln'))
    output:
        out_fam_list='aln_to_concatenate'
    run:
        with open(output.out_fam_list, 'w') as out_fh:
            for f in input.out_seq:
                out_fh.write(f+'\n')

rule alignment:
    input:
        os.path.join(families_seq_dir, '{fam}.fa')
    output:
        os.path.join(families_seq_dir, '{fam}.aln')
    conda: 'mafft_7_310_env.yml'
    shell:
        '''
        mafft --nuc --maxiterate 0 --retree 2 --parttree \
{input} > {output}
        '''

rule sort:
    input:
        seq_list_f=seq_list
    output:
        out_seq=dynamic(os.path.join(families_seq_dir, '{fam}.fa'))
    conda: 'concatenate_seq_env.yml'    
    params:
        out_stats='fam_stats.txt',
        out_dir= families_seq_dir
    shell:
        '''
        geneStats.py --l {input.seq_list_f} --d {params.out_dir} \
--o {params.out_stats}
        '''

rule list_cons_sequences:
    input:
        expand('{tmp_d}/{strain}/target_regions.fa', tmp_d= results_dir, strain=list(dna_reads.keys()))
    output:
        seq_list_f=seq_list
    params: 
        tmp_d= results_dir,
        strains= list(dna_reads.keys())
    shell:
        '''
        parallel 'echo {{}}"\t"{params.tmp_d}/{{}}/target_regions.fa' ::: {params.strains} \
> {output.seq_list_f}
        '''

rule consensus_seqs:
    input:
        ref_region_seqs=REF_GFF+'.gene_regions.fa',
        vcf_gz='{tmp_d}/{strain}/bwa.vcf.gz',
        vcf_gz_index= '{tmp_d}/{strain}/bwa.vcf.gz.tbi'
    output:
        cons_fa='{tmp_d}/{strain}/target_regions.fa'
    conda:'bcftools_1_6_env.yml'
    shell:
        '''
        bcftools consensus -f {input.ref_region_seqs} \
-o {output.cons_fa} {input.vcf_gz} 
        '''
        

rule ref_regions: 
    input:
        ref_gff=REF_GFF,
        ref_fa=REF_FA
    output:
        ref_regions=REF_GFF+".gene_regions",
        ref_region_seqs=REF_GFF+".gene_regions.fa"
    conda: 'py27.yml'
    shell:
        '''
        make_GeneRegions.py --g {input.ref_gff} --f gene > {output.ref_regions}
        parallel -j 1 \'samtools faidx {input.ref_fa} {{}}\' \
 :::: {output.ref_regions} > {output.ref_region_seqs}
        '''

rule index_ref:
    input:
        REF_FA  
    output:
        REF_FA+".bwt",
        REF_FA+".fai"
    threads:1
    shell:
        """
        samtools faidx {input}
        bwa index {input}
        """

rule mapping:
    input:
        ref=REF_FA,
        ref_index=REF_FA+".bwt",
        FQ1= os.path.join(new_reads_dir, '{strain}.cleaned.1.fq.gz'),
        FQ2= os.path.join(new_reads_dir, '{strain}.cleaned.2.fq.gz')
    output:
        sam="{tmp_d}/{strain}/bwa.sam",
        bam="{tmp_d}/{strain}/bwa.bam",
        sorted_bam="{tmp_d}/{strain}/bwa.sorted.bam",
        sorted_bam_index="{tmp_d}/{strain}/bwa.sorted.bam.bai"
    threads:1
    shell:
        """
        bwa mem -v 2 -M -R \'@RG\\tID:snps\\tSM:snps\' -t {threads} \
 {input.ref} {input.FQ1} {input.FQ2} > {output.sam}
        samtools view -b -S {output.sam} -o {output.bam} -@ {threads}
        bamtools sort -in {output.bam} -out {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

rule redirect_and_preprocess_reads:
    input: 
        infile1=lambda wildcards: dna_reads[wildcards.strain][0],
        infile2=lambda wildcards: dna_reads[wildcards.strain][1]
    output:
        log_f= os.path.join(new_reads_dir, '{strain}.log'),
        f1= os.path.join(new_reads_dir, '{strain}.cleaned.1.fq.gz'),
        f2= os.path.join(new_reads_dir, '{strain}.cleaned.2.fq.gz')
    params:
        adaptor_f= adaptor_f,
        tmp_f1= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.1.fq'.format(wildcards.strain)),
        tmp_f2= lambda wildcards: os.path.join(
            new_reads_dir, '{}.cleaned.2.fq'.format(wildcards.strain))
    shell:
         '''
        if [ -e "{params.adaptor_f}" ]
        then
            fastq-mcf -l 50 -q 20 {params.adaptor_f} {input.infile1} {input.infile2} \
-o {params.tmp_f1} -o {params.tmp_f2} > {output.log_f}
            gzip -9 {params.tmp_f1}
            gzip -9 {params.tmp_f2}
        else
            echo 'Reads not trimmed'
            echo 'No trimming' > {output.log_f}
            echo $(readlink {input.infile1}) >> {output.log_f}
            echo $(readlink {input.infile2}) >> {output.log_f}
            ln  -s {input.infile1} {output.f1}
            ln  -s {input.infile2} {output.f2}
        fi
        '''

rule call_var:
    input:
        ref=REF_FA,
        ref_index=REF_FA+".fai",
        sorted_bam="{tmp_d}/{strain}/bwa.sorted.bam",
        sorted_bam_index="{tmp_d}/{strain}/bwa.sorted.bam.bai"
    output:
        vcf='{tmp_d}/{strain}/bwa.vcf',
        vcf_gz='{tmp_d}/{strain}/bwa.vcf.gz',
        vcf_gz_index= '{tmp_d}/{strain}/bwa.vcf.gz.tbi'
    params:
        freebayes_params= '-p 1'
    threads: 16
    conda:'freebayes_1_1_env.yml'
    shell:
        """
        freebayes-parallel <(fasta_generate_regions.py {input.ref_index} 100000) \
 {threads} -f {input.ref} {params.freebayes_params} {input.sorted_bam} \
 > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf_gz}
        tabix -p vcf {output.vcf_gz}
        """
