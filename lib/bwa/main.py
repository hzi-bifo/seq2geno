'''
rule ng_paired_read_bwa_mapping:
    input:
        FQ1=lambda wildcards: DNA_READS[wildcards.strain][0],
        FQ2=lambda wildcards: (DNA_READS[wildcards.strain][1] if
len(DNA_READS[wildcards.strain]) > 1 else None),
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        P_BWA_SAI1= temp('{TMP_D}/{strain}/stampy/bwa1.sai'),
        P_BWA_SAI2= temp('{TMP_D}/{strain}/stampy/bwa2.sai'),
        BWA_BAM= temp('{TMP_D}/{strain}/stampy/pre_bwa.bam')

    params:
        bwa_exe= 'bwa',
        samtools_exe= 'samtools',
        result_dir= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'strain'),
        BWA_OPT='-q10',
        CORES=CORES 
    run:
'''
def aln(cores, opt, ref, fq, bwa_sai):
    cmd= [BWA_BIN, 'aln', opt, '-t', str(cores), 
    ref, fq]
    with open(bwa_sai, 'w') as bwa_sai_fh:
        subprocess.run(cmd, stdout= bwa_sai_fh)


BWA_BIN=snakemake.params[bwa_exe]
# 'bwa aln' the first part
aln(snakemake.params['CORES'],
     snakemake.params['BWA_OPT'],
    snakemake.input['REF'],
    snakemake.input['FQ1'],
    snakemake.output['P_BWA_SAI1']
)
if not (snakemake.input['FQ2'] is None):
    aln(snakemake.params['CORES'],
         snakemake.params['BWA_OPT'],
        snakemake.input['REF'],
        snakemake.input['FQ2'],
        snakemake.output['P_BWA_SAI2']
    )

        cmd1= [params.bwa_exe, 'aln', params.BWA_OPT, '-t', str(params.CORES), 
        input.REF, input.FQ1]
        with open(output.P_BWA_SAI1, 'w') as bwa_sai1_fh:
            subprocess.run(cmd1, stdout= bwa_sai1_fh)

        if not (input.FQ2 is None):
            cmd2= [params.bwa_exe, 'aln', params.BWA_OPT, '-t', str(params.CORES),
            input.REF, input.FQ2]
            with open(output.P_BWA_SAI2, 'w') as bwa_sai2_fh:
                subprocess.run(cmd2, stdout= bwa_sai2_fh)
            cmd_tosam=[
                params.bwa_exe, 'sampe', '-r',
                "'@RG\tID:{}\tSM:{}'".format(wildcards.strain, wildcards.strain),
                input.REF, output.P_BWA_SAI1, output.P_BWA_SAI2,
                input.FQ1, input.FQ2]
            cmd_tobam=[
                params.samtools_exe, 'view', '-bS', '-@', str(params.CORES)]
            with open(output.BWA_BAM, 'w') as bwa_paired_fh:
                p_tosam= subprocess.Popen(cmd_tosam, stdout= subprocess.PIPE)
                p_tobam= subprocess.Popen(cmd_tobam, stdin= p_tosam.stdout,
                    stdout= bwa_paired_fh)
                p_tosam.stdout.close()
        else:
            with open(output.P_BWA_SAI2, 'w') as bwa_sai2_fh:
                bwa_sai2_fh.write('')
            cmd_tosam= [
                params.bwa_exe, 'samse', '-r',
                "'@RG\tID:{}\tSM:{}'".format(wildcards.strain, wildcards.strain),
                input.REF, input.P_BWA_SAI1, input.FQ1]
            cmd_tobam=[
                params.samtools_exe, 'view', '-bS', '-@', str(params.CORES)]
            with open(params.BWA_BAM, 'w') as bwa_paired_fh:
                p_tosam= subprocess.Popen(cmd_tosam, stdout= subprocess.PIPE)
                p_tobam= subprocess.Popen(cmd_tobam, stdin= p_tosam.stdout,
                    stdout= bwa_paired_fh)
                p_tosam.stdout.close()
