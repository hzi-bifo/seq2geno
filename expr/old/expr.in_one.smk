import pandas as pd

list_f= config['list_f']
rna_reads= {}
with open(list_f, 'r') as list_fh:
    for l in list_fh:
        d=l.strip().split('\t')
        rna_reads[d[0]]= d[1]

strains=list(rna_reads.keys())
out_table= config['out_table']
out_log_table= config['out_log_table']
ref_fasta= config['ref_fasta']
annot_tab= config['annot_tab']
r_annot= config['r_annot']

rule:
    input:
        out_table,
        out_log_table

rule collect_rpg:
    input:
        rpg_files= expand('{strain}.rpg',strain= strains)
    output:
        rpg_tab=out_table,
        log_rpg_tab=out_log_table
    wildcard_constraints:
        strain='^[^\/]+$'
    run:
        import pandas as pd
        import re
        # import all files
        series_to_collect=[]
        for f in input.rpg_files:
            s= re.search('([^\/]+).rpg', f).group(1)
            df=pd.read_csv(f, sep='\t', header= None)
            col_dict={0: 'locus', 1: s, 2: 'antisense_count'}
            df=df.rename(columns=col_dict).set_index('locus')
            series_to_collect.append(df.loc[:,s])
        rpg_df=pd.DataFrame(series_to_collect)
        rpg_df.to_csv(output[0], sep='\t') 

        log_rpg_df= pd.np.log(rpg_df+1)
        log_rpg_df.to_csv(out_log_table, sep = "\t")    
        
rule my_stampy_pipeline:
    input:
        #infile=strain_fq,
        infile=lambda wildcards: rna_reads[wildcards.strain],
        reffile=ref_fasta,
        ref_index_stampy=ref_fasta+'.stidx',
        ref_index_bwa=ref_fasta+'.bwt',
        annofile=annot_tab,
        Rannofile=r_annot
    output:
        sam='{strain}.sam',
        art='{strain}.art',
        sin='{strain}.sin',
        flatcount='{strain}.flatcount',
        rpg='{strain}.rpg',
        stat='{strain}.stats'
    threads:1
    conda: 'snps_tab_mapping.yml'
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/perl5/site_perl/5.22.0:\
$CONDA_PREFIX/lib/perl5/5.22.2:\
$CONDA_PREFIX/lib/perl5/5.22.2/x86_64-linux-thread-multi/:\
$PERL5LIB
        my_stampy_pipeline {wildcards.strain} {input.infile} \
{input.reffile} {input.annofile} {input.Rannofile} 2> {wildcards.strain}.log
        """

rule stampy_index_ref:
    input:
        reffile=ref_fasta
    output:
        ref_fasta+'.bwt',
        ref_fasta+'.stidx',
        ref_fasta+'.sthash'
    shell:
        '''
        stampy.py -G {input.reffile} {input.reffile}
        stampy.py -g {input.reffile} -H {input.reffile}
        bwa index -a bwtsw {input.reffile}
        '''
