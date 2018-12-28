import pandas as pd

        
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
    threads: 1
    conda: ENV_FILES_POOL.find_yaml('spades_3_10_env')
    script: 'run_spades.py'
#    run:
#        import subprocess
#        import sys
#        out_dir= params.spades_outdir
#        cmd= []
#        try:
#            if len(input.READS) == 2:
#                cmd= [params.SPADES_BIN, params.SPADES_OPT, '--threads',
#                    str(threads), '-o', out_dir, 
#                    '-1', input.READS[0], 
#                    '-2', input.READS[1]]
#            elif len(input.READS) == 1:
#                cmd= [params.SPADES_BIN, params.SPADES_OPT, '--threads',
#                    str(threads), '-o', out_dir, 
#                    '-s', input.READS[0]]
#            else:
#                raise Exception('The reads are neither single nor paired')
#            subprocess.run(cmd)
#        except Exception as e:
#            print(e)
#            sys.exit()
