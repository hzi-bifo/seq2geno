import sys
sys.path.insert(0,
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/lib')
import InternalVar
import os

data= {'ref_fa': 'ref.fa',
        'ref_gbk': 'ref.gbk',
        'dna_reads':
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.dna.tsv',
        'rna_reads':
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.rna.tsv',
        'phe_table': 
        'data/strains/pheno/phenotypes.edit.mat'}
main_var= InternalVar.reader_for_main(data)
main_var.load_materials()
main_var.load_reference()
main_var.load_library()

import ExternalSoftware as es
software_pool= es.SoftwarePool(env_dir=
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/env')
print(software_pool)
software_pool.find_software('roary')

os.environ["PATH"] = '{}:{}'.format('/home/thkuo/miniconda3/envs/seq2geno/bin', os.environ["PATH"])
print(os.environ["PATH"]) 
print(software_pool)
software_pool.find_software('roary')
software_pool.find_software('roary', target_env=
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/env/roary_env')
software_pool.find_software('roary', target_dir=
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/lib/roary-3.8.2/bin')

import EnvironmentManager 
env_pool= EnvironmentManager.Pool('/home/thkuo/miniconda3/envs/')
env_pool.activate_env_cmd('roary_env_backup')
env_pool.update_variables_cmd('PER5LIB',
        '/home/thkuo/miniconda3/envs/roary_env_backup/lib/perl5/site_perl/5.22.0/:$PERL5LIB')
env_pool.update_variables_cmd('whatever',
        '/home/thkuo/miniconda3/envs/roary_env_backup/lib/perl5/site_perl/5.22.0/:$PERL5LIB',
        check_existence= True)
