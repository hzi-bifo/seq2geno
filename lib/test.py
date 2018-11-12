import sys
sys.path.insert(0,'/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/lib')
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
software_pool= es.SoftwarePool(env_sensitive= True)
print(software_pool)
software_pool.find_software('samtools', include_env= True)

os.environ["PATH"] = '{}:{}'.format('/home/thkuo/miniconda3/envs/seq2geno/bin', os.environ["PATH"])
print(os.environ["PATH"]) 
print(software_pool)
software_pool.find_software('samtools', include_env= True)
