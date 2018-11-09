import sys
sys.path.insert(0,'/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v9/lib')
import InternalVar

data= {'dna_reads':
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.dna.tsv',
        'rna_reads':
        '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/merged/data/samples.3.rna.tsv',
        'phe_table': 
        'data/strains/pheno/phenotypes.edit.mat'}
main_var= InternalVar.reader_for_main()
main_var.load_materials(data)
