class Samples:
    def __init__(self, user_opt):
        import sys
        import os
        import InternalVar.ParseSamplesTab as pst

        # dna
        self.DNA_READS= pst.read_sampletab(user_opt['dna_reads'])
        self.DNA_READS_NUM= self.DNA_READS.shape[0]
        # rna
        self.RNA_READS= pst.read_sampletab(user_opt['rna_reads'])
        # phenotypes
        self.PHE_TABLE_F= user_opt['phe_table']

    def __str__(self):
        return('{} dna samples\n{} rna sample\nphenotype file: {}'.format(
            str(self.DNA_READS.shape[0]), 
            str(self.RNA_READS.shape[0]),
            self.PHE_TABLE_F))

    def call_dna_reads(self):
        return(self.DNA_READS)

    def call_rna_reads(self):
        return(self.RNA_READS)

    def call_phenotype_file(self):
        return(self.PHE_TABLE_F)

