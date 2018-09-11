'''
Separate the indels and SNPs. For the indels, filter the records by quality.
For the SNPs, detect the syn and non-syn variants

The syn and non-syn distinguishment require the sequence, so it seems better to
change the IO.
'''
import subprocess
import pandas as pd
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
import re

def is_syn_mut(vcf_record, regional_ref_seq, fea):
    '''
    Purpose: 
    Distiguish syn and non-syn SNPs
    
    Inputs:
    vcf_record: one row in the vcf table, column names are available
    regional_ref_seq: the Seq object of region described by fea
    fea: the target region described in the input genbank file
    '''
    # Move the coordinates to have the caculation region-based
    relative_pos= int(vcf_record['POS']) - (fea.location.start+1)

    # mutate the sequence
    tmp=list(str(regional_ref_seq))
    tmp[relative_pos]= vcf_record['ALT']
    regional_alt_seq=Seq(''.join(tmp), IUPACAmbiguousDNA())

    # the outcome of mutation
    ref_aa= str(regional_ref_seq[(relative_pos//3
        * 3):((relative_pos//3+1)*3)].translate())
    alt_aa= str(regional_alt_seq[(relative_pos//3
        * 3):((relative_pos//3+1)*3)].translate())
    is_syn= (True if ref_aa == alt_aa else False)

    return(pd.DataFrame({'ref_aa': [ref_aa], 'alt_aa': [alt_aa], 'is_syn':
        [is_syn]}))



def parse_bcftools_output(bcftools_out):
    '''
    Purpose:
    Parse the result of bcftools

    Inputs:
    The raw outcome piped from bcftools, which is bytes instead of
    strings
    '''
    bcftools_out_str=bcftools_out.stdout.decode("utf-8")
    regional_vars_list= [row.split('\t') for row in bcftools_out_str.split('\n')
            if ((re.search('^#', row) is None) and (re.search('\w', row)))]
    if len(regional_vars_list) < 1:
       return(pd.DataFrame()) 

    vcf_cols= ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
            'FORMAT']
    vcf_rename_dict={n: vcf_cols[n] for n in range(len(vcf_cols))}
    regional_vars_dict= {n: regional_vars_list[n] for n in
            range(len(regional_vars_list))}
    regional_vars_df= pd.DataFrame(data=regional_vars_dict).transpose()
    regional_vars_df.rename(columns= vcf_rename_dict, inplace=True)
    return(regional_vars_df)

gbk_f= '../data/reference/RefCln_UCBPP-PA14.gbk'
merged_vcf= 'vcf.gz'

# read the gbk
# for the coordinates and the sequences
rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
chromosome= rec.id
# filter target regions types (gene, CDS...etc)
features= [fea for fea in rec.features if fea.type == 'gene']

# coordinates in bcftools are '1-based' and inclusive
for fea in features:
    gene_name=fea.qualifiers['name'][0] if 'name' in fea.qualifiers else '.'
    gene_id=fea.qualifiers['locus_tag'][0] if 'locus_tag' in fea.qualifiers else '.'
    start= int(re.sub('\W','',str(fea.location.start+1)))
    end= int(re.sub('\W','', str(fea.location.end)))
    coord_str= '{}:{}-{}'.format(chromosome, start, end)
    bcftools_cmd= ['bcftools', 'filter', '-r', coord_str, merged_vcf]
    regional_vars= subprocess.run(bcftools_cmd,
            stdout=subprocess.PIPE) 
    regional_vars_df= parse_bcftools_output(regional_vars)
    regional_ref_seq= rec.seq[fea.location.start:fea.location.end]

    if regional_vars_df.empty:
        continue

    ## extract snps
    snp_mask= regional_vars_df['REF'].apply(lambda x: (len(x) == 1) & (x != '.'))
    snp_mask= snp_mask & regional_vars_df['ALT'].apply(lambda x: (len(x) == 1) & (x != '.'))
    snps_df= regional_vars_df[snp_mask]
    indel_df= regional_vars_df[snp_mask.apply(lambda x: not x)]

    ## translate the codon and
    ## distinguish syn and non-syn snps
    conv_out= pd.concat(snps_df.apply(lambda r: is_syn_mut(r, regional_ref_seq,
        fea),axis=1).tolist()).reset_index()
    snps_df= pd.concat([snps_df, conv_out], axis= 1)
    ## set the feature names
    snps_df['name']= snps_df.apply(lambda r: '|'.join([gene_name, r['POS'], r['REF'], r['ALT'],
        r['ref_aa'], r['alt_aa']]), axis= 1)
    snps_df.set_index('name', inplace= True)
    syn_snps_df= snps_df[snps_df['is_syn']]
    nonsyn_snps_df= snps_df[snps_df['is_syn'].apply(lambda x: not x)]
    

'''
# read the vcf
# directly piped from bcf outcome instead of reading the file
fam_vcf=''
vcf_columns= ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
        'FORMAT', 'unknown']
fam_df= pd.read_csv(fam_vcf, sep= '\t', comment= '#', dtype= str, 
        names=vcf_columns)

# extract snps
# Single substitution
snp_mask= fam_df['REF'].apply(lambda x: (len(x) == 1) & (x != '.'))
snp_mask= snp_mask & fam_df['ALT'].apply(lambda x: (len(x) == 1) & (x != '.'))
snps_df= fam_df[snp_mask]
indel_df= fam_df[snp_mask.apply(lambda x: not x)]

# For snps, find syn mutations
table = CodonTable.ambiguous_dna_by_id[1]
'''
