'''
Separate the indels and SNPs. For the indels, filter the records by quality.
For the SNPs, detect the syn and non-syn variants

The syn and non-syn distinguishment require the sequence, so it seems better to
change the IO.
'''
import gzip
import subprocess
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re
import numpy as np
import datetime


def vcf_to_feat_table(vcf_df, target_samples):
    feat_tab= vcf_df.apply(lambda r: convert_to_bin(r, target_samples),
            axis=1).transpose()
    return(feat_tab)

def convert_to_bin(vcf_row, strains):
    '''
    Purpose:
    Convert the vcf data into binary features

    Input:
    One row in the vcf table, and it should be a series
    '''
    name= vcf_row.name
    keys=str(vcf_row['FORMAT']).split(':')
    bin_d= {}
    for strain in strains:
        values=str(vcf_row[strain]).split(':') 
        d= {keys[n]: str(values[n]) for n in range(len(values))}
        out= '.'
        if d['GT'] == '.':
            out= np.nan
        elif d['GT'] == '0':
            out= '0'
        else:
            out= '1'
        bin_d[strain]=out
    return(pd.Series(data= bin_d, name= name))

def create_bed3_string(fea):
    '''
    Purpose:
    Create the bed-style coordinate string

    Input:
    The feature object of genbank, which should include the qualifiers "name"
    and "locus_tag"
    '''
    gene_name=fea.qualifiers['name'][0] if 'name' in fea.qualifiers else '.'
    gene_id=fea.qualifiers['locus_tag'][0] if 'locus_tag' in fea.qualifiers else '.'
    ## In bed format, it's 0-based starting for starting site while 1-based for
    ## ending site, which is exactly the parsing result of biopython
    start= re.sub('\W','',str(fea.location.start))
    end= re.sub('\W','', str(fea.location.end))
    coord_str= '\t'.join([chromosome, start, end])
    return(coord_str)

def extract_snps(regional_vars_df, reset_index):
    '''
    Purpose:
    Separate snp and non-snps

    Input:
    The vcf-like dataframe 
    '''
    snp_mask= regional_vars_df['REF'].apply(lambda x: (len(x) == 1) & (x != '.'))
    snp_mask= snp_mask & regional_vars_df['ALT'].apply(lambda x: (len(x) == 1) & (x != '.'))
    snps_df= regional_vars_df[snp_mask]
    indel_df= regional_vars_df[snp_mask.apply(lambda x: not x)]
    if reset_index:
        snps_df.reset_index(drop=True, inplace= True)
        indel_df.reset_index(drop=True, inplace= True)
    return(snps_df, indel_df)

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
    regional_alt_seq=Seq(''.join(tmp), IUPAC.unambiguous_dna)

    # the outcome of mutation
    ref_aa= str(regional_ref_seq[(relative_pos//3
        * 3):((relative_pos//3+1)*3)].translate())
    alt_aa= str(regional_alt_seq[(relative_pos//3
        * 3):((relative_pos//3+1)*3)].translate())
    is_syn= (True if ref_aa == alt_aa else False)

    return(pd.DataFrame({'ref_aa': [ref_aa], 'alt_aa': [alt_aa], 'is_syn':
        [is_syn]}))

def create_vcf_dataframe(vcf_rows):
    '''
    Purpose:
    Igenore the headers but retain the header to create a vcf dataframe
    '''
    vcf_col_str= ''
    vcf_list= []
    for line in vcf_rows:
        if (re.search('^#',line)):
            vcf_col_str= line
        elif (re.search('^#',line) is None) and re.search('\w', line):
            vcf_list.append(line.strip().split('\t'))
        else:
            continue

    vcf_col_str=re.sub('^#+', '', vcf_col_str)
    vcf_col= vcf_col_str.split('\t')
    vcf_col_dict= {n:vcf_col[n] for n in range(len(vcf_col))}
    vcf_df= pd.DataFrame(data=
            {m: vcf_list[m] for m in range(len(vcf_list))}).transpose()
    vcf_df.rename(columns= vcf_col_dict, inplace= True)

    return(vcf_df)


def parse_bcftools_output(bcftools_out):
    '''
    Purpose:
    Parse the result of bcftools

    Inputs:
    The raw outcome piped from bcftools, which is bytes instead of
    strings
    '''
    bcftools_out_str=bcftools_out.stdout.decode("utf-8")
    vcf_like_df=create_vcf_dataframe(bcftools_out_str.split('\n'))
    return(vcf_like_df)

gbk_f= snakemake.input['ref_gbk']
coding_vcf=snakemake.input['coding_vcf_gz']
igr_vcf=snakemake.input['igr_vcf_gz'] 
strains= snakemake.params['strains']
bcftools_bin= snakemake.params['bcftools_bin']
syn_out= snakemake.output['syn_snps_table']
nonsyn_out=snakemake.output['nonsyn_snps_table']


#####
## for coding regions
#####
# read the gbk
# to find the coordinates and the sequences of genes
rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
chromosome= rec.id
# filter by target regions
features= [fea for fea in rec.features if fea.type == 'gene']

full_syn_snps_df= pd.DataFrame({})
full_nonsyn_snps_df= pd.DataFrame({}) 
for fea in features:
    # mask variants by coding regions
    gene_name=fea.qualifiers['name'][0] if 'name' in fea.qualifiers else '.'
    gene_id=fea.qualifiers['locus_tag'][0] if 'locus_tag' in fea.qualifiers else '.'
    # be careful about the coordinates (i.e. 1-baed vs 0-based)
    start= int(re.sub('\W','',str(fea.location.start+1)))
    end= int(re.sub('\W','', str(fea.location.end)))
    coord_str= '{}:{}-{}'.format(chromosome, start, end)
    bcftools_cmd= [bcftools_bin, 'filter', '-r', coord_str, coding_vcf]
    regional_vars= subprocess.run(bcftools_cmd, stdout=subprocess.PIPE) 
    regional_vars_df= parse_bcftools_output(regional_vars)
    regional_ref_seq= rec.seq[fea.location.start:fea.location.end]

    if regional_vars_df.empty:
        continue

    # extract snps
    snps_df, indel_df= extract_snps(regional_vars_df, reset_index= True)
    if snps_df.shape[0] == 0:
        continue

    # translate the codon and
    # distinguish syn and non-syn snps
    conv_out= pd.concat(snps_df.apply(lambda r: is_syn_mut(r, regional_ref_seq,
        fea),axis=1).tolist()).reset_index(drop= True)
    snps_df= pd.concat([snps_df, conv_out], axis= 1)

    # assign the feature names
    snps_df['name']= snps_df.apply(lambda r: '|'.join([gene_name,str(r['POS']), r['REF'], r['ALT'],
        r['ref_aa'], r['alt_aa']]), axis= 1)
    snps_df.set_index('name', inplace= True)
    syn_snps_df= snps_df[snps_df['is_syn']]
    nonsyn_snps_df= snps_df[snps_df['is_syn'].apply(lambda x: not x)]

    full_syn_snps_df=pd.concat([full_syn_snps_df, syn_snps_df], axis= 0)
    full_nonsyn_snps_df=pd.concat([full_nonsyn_snps_df, nonsyn_snps_df], axis= 0)
    
syn_bin_tab= vcf_to_feat_table(full_syn_snps_df, strains) 
nonsyn_bin_tab= vcf_to_feat_table(full_nonsyn_snps_df, strains) 

#####
## for intergenic regions
#####
# read the intergenic vcf file (gz)
igr_vcf_h= gzip.open(igr_vcf, 'r')
igr_vcf_df= create_vcf_dataframe(igr_vcf_h.read().decode("utf-8").split('\n'))
# extract snps
igr_snps_df, igr_indel_df= extract_snps(igr_vcf_df, reset_index= True)
igr_snps_bin_tab= pd.DataFrame({})
# assign the feature name
if igr_snps_df.shape[0] > 0:
    region_name= 'intergenic'
    igr_snps_df['name']= igr_snps_df.apply(lambda r: '|'.join([region_name, r['POS'],
        r['REF'],r['ALT']]), axis=1)
    igr_snps_df.set_index('name', inplace= True)
    igr_snps_bin_tab= vcf_to_feat_table(igr_snps_df, strains)


#####
## clean and print the results
#####
## syn mutants include those in the intergenic regions
syn_bin_tab= pd.concat([syn_bin_tab, igr_snps_bin_tab], axis= 1, join= 'outer')
syn_bin_tab=syn_bin_tab.loc[strains, :].fillna('NA')
nonsyn_bin_tab=nonsyn_bin_tab.loc[strains, :].fillna('NA')

syn_bin_tab.to_csv(syn_out, 
    sep= '\t', 
    header=True, index= True)

nonsyn_bin_tab.to_csv(nonsyn_out, 
    sep= '\t', 
    header=True, index= True)
