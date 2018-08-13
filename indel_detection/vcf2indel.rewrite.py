'''
Feature: converting the table to put each allele type in one row
'''

import pandas as pd

def is_del(vcf_row):
    # ALT should contain only one allele type
    ref= vcf_row['REF']
    alt= vcf_row['ALT']
    return((len(ref) > len(alt)) | ((ref != '.')&(alt == '.')))

def is_ins(vcf_row):
    # ALT should contain only one allele type
    ref= vcf_row['REF']
    alt= vcf_row['ALT']
    return(len(ref) < len(alt))

def find_genotype(f, values_str):
    fields_list= f.split(':')
    values_list= values_str.split(':')
    fields_dict= {fields_list[n]: values_list[n] for n in range(len(values_list))}
    return(fields_dict['GT'])

def is_gt(gt, f, vcf_status_row):
    geno_types= vcf_status_row.apply(lambda x: find_genotype(f, x))
    return(geno_types == gt)

def row_expanded_by_alt(vcf_row):
    # First list all genotypes, 
    # and then find the strains with each certain genotype.
    # In the outcome, the rows are all the genotypes, which 
    # should be at least 2 for one REF and one ALT, and the 
    # columns includes all the strains. The values of each cell is binary. The
    # value of (M, N) describes whether an allele type M was found in strain N.
    status_format=vcf_row['FORMAT']
    gts_nuc= [str(vcf_row['REF'])]+ vcf_row['ALT'].split(',')
    gts_status= [('{}/{}'.format(str(n), str(n)) if gts_nuc[n] != '.' else './.') for n in range(len(gts_nuc))] 
    status_nuc_dict= {gts_status[n]: gts_nuc[n] for n in range(len(gts_nuc))}
    expanded_status_df= pd.DataFrame(data= {gt: is_gt(gt, status_format,
        vcf_row[9:]) for gt in gts_status}).transpose()
    left= pd.DataFrame([vcf_row[:9]]*len(gts_nuc))
    left['ALT']=[status_nuc_dict[n] for n in expanded_status_df.index.values.tolist()]
    left.reset_index(drop= True, inplace= True)
    expanded_status_df.reset_index(drop= True, inplace= True)
    return(pd.concat([left, expanded_status_df], axis= 1))

def df_expanded_by_alt(vcf_df):
    # Apply the expansion (row_expandedby_alt) to each row
    new_df=pd.concat([row_expanded_by_alt(vcf_df.iloc[n,]) for n in
        range(vcf_df.shape[0])],axis=0)
    return(new_df.reset_index(drop= True))

#def vcf2indel(vcf, gene_name, out, strains_perc_cutoff= 0.5, len_cutoff= 8):
def vcf2indel(vcf, gene_name, out, len_cutoff= 8):
    # Import a vcf -> convert the format -> filter and detect indels -> find
    # the strains having the indels

    # read vcf
    vcf_df= pd.read_csv(vcf, header= 0, index_col= None, sep= '\t')

    # one ALT type, one row
    # convert the genotypes to binary values in the status 
    exp_vcf_df= df_expanded_by_alt(vcf_df)

    # filter variants and find indels
    # case I. ALT is missing value
    # The class 1 are always the minority
    case1_mask= (exp_vcf_df['ALT'] ==
            '.')&(exp_vcf_df.iloc[:,9:].sum(axis=1)>(strain_num/2))
    exp_vcf_df.loc[case1_mask, 9:]= exp_vcf_df.loc[case1_mask,
            9:].applymap(lambda x: not x)
    # case II. ALT is not missing value
    # If the minority has missing values, include them into class 1


    indel_df= exp_vcf_df[exp_vcf_df.apply(is_del, axis= 1) | exp_vcf_df.apply(is_ins, axis= 1)] 
    # by length
    indel_df= indel_df[abs(indel_df['REF'].apply(lambda x: len(x.replace('.',
        ''))) - indel_df['ALT'].apply(lambda x: len(x.replace('.', ''))))
        > len_cutoff]# the dot (.) has length 1 but it's not adequate to 1 nucleotide
    # missing values are included only when it's the minority
    strain_num= indel_df.shape[1]-9
    indel_df= indel_df.drop(indel_df[(indel_df['ALT'] == '.')
        & ((indel_df.iloc[:,9:].sum(axis=1)/strain_num) >= 0.5)].index)
    # remove redundant rows
    indel_df= indel_df[indel_df.iloc[:,9:].sum(axis= 1) > 0]
    
    # find strains with the indels
    status_df= indel_df.iloc[:, 9:]
    outcomes= status_df.sum().apply(lambda x: 1 if x > 0 else 0)

    # write the results
    pd.DataFrame({gene_name: outcomes}).to_csv(out, sep='\t', encoding= 'utf-8')

vcf2indel(snakemake.input[0], snakemake.wildcards['fam'],
        snakemake.output['FAM_INDELS_TXT'],
        strains_perc_cutoff=snakemake.params['strains_perc_cutoff'],
        len_cutoff=snakemake.params['len_cutoff'] )
