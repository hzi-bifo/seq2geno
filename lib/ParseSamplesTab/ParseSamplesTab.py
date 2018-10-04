'''
Input:
    The input should be a two-column tsv table, where the left one includes
 samples/strains and the right one lists the seqeuncing reads. 

Output:
    A pandas series
'''
import pandas as pd
import re

def parse_list(f, reads_source):
    '''
    Read the file and create a pandas series
    '''
    fh= open(f, 'r')
    files_dict= {}
    for l in fh:
        if re.search('\w', l) is None:
            continue
        d= l.strip().split('\t')
        files_dict[d[0]]= d[1].split(',')
    s= pd.Series(files_dict, name=reads_source)
    return(s)

def exclude_rule(x):
    outcome= False

    # by number of files
    if len(x) > 2:
        outcome= True
    # by file name
    if any(re.search('\w', f) is None for f in x ):
        outcome= True
    return(outcome)


def read_sampletab(f, reads_source= 'dna'):
    #f= 'test.samples'
    #reads_source= 'dna'

    ## parse the file
    s= parse_list(f, reads_source)

    ## check the files
    bad_samples= s[s.apply(exclude_rule)].index.values.tolist()
    print('File checking ({})...'.format(f))
    if len(bad_samples) > 0 :
        print('[ERROR] Please check {}'.format(','.join(bad_samples)))
        exit()
    else:
        print('done')
        return(s)
