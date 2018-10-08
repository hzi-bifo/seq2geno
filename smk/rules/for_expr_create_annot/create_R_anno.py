
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re

'''
Only the column "gene_name" will be used in the subsequently steps...

'''

gbk_f= snakemake.input['ref_gbk']
#species_strain= snakemake.params['species']
out_f=snakemake.output['R_anno_f']
#gbk_f= 'data/reference/RefCln_UCBPP-PA14.gbk'
#species_strain= 'Paeru_PA14'
#out_f='test.tab'

# read the gbk
rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
chr_len= str(len(rec.seq))
acc= rec.id

with open(out_f, 'w') as out_fh:
    columns=['locus', 'ID', 'gene_name', 'type', 'Start', 'End']
    out_fh.write('\t'.join(columns)+'\n')


    for fea in [fea for fea in rec.features if ((fea.type == 'gene') or
        (not (re.search('RNA', fea.type) is None)))]:
        if not ('locus_tag' in fea.qualifiers):
            continue
        acc= fea.qualifiers['locus_tag'][0]
        gene_name= (fea.qualifiers['name'][0] if 'name' in fea.qualifiers else
                '')
        locus= ','.join([acc, gene_name])
        # be careful about the coordinates
        d=[locus, acc, gene_name, str(fea.type), str(fea.location.start+1),
                str(fea.location.end)]
        out_fh.write('\t'.join(d)+'\n')

