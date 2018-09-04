
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import re


gbk_f= snakemake.input['ref_gbk']
species_strain= snakemake.params['species']
out_f=snakemake.output['anno_f']
#gbk_f= 'data/reference/RefCln_UCBPP-PA14.gbk'
#species_strain= 'Paeru_PA14'
#out_f='test.tab'

# read the gbk
rec= SeqIO.read(open(gbk_f, 'r', encoding='windows-1252'), 'gb')
chr_len= str(len(rec.seq))
acc= rec.id
sequence_type= 'Chromosome'

with open(out_f, 'w') as out_fh:
    header_line= '@'+'|'.join([species_strain, acc, chr_len])
    columns=['@Strain', 'Refseq_Accession', 'Replicon', 'Locus_Tag', 'Feature_Type',
    'Start', 'Stop', 'Strand', 'Gene_Name']
    out_fh.write(header_line+'\n')
    out_fh.write('\t'.join(columns)+'\n')

    target_feature= 'gene'
    for fea in [fea for fea in rec.features if fea.type == target_feature]:
        if not ('locus_tag' in fea.qualifiers):
            continue
        d=[species_strain, acc, sequence_type, fea.qualifiers['locus_tag'][0],
                target_feature, str(fea.location.start), str(fea.location.end),
                '+' if fea.strand == 1 else '-', fea.qualifiers['name'][0] if
                'name' in fea.qualifiers else '.']
        out_fh.write('\t'.join(d)+'\n')

