
one_big_aln= AlignIO.read(alignments[0], 'fasta')
one_big_aln.sort()
for f in alignments:
    aln= AlignIO.read(f, 'fasta')
    aln.sort()
    one_big_aln= one_big_aln+aln

with open(one_big_aln_f, 'w') as out_fh:
    AlignIO.write(one_big_aln, out_fh, 'fasta')
