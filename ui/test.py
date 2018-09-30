import ParseSamplesTab as pst

f= '../test.samples'
reads_source= 'dna'
dna_s= pst.read_samtab(f, reads_source)
