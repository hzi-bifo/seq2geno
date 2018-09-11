
fam='dnaA'
input_vcf='./vcf.gz'
input_vcf_index='./vcf.gz.tbi'
bcftools filter -R regions.bed $input_vcf > $fam'.vcf'
