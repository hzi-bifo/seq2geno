<?xml version="1.0" ?>
<project name="geno2pheno" output="geno2pheno_out">
  <genotype>
    <bin_tables normalization="binary" path="b" transpose="False">geno2pheno</bin_tables>
    <con_tables normalization="numeric" path="c" transpose="False">geno2pheno</con_tables>
    <sequence kmer="6" path="a">geno2pheno</sequence>
  </genotype>
  <phenotype path="p">
</phenotype>
  <phylogentictree path="t">
</phylogentictree>
  <predict name="classX_vs_classY">
    <optimize>scores_f1_1</optimize>
    <eval folds="10" test="0.1">rand</eval>
    <label value="1">resistant</label>
    <label value="0">sensitive</label>
    <label value="2">sensitive</label>
    <model>
      <svm/>
      <rf/>
    </model>
  </predict>
</project>
