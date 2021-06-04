#!/bin/bash

seq2geno_zip=$1
cd $( dirname $seq2geno_zip )
unzip $seq2geno_zip
cd $( basename $seq2geno_zip | sed 's#.zip##' ) 
cat files/_seq2geno_inputs.yml | sed "s#\(\s\+\)__\/#\1$(realpath ./)\/#g" \
  > seq2geno_inputs.yml 

cat files/_dna_list | sed "s#__#$(realpath ./)\/#g" \
  > dna_list

cat files/_rna_list | sed "s#__#$(realpath ./)\/#g" \
  > rna_list
