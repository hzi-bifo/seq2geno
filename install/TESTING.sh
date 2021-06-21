#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# This script test Seq2Geno with the example dataset.
# In the processes, procedure-specific dependencies will also be installed if needed
# This script DOESN'T activate any conda environment as it evaluates the current environment

check_conda_channels () {
	## ensure conda channels
	echo '+check conda channels...'
	for c in hzi-bifo conda-forge/label/broken bioconda conda-forge defaults; do 
		echo '+'$c
		if [ $(conda config --get channels | grep $c | wc -l) -eq 0 ]; then
			conda config --add channels $c
		fi
	done
	cd $SEQ2GENO_HOME
}

set_roary_dependencies () {
	## Roary dependencies
	cd $SEQ2GENO_HOME/denovo/lib/Roary
	export PERL_MM_USE_DEFAULT=1
	export PERL5LIB=$( realpath . )/lib:$PERL5LIB
	./install_dependencies.sh
	cd $SEQ2GENO_HOME
}

download_proc_specific_env () {
	## decompress the example dataset to install the process-specific environments
  	echo $( realpath . )
	echo '+extract example dataset'
	tar -zxvf example_sg_dataset.tar.gz 
	cd $SEQ2GENO_HOME/example_sg_dataset/
	./CONFIG.sh
	echo '+test the procedures in dryrun mode with the example dataset'
	$SEQ2GENO_HOME/main/seq2geno -f ./seq2geno_inputs.yml || return false
}

#>>>
check_conda_channels || { echo "Errors in setting conda channels"; exit; }
set_roary_dependencies || { echo "Errors in installation of Roary dependecies"; exit; }
download_proc_specific_env || { echo "Errors in installation of the process-specific environments failed"; exit; }

## Finalize
echo 'Finished. You might want to go to '$SEQ2GENO_HOME'/example_sg_dataset/ and try'
