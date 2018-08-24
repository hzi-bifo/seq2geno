	my_stampy_pipeline
	    - sam2art.pl 
		-s 2 $prefix.sam > $prefix.art
		-s 2 -l $prefix.sam > $prefix.sin
		-f -s 2 $prefix.sam > $prefix.flatcount
	    - art2genecount.pl 
		-b -a $prefix.sin -t tab -r $annotation > $prefix.rpg
		- annotation: path to annotation file for R
	    - sam_statistics.pl 
		-r $prefix.sam > $prefix.stats
	    - genes_statistics.R 
		$prefix.rpg $Rannotation $prefix.stats $prefix
	art2cov.py 
		-f dictionary -s SE -c 1 -o coverage_cut1.txt
		- dictionary: prefix of ".art" and ".rstats"
