#!/usr/bin/perl -w
# art2genecount.pl
# AUTHOR: Andreas Dötsch, Denitsa Eckweiler
# LAST REVISED: May 2012
# 

## reads art file and writes gene count for DESeq
use strict;
use Getopt::Std;

#variables

my $usage = "\n\nusage: $0 [-t  <reference type> -c] -a <art file> -r <reference file>\n".
            "Reads a reference file of the specified format and and art file and write gene counts for expression analysis.\n\n".
	    "-a <art file>      \tart (artemis readable pileup) file containing information on expression of both strands. Should be created with \"sinister\" option.\n".
	    "-t <reference type>\ttype of reference file used. Paste -t tab for tabulated file.\n".
	    "-r <reference file>\tFile containing the reference genes. Format is specified with the -t option.\n".
	    "-b\t\tBlank mode. Only readcounts are printed in the output without gene IDs and anti-sense counts.\n".
	    "-c                 \tif this flag is set, only coding sequences (no ncRNA, tRNA, rRNA) are considered from the annotation.\n\n";
our($opt_a,$opt_t,$opt_r,$opt_c,$opt_b);
getopts('a:t:r:cb') or die $usage;

if(!defined($opt_a)){die $usage}
if(!defined($opt_t)){die $usage}
if(!defined($opt_r)){die $usage}
if(!defined($opt_b)){ $opt_b = 0 } else { $opt_b = 1 }
my $art_file	= $opt_a; 
my $type 	= $opt_t; 
my $ref_file	= $opt_r;
my $blankmode   = $opt_b;
my $cds_only;
if(defined($opt_c)){$cds_only = 1} else {$cds_only = 0}
my ($genomesize,$tmp,$i,$line,$genetype,$genetxt,$geneID,$n_genes); #skalar variables
my (@forward_left,@forward_right,@reverse_left,@reverse_right,@parts); #arrays
my (%genecount,%genestart,%geneend,%genestrand,%antisensecount); #hashes

open(ART,$art_file)
	  or die "Error reading file: $!\n";
open(REF,$ref_file)
	  or die "Error reading file: $!\n";

#$tmp = `wc -l $art_file`;
#@parts = split(/\s/,$tmp);
#$genomesize = $parts[0]; #genomesize is now included in the tab annotation file header
#print "$genomesize\n";

#### initalize vectors ####
#check for header in annotation
$line = <REF>;
if(substr($line,0,1) eq "@"){
	#header found, first line should contain genome size
	@parts = split(/\|/,$line);
	$genomesize = $parts[2];
	chomp($genomesize)
#	print STDERR "Genomesize found in annotation header: $genomesize\n";
}

        $genomesize=6537648;
for($i=1; $i<=$genomesize;$i++){
	$forward_left[$i] = 0; 
	$forward_right[$i] = 0; 
        $reverse_left[$i] = 0; 
        $reverse_right[$i] = 0; 
}

#### parse ART ######
while(defined($line = <ART>)){
	#skip header
	if(substr($line,0,1) eq "#"){next}

	#read strand coverage	
	chomp($line);
	@parts = split(/\s/,$line);
	$i = $parts[0];
	$forward_left[$i] = $parts[1];
	$reverse_left[$i] = $parts[2];
	$forward_right[$i] = $parts[1];
	$reverse_right[$i] = $parts[2];
}

##### parse annotation #####
$n_genes = 0;
while(defined($line = <REF>)){
	chomp($line);
	#parse current gene
	if($type eq "gtf"){
		#use GTF format
		@parts = split(/\t/,$line);
#		$genetype      = ???, TODO
		$geneID    = $parts[1];
		$genestart{$geneID} = $parts[3];
		$geneend{$geneID}   = $parts[4];
		$genestrand{$geneID}= $parts[6];
		$genetxt   = $parts[8];
		@parts = split(/"/,$genetxt);

	}elsif($type eq "tab"){
		#skip header
		if(substr($line,0,1) eq "@"){next}

		#use tab separated format from pseudomonas.com
		@parts = split(/\t/,$line);
		if(defined($parts[8])){
			$geneID = "$parts[3],$parts[8]";
		}else{
			$geneID = $parts[3];
		}
		$genetype   = $parts[4];
		$genestart{$geneID}  = $parts[5];
		$geneend{$geneID}    = $parts[6];
		$genestrand{$geneID} = $parts[7];

	}else{
		print STDERR "No format description found!\n";
	}
	$n_genes++;

	#skip non-coding genes, if this option is set
	if($cds_only&&($genetype ne "CDS")){next}

	#check annotation for errors
	if(!( ($genestart{$geneID} =~ /^[0-9]+$/) & ($geneend{$geneID} =~ /^[0-9]+$/) )){
		print STDERR "!!! Annotation error for gene $geneID !!!\n";
	}

	#count reads
       $genecount{$geneID}=0;
	for($i = $genestart{$geneID}; $i <= $geneend{$geneID}; $i++){

		if($i > $genomesize){
			print STDERR "Out of range values for gene $geneID\n";
			last;
		}
		if($genestrand{$geneID} eq "+"){
			if(($forward_left[$i]>0) or ($forward_right[$i]>0)) { $genecount{$geneID} ++;}

		}elsif($genestrand{$geneID} eq "-"){
			if(($reverse_left[$i]>0) or ($reverse_right[$i]>0)) {$genecount{$geneID} ++;}

		}else{
			print STDERR "Hmmm, inconsistent strand information for gene $geneID ($genestart{$geneID}:$geneend{$geneID}).\n";
		}


	}
$genecount{$geneID}=100*$genecount{$geneID}/($geneend{$geneID}-$genestart{$geneID}+1)
}

foreach $geneID (sort keys %genecount){
		print "$geneID\t$genestart{$geneID}\t$geneend{$geneID}\t$genestrand{$geneID}\t$genecount{$geneID}\n";

}

close ART;
close REF;
