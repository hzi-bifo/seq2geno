#!/usr/bin/perl

use Getopt::Std;
use Statistics::Lite qw(:all);
use POSIX;

our($opt_p,$opt_r);
getopts('pr');

my $input_sam = @ARGV[$#ARGV];

if(!defined($opt_p)) { $opt_p = 0 } else { $opt_p = 1 }
if(!defined($opt_r)) { $opt_r = 0 } else { $opt_r = 1 }

#read input file
open(INPUT,$input_sam)
  or die "Error reading file: $!\n";

#parameters
my $min_mapqual = 20;

#preallocate variables
my $paired		= $opt_p;
my $r_mode		= $opt_r;	# write output in format that can be read by the genes_statistics.R script
my $n_dimer		= 0;
my $n_unmapped		= 0;
my $n_reads		= 0;
my $n_lowmapqual	= 0;
my ($line,$unmapped_perc, $mapped_perc, $highmapq_perc, $dimer_perc);
my (@parts, @reads, @mism, @inserts);
for(my $i = 0; $i <= $genomesize; $i++){
	$reads[$i] = 0;
}

#parse reads
while(defined($line = <INPUT>)){
	@parts = split(/\t/,$line);
	#skip header
	if(substr($parts[0],0,1) eq "@"){next}

	#count unmapped reads
	if($parts[1] & 4){$n_unmapped++; next}

	#count adapter dimer hits
	if($parts[2] eq 'RNAseq_PCR-adapter-dimer'){$n_dimer++; next};

	#count low mapping quality hits
	if($parts[4] < $min_mapqual){$n_lowmapqual++};

	#count mapped reads
	$n_mapped++;

	#count mismatches
	if(substr($parts[$#parts],0,2) eq "NM"){
		$mism[substr($parts[$#parts],5,2)]++;
	}

	#paired end options
	if($paired && ($parts[1] & 32)){ # first read in paired end mode
		push(@inserts,$parts[8])
	}
}
$n_reads = $n_dimer + $n_unmapped + $n_mapped;
$unmapped_perc = $n_unmapped/$n_reads * 100;
$mapped_perc   = $n_mapped/$n_reads * 100;
$dimer_perc    = $n_dimer/$n_reads * 100;
$highmapq_perc  = $mapped_perc - ($n_lowmapqual/$n_reads * 100);

if($r_mode){
	print "Total_reads:\t$n_reads\n";
	print "Mapped_reads:\t$n_mapped\n";
	print "HighMapQ_reads:\t".($n_mapped - $n_lowmapqual)."\n";
	if($paired){ # calculate the median insert size
		print "Median_insert:\t".median(@inserts)."\n";
	}
} else {
	print "Total reads: $n_reads \n";
	print "Mapped reads: \t\t$mapped_perc %\n";
	print "Unmapped/repeat reads:\t$unmapped_perc %\n";
	print "Adapter dimers: \t$dimer_perc %\n";
	print "Mapped reads with MQ >= $min_mapqual: $highmapq_perc %\n";
	if($paired){ # calculate the median insert size
		print "Median insert size:\t".median(@inserts)."\n";
	}
}
close(INPUT);
