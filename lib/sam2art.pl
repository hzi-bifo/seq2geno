#!/usr/bin/perl -w
# sam2art.pl
# AUTHOR: Andreas Dötsch
# LAST REVISED: Aug 2012
# 

## convert SAM (.sam) of single and mate pair reads to Artemis readable (.art) format

use strict;
use Getopt::Std;
use POSIX;

#variables

my $usage = "\n\nusage: $0 -s <strands> [-l -p -w <window>] <sam file>\n".
            "Convert SAM (.sam) to Artemis readable (.art) format.\n\n".
	    "-s <strands>\t1 for unstranded or single stranded, 2 for double stranded \n".
	    "-l\t\t'sinister' profile (5'-coverage, Filiatrault et al.). When -p is also used, both 5' and 3'- ends are simulateously detected for both plus and minus strand. Default is unset (output of read depth).\n".
	    "-d\t\t'dexter' profile. Reads only 3'-ends of reads (not to be confused with 3'-ends of fragments).\n".
	    "-p\t\tpaired end (mate pairs) mode.\n".
	    "-4\t\tprint values of mates independently (paired end mode only.\n".
	    "-h <threshold>\tmap 'holes' of coverage (coverage is less than <threshold>).\n".
	    "-q <threshold>\tminimum mapping quality (default: 0).\n".
	    "-f\t\tFlat mode. Prints a flat file for fasta database access. Readcounts are printed in 6-digit fields for subsequent positions througout the genome.\n".
	    "-c <chromosome>\tname of contig/chromosome to analyse. If the reference contains multiple contigs or chromosomes, this option should be used to analyse each contig independently.\n".
	    "-w <window>\tIf set, read counts will be calculated for windows of the specified size. 'Sinister' option will be set automatically for w > 1. Default is 1 (base-wise).\n\n";
our($opt_s,$opt_l,$opt_w,$opt_p,$opt_h,$opt_q,$opt_c,$opt_4,$opt_f,$opt_d);
getopts('s:lpw:q:c:4fd') or die $usage;
if (!defined($opt_s) ) {$opt_s = 2; print STDERR "Using default strand setting: 2-stranded.\n";}
if (!defined($opt_w) ) {$opt_w = 1}
if (!defined($opt_l) ) {$opt_l = 0}
if (!defined($opt_q) ) {$opt_q = 0}
if (!defined($opt_p) ) {$opt_p = 0;$opt_4 = 0}
if (!defined($opt_4) ) {$opt_4 = 0}		
if (!defined($opt_f) ) {$opt_f = 0} else {$opt_f = 1}
if (!defined($opt_c) ) {$opt_c = ""}
if (!defined($opt_d) ) {$opt_d = 0}

if (($#ARGV + 1) < 1) {die $usage;}

my $input_file	= $ARGV[$#ARGV];	# SAM file
my $strands	= $opt_s;		# stranded option
my $sinister	= $opt_l;		# sinister option (forces -4 in paired end mode)
my $dexter      = $opt_d;		# dexter option
my $window	= $opt_w;		# windowsize
my $paired	= $opt_p;		# paired end option
my $minmapq	= $opt_q;		# minimum mapping quality
my $pe_sepout   = $opt_4;		# independent output of mates (in paired end mode)
my $flatmode    = $opt_f;		# Flat mode. Prints a flat file for fasta database access. Readcounts are printed in 6-digit fields for subsequent positions througout the genome.
my $usecontig	= $opt_c;		# processed contig/chromosome

if(!$pe_sepout && $sinister && $paired) { 
	$pe_sepout = 1;
	print STDERR "-l and -p were selected, forcing -4.\n";
}
if(!$paired && $pe_sepout) {
	$pe_sepout = 0;
	print STDERR "-4 was selected without -p, ignoring -4.\n";
}
if($window > 1 && $sinister == 0){
	$sinister = 1;
	print STDERR "You selected a window size > 1. Sinister option has been set automatically.\n";
}

# variables
my (@plus_left, @minus_left, @plus_right, @minus_right, @parts);
my ($i, $line, $flag, $pos, $mapq, $maxpos, $readlength, $pos_tmp, $n_minus, $n_plus, $readcount, $mappedcount, $mappedpercent, $basecount, $tmp1, $tmp2, $contig);

#print $input_file."\n";

#!!! Overflow workaround !!!
#This trigger prevents out-of-genome hits by cutting after the largest mapped position.
#This may lead to a limited loss of information at the 3' end of the sequence
my $no_overflow = 1;
my $highest_pos_found = 0;

open(INPUT,$input_file)
	  or die "Error reading file $input_file: $!\n";

$maxpos = 8000000; $pos_tmp = 0; 
$readcount = 0; $mappedcount = 0;$basecount = 0;

while(defined($line = <INPUT>)){
	#skip header lines
	if(substr($line,0,1) eq '@'){ next }
	
	#parse read
	@parts	= split(/\t/,$line);
	$flag	= $parts[1];
	$contig	= $parts[2];
	$pos	= $parts[3];
	$mapq	= $parts[4];
	$readlength = length($parts[9]);
#	$basecount += $readlength;

	$readcount++;

	#if several contigs/chromosomes are present, skip if the wrong one was hit
	if( ($usecontig ne "") && ($usecontig ne $contig) ){next}

	#skip read if mapping quality is below threshold (i.e. unmapped reads or non-unique hits)
	if($mapq < $minmapq){next}
	$mappedcount++;

	#count strand hits, paired end mode
	if($paired){
		if($sinister){
			if($window == 1){
				if($flag == 99)  {$plus_left[$pos]++} 					# left (5') read of a plus strand fragment
				if($flag == 147) {$plus_right[$pos+$readlength-1]++}			# right (3') read of a plus strand fragment
				if($flag == 163) {$minus_right[$pos]++}					# right (3') read of a minus strand fragment
				if($flag == 83)  {$minus_left[$pos+$readlength-1]++}			# left (5') read of a minus strand fragment
			} else {
				if($flag == 99)  {$plus_left[ceil($pos/$window)]++}			# left (5') read of a plus strand fragment
				if($flag == 147) {$plus_right[ceil(($pos+$readlength-1)/$window)]++}	# right (3') read of a plus strand fragment
				if($flag == 163) {$minus_right[ceil($pos/$window)]++}			# right (3') read of a minus strand fragment
				if($flag == 83)  {$minus_left[ceil(($pos+$readlength-1)/$window)]++}	# left (5') read of a minus strand fragment
			}
		} elsif($dexter){
			if($window == 1){
				if($flag == 99)  {$plus_left[$pos+$readlength-1]++} 			# left (5') read of a plus strand fragment
				if($flag == 147) {$plus_right[$pos]++}					# right (3') read of a plus strand fragment
				if($flag == 163) {$minus_right[$pos+$readlength-1]++}			# right (3') read of a minus strand fragment
				if($flag == 83)  {$minus_left[$pos]++}					# left (5') read of a minus strand fragment
			} else {
				if($flag == 99)  {$plus_left[ceil(($pos+$readlength-1)/$window)]++}	# left (5') read of a plus strand fragment
				if($flag == 147) {$plus_right[ceil($pos/$window)]++}			# right (3') read of a plus strand fragment
				if($flag == 163) {$minus_right[ceil(($pos+$readlength-1)/$window)]++}	# right (3') read of a minus strand fragment
				if($flag == 83)  {$minus_left[ceil($pos/$window)]++}			# left (5') read of a minus strand fragment
			}
		} else {
			if($flag == 99){								# left (5') read of a plus strand fragment
				for($i = $pos; $i < $pos + $readlength; $i++){ $plus_left[$i]++   }
				$basecount += $readlength;
			}
			if($flag == 147) {								# right (3') read of a plus strand fragment
				for($i = $pos; $i < $pos + $readlength; $i++){ $plus_right[$i]++  }
				$basecount += $readlength;
			}
			if($flag == 163) {								# right (3') read of a minus strand fragment
				for($i = $pos; $i < $pos + $readlength; $i++){ $minus_right[$i]++ }
				$basecount += $readlength;
			}
			if($flag == 83) {								# left (5') read of a minus strand fragment
				for($i = $pos; $i < $pos + $readlength; $i++){ $minus_left[$i]++  }
				$basecount += $readlength;
			}
		}
	} else {	
	#count strand hits, single end mode
		if($sinister){
			if($window == 1){
				if($flag == 0) {$plus_left[$pos]++} #plus strand
				if($flag == 16) {$minus_left[$pos + $readlength - 1]++} #minus strand
			} else {
				if($flag == 0) {$plus_left[ceil($pos/$window)]++} #plus strand
				if($flag == 16) {$minus_left[ceil(($pos + $readlength - 1)/$window)]++} #minus strand
			}
		} elsif($dexter){
			if($window == 1){
				if($flag == 0) 	{$plus_left[$pos + $readlength - 1]++}  #plus strand
				if($flag == 16)	{$minus_left[$pos]++}	#minus strand
			} else {
				if($flag == 0) {$plus_left[ceil(($pos + $readlength - 1)/$window)]++} #plus strand
				if($flag == 16) {$minus_left[ceil($pos/$window)]++} #minus strand
			}
		} else {
			if($flag == 0) { #plus strand
				for($i = 0; $i < $readlength; $i++){
					$plus_left[$i+$pos]++;
				}
			}
			if($flag == 16) { #minus strand
				for($i = 0; $i < $readlength; $i++){
					$minus_left[$i+$pos]++;
				}
			}
		}
	}
	if($highest_pos_found < $pos){$highest_pos_found = $pos}
}
close(INPUT);

$mappedpercent = $mappedcount/$readcount * 100;

#print header
if($flatmode){ print "      " }
elsif($pe_sepout) {
	print "# base forward_left reverse_left forward_right reverse_right  # high quality (MQ >= $minmapq): $mappedcount/$readcount ($mappedpercent %)\n"; 	
	print "# colour 255:0:0 0:255:0 255:0:255 0:255:255\n";
}elsif($strands == 1) {
	print "# base coverage  # high quality (MQ >= $minmapq): $mappedcount/$readcount ($mappedpercent %)\n"; 	
	print "# colour 0:0:0\n";
}elsif($strands == 2) {
	print "# base forward reverse # high quality (MQ >= $minmapq): $mappedcount/$readcount ($mappedpercent %)\n"; 	
	print "# colour 255:0:0 0:255:0\n";
}

#workaround for out-of-genome hits
if($no_overflow > 0){$maxpos = $highest_pos_found}

#print output
for($i = 1; $i <= $maxpos; $i++){
	if($flatmode) {
		if(!(defined($plus_left[$i])))  { $plus_left[$i] = 0  }
		if(!(defined($minus_left[$i]))) { $minus_left[$i] = 0 }
		if(!(defined($plus_right[$i])))  { $plus_right[$i] = 0  }
		if(!(defined($minus_right[$i]))) { $minus_right[$i] = 0 }
		print sprintf("%06d",$plus_left[$i] + $minus_left[$i] + $plus_right[$i] + $minus_right[$i]);
#         print "\n";
	} else {
		# skip empty positions
		if(!(defined($plus_left[$i])))  { $plus_left[$i] = 0  }
		if(!(defined($minus_left[$i]))) { $minus_left[$i] = 0 }
		if($paired) {
			if(!(defined($plus_right[$i])))  { $plus_right[$i] = 0  }
			if(!(defined($minus_right[$i]))) { $minus_right[$i] = 0 }
			if(($plus_left[$i] + $minus_left[$i] + $plus_right[$i] + $minus_right[$i]) == 0) { next }
		} elsif(($plus_left[$i] + $minus_left[$i]) == 0) { next }
	
		# determine position in windowed mode
		$pos = floor($i*$window-($window-1)/2);
	
		# print values
		if($paired) {
			if($pe_sepout){
				print "$pos\t$plus_left[$i]\t$minus_left[$i]\t$plus_right[$i]\t$minus_right[$i]\n";
			} else {
				$tmp1 = $plus_left[$i] + $plus_right[$i];
				$tmp2 = $minus_left[$i] + $minus_right[$i];
				print "$pos\t$tmp1\t$tmp2\n";
			}
		} elsif($strands == 1) {
			print "$pos\t".($plus_left[$i] + $minus_left[$i])."\n";
		} elsif($strands == 2) {
			print "$pos\t$plus_left[$i]\t$minus_left[$i]\n";
		}
	}
}
