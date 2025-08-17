#!/usr/bin/perl -w 
use strict;

my $file = shift; 
unless ($file and  $file=~/bam$/) {
	print "please input bam file\n" ; 
	exit ; 
}
open IN, "samtools view -h -f 0x2  $file chrM |" ||die$! ; 
while (<IN>) {
	chomp;
	if (/^\@/) { 
		next if  /\@SQ\tSN:chr[\d+|X|Y|Un]/ ;  
		print "$_\n" ;
	}else {
		my @aa = split /\t/,$_; 
#		print "$_\n" ;  
		my $pos = $aa[3] ; 	  
#		print "$pos\t$_\n" ; 
#		if ($pos >=500 and $pos <= 16569-500 ) {
#			next if $aa[4] <20; 
#		} 
		print "$_\n" ; 
	}
}
close IN; 
