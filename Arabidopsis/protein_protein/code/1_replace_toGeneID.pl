#!/usr/bin/perl
use strict;
use warnings;

#use= mapping.txt
#

#script that read list of pp interactions and convert de to Gene IDs using file mapping

my %dict; #hash containing mapping

#load gene names and EC numbers of uniprot

open my $MAP, "< $ARGV[0]" or die "ERROR - unable to open $ARGV[0] mapping file\n";
print "Reading mapping file\n";
while ( my $line = <$MAP> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	if (length($columns[0])>0 && length($columns[1])>0) {
		$dict{$columns[0]} = $columns[1];	
	}
} #end while
close $MAP;

#replace for alternative name
my %mapped; #hash containing lines of filtered net
my $count = 0;
open my $PPI, "< $ARGV[1]" or die "ERROR - unable to open $ARGV[1] PPI file\n";
print "Reading $ARGV[1] PPI file\n";
while ( my $line = <$PPI> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	if ($#columns ==2) {
		if ( exists $dict{$columns[0]} ) { 
			$columns[0] = $dict{$columns[0]}; 
			$count++;
		}
		if ( exists $dict{$columns[0]} ) { 
			$columns[1] = $dict{$columns[1]}; 
			$count++;
		}
		if (length($columns[0])>0 && length($columns[1])>0 && length($columns[2])>0) {
			my $entry = "$columns[0]\t$columns[1]\t$columns[2]";			
			$mapped{$entry} = 1;
		}


	}
} #end while
close $PPI;

#print filtered PPI file
print "$count replaces\nPrinting $ARGV[1] mapped PPI file\n";
my $nameout  = "../data/mapped_ppi.txt";
open my $OUT, '>', "$nameout";
foreach my $key (sort keys %mapped) {
	print $OUT "$key\n";	
}
close $OUT;

exit;

