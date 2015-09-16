#!/usr/bin/perl
use strict;
use warnings;

#use= 
# ./7_filter_bignet_GS.pl ../Data/arab_GS_EC-gene_GO.txt ../Results/Gen-Gen.edgelist_GeneID.txt 
# ./7_filter_bignet_GS.pl ../Data/arab_GS_EC-gene_KEGG.txt ../Results/Gen-Gen.edgelist_GeneID.txt 

my %gold; #hash containing names and EC numbers of the gold standard net

#load gene names and EC numbers of Gold Standards

open my $GOLDS, "< $ARGV[0]" or die "ERROR - unable to open $ARGV[0] Gold Standard file\n";
print "Reading Gold Standards\n";
while ( my $line = <$GOLDS> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	$gold{$columns[1].";".$columns[0]} = 1;
} #end while
close $GOLDS;

#filter for bignet, only edges with one GS are keeped
my @filterednet; #array containing lines of filtered net
my $count = 0;
open my $BIGNET, "< $ARGV[1]" or die "ERROR - unable to open $ARGV[1] network file\n";
print "Reading $ARGV[1] network file\n";
while ( my $line = <$BIGNET> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	if ( exists $gold{$columns[0]} || exists $gold{$columns[1]} ) {
		$filterednet[$count] = $line;
		$count++;
	}
} #end while
close $BIGNET;

#print filtered net
print "Printing $ARGV[1] filtered network file\n";
my $nameout  = "../Results/Gen-Gen_GS_filtered.txt";
open my $OUT, '>', "$nameout";
for my $i (0 .. $#filterednet) {
	print $OUT "$filterednet[$i]\n";	
}
close $OUT;

exit;


