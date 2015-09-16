#!/usr/bin/perl
use strict;
use warnings;

#this script read a GS gene list with EC number annotation and a network in edgelist format and returns all the nodes that are not gold standards
#use= 
# ./8_unannotated_GenesID.pl ../Data/arab_GS_EC-gene_GO.txt ../Results/Gen-Gen.edgelist_GS-GO_filtered.txt 
# ./8_unannotated_GenesID.pl ../Data/arab_GS_EC-gene_KEGG.txt ../Results/Gen-Gen.edgelist_GS-KEGG_filtered.txt 

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
my %unannotated; #hash containing unannotated genes IDs
my $count = 0;
open my $BIGNET, "< $ARGV[1]" or die "ERROR - unable to open $ARGV[1] network file\n";
print "Reading $ARGV[1] network file\n";
while ( my $line = <$BIGNET> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	if ( exists $gold{$columns[0]} ) {
	} else {
		$unannotated{$columns[0]} = 1;		
	}
	if ( exists $gold{$columns[1]} ) {
	} else {
		$unannotated{$columns[1]} = 1;		
	}
} #end while
close $BIGNET;

#print filtered net
print "Printing $ARGV[1] unannotated GeneID list\n";
my $nameout  = "../Results/unannotated_GeneID_list.txt";
open my $OUT, '>', "$nameout";

foreach my $key (sort keys %unannotated) {
	print $OUT "$key\n";		
}
close $OUT;

exit;


