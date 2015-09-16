#!/usr/bin/perl
use strict;
use warnings;

#use= 
#./1_filter_blast_results_dictionary.pl ../Data/uniprot_EC_dictionary.txt ../Results/blasted_query_gene.txt 

#script that read a formated blast result (by the BLAST.R script) and the ec number geneID of uniprot equivlaences and print the blasted genes list

my %gold; #hash containing names and EC numbers of the uniprot genes dictionary

#load gene names and EC numbers of uniprot

open my $GOLDS, "< $ARGV[0]" or die "ERROR - unable to open $ARGV[0] Gold Standard file\n";
print "Reading Uniprot Dictionary\n";
while ( my $line = <$GOLDS> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	$gold{$columns[0]} = $columns[1];
} #end while
close $GOLDS;

#filter for blastp result, query name is replaced by EC number
#my @filteredblast; #array containing lines of filtered net
my %filteredblast; #hash containing lines of filtered net
my $count = 0;
open my $BIGBLAST, "< $ARGV[1]" or die "ERROR - unable to open $ARGV[1] blast result file\n";
print "Reading $ARGV[1] blast result file\n";
while ( my $line = <$BIGBLAST> ) {
	chomp $line;
	my @columns = split (/\t/, $line);
	my $entry = "$gold{$columns[1]}\t$columns[0]";
	$filteredblast{$entry} = 1;
	#$filteredblast[$count] = $entry;
	$count++;
} #end while
close $BIGBLAST;

#print filtered net
print "Printing $ARGV[1] filtered blast file\n";
my $nameout  = "../Results/Blasted_genes.list";
open my $OUT, '>', "$nameout";
#for my $i (0 .. $#filteredblast) {
#	print $OUT "$filteredblast[$i]\n";	
#}
foreach my $key (sort keys %filteredblast) {
	print $OUT "$key\n";	
}
close $OUT;

exit;

