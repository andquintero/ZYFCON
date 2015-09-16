#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./Gen_x_activities.pl ../Results/Blasted_genes.list ../Data/gene-short
#creates a matrix of genes x activities using the file crated by 1_Fetch_genes.pl and a list of genes
#out: ../Results/Gen_x_activities.matrix

	my %genes_and_ec = create_hashofhashes_second("< $ARGV[0]");
	my %ec_and_genes = create_hashofhashes_first("< $ARGV[0]");

	#print_hash_of_hashes(\%ec_and_genes);
	#print_hash_of_hashes(\%genes_and_ec);
	#print_hash_of_hashes(\%ec_and_ec);
	
	my @ec;
	foreach my $mk_ec_genes (sort keys %ec_and_genes) {
		push (@ec, $mk_ec_genes);
		
	}
	
	my @genes;
	open my $GENES, "< $ARGV[1]";
	while ( my $line = <$GENES> ) {
		chomp $line;
		push (@genes, $line);
	}
	
	#print @genes;
	#print "$genes[1]\n";


	open my $OUT, "> ../Results/Gen_x_activities.matrix";
	
	print $OUT "Genes/Activities\t"; 
	for my $i (0 .. $#ec-1){
		print $OUT "$ec[$i]\t"; 
		
	}
	print $OUT "$ec[$#ec]\n"; 
		
	for my $i (0 .. $#genes) {
		#my %hash = $genes_and_ec{$mk_genes_ec};
		my @activity;
		#print "$genes[$i]\n";
		
		for my $j (0 ..$#ec){
			if ( exists $genes_and_ec{$genes[$i]}{$ec[$j]}) {
				push (@activity, "1"); 
				print "$genes_and_ec{$genes[$i]}{$ec[$j]}\n";
			} else {
				push (@activity, "0");
				#print "$genes_and_ec{$genes[$i]}{$ec[$j]}\n";
				
			}
		}
		
		
		print $OUT "$genes[$i]\t"; 
		for my $j (0 .. $#activity-1){
			print $OUT "$activity[$j]\t"; 
			
		}
		print $OUT "$activity[$#activity]\n"; 
		
		
	}
	
	close $OUT;


exit 0;


















#________________________________Subroutines___________________________________________________________

#Creates a hash of hashes, Firts column MasterKeys:
sub create_hashofhashes_first {
	my ($FILE) = @_;	
	my %masterhash; #hash of hashes with MasterKeys
	
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	while ( my $line = <DATA> ) {
		chomp $line;
		my @columns = split (/\t/, $line);
		my %subhash;
		
		if ( exists $masterhash{$columns[0]} ) {
			%subhash = %{$masterhash{$columns[0]}};
		} else {
			%subhash = ();
		}
		
		$subhash{$_}++ for ($columns[1]);
		$masterhash{$columns[0]} = {%subhash}; # se hace una hash de hashes,
	}
	close DATA;
	return %masterhash;
}

#________________________________________________________________________________________________________

#Creates a hash of hashes, Second column MasterKeys:
sub create_hashofhashes_second {
	my ($FILE) = @_;	
	my %masterhash; #hash of hashes with MasterKeys
	
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	while ( my $line = <DATA> ) {
		chomp $line;
		my @columns = split (/\t/, $line);
		my %subhash;
		
		if ( exists $masterhash{$columns[1]} ) {
			%subhash = %{$masterhash{$columns[1]}};
		} else {
			%subhash = ();
		}
		
		$subhash{$_}++ for ($columns[0]);
		$masterhash{$columns[1]} = {%subhash}; # se hace una hash de hashes,
	}
	close DATA;
	return %masterhash;
}
