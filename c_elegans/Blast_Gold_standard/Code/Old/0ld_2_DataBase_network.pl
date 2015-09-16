#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./2_DataBase_network.pl ../Results/Blasted_genes.list ../../Metabolic_network/Results/network.reactions_slim.filter_40


	my %genes_and_ec = create_hashofhashes_second("< $ARGV[0]");
	my %ec_and_genes = create_hashofhashes_first("< $ARGV[0]");
	my %ec_and_ec = create_hashofhashes_first("< $ARGV[1]");

	#print_hash_of_hashes(\%ec_and_genes);
	#print_hash_of_hashes(\%genes_and_ec);
	#print_hash_of_hashes(\%ec_and_ec);

	open my $OUT, ">> ../Results/Gen-Gen.list";
		
	foreach my $mk_genes_ec (sort keys %genes_and_ec) {
		foreach my $k_genes_ec (keys %{$genes_and_ec{$mk_genes_ec}}) {
			if ( exists $ec_and_ec{$k_genes_ec} ) {
				foreach my $k_ec_ec (keys %{$ec_and_ec{$k_genes_ec}}) {
					if ( exists $ec_and_genes{$k_ec_ec} ) {
						foreach my $k_ec_genes (keys %{$ec_and_genes{$k_ec_ec}}) {
							print $OUT "$mk_genes_ec\t$k_ec_genes\n"; 
						}
					}
				}
			}
			
		}
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

#________________________________________________________________________________________________________

#Prints a hash oh hashes

sub print_hash_of_hashes {
  my $master_reference = shift;
  my %masterhash = %$master_reference;
  
	foreach my $masterkey (sort keys %masterhash) {
		print "$masterkey\t"; #imprimir EC number
		foreach my $key (keys %{$masterhash{$masterkey}}) {
			print "$key\t"; #imprimir metabolitos
		}
		print "\n"; 
	}
  
}
