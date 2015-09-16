#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./2_DataBase_network_no_duplicates_nodot_plusEC.pl ../Results/Blasted_genes_plusGS.list ../../../Metabolic_network/Results/network.reactions_slim.filter_40


	my %genes_and_ec = create_hashofhashes_second("< $ARGV[0]");
	my %ec_and_genes = create_hashofhashes_first("< $ARGV[0]");
	my %ec_and_ec = create_hashofhashes_first("< $ARGV[1]");

	#print_hash_of_hashes(\%ec_and_genes);
	#print_hash_of_hashes(\%genes_and_ec);
	#print_hash_of_hashes(\%ec_and_ec);

	
	
	my %genes_and_genes = ();
	
	my $count = 0;	
	my $totalgenes= scalar keys %genes_and_ec;
	foreach my $mk_genes_ec (sort keys %genes_and_ec) {
		$count++;
		print "searching gen $count of $totalgenes\n";
		foreach my $k_genes_ec (keys %{$genes_and_ec{$mk_genes_ec}}) {
			if ( exists $ec_and_ec{$k_genes_ec} ) {
				foreach my $k_ec_ec (keys %{$ec_and_ec{$k_genes_ec}}) {
					if ( exists $ec_and_genes{$k_ec_ec} ) {
						foreach my $k_ec_genes (keys %{$ec_and_genes{$k_ec_ec}}) {
							#print $OUT "$mk_genes_ec\t$k_ec_genes\n"; 
							if (exists $genes_and_genes{$mk_genes_ec.";".$k_genes_ec}) {
								$genes_and_genes{$mk_genes_ec.";".$k_genes_ec}{$k_ec_genes.";".$k_ec_ec} = 1;
								#my %subhash = %{$genes_and_genes{$mk_genes_ec}};
								#$subhash{$k_ec_genes} = 1;
								###$subhash{$_}++ for $k_ec_genes;
								#$genes_and_genes{$mk_genes_ec} = {%subhash}; # se hace una hash de hashes,
								#print "$genes_and_genes{$mk_genes_ec}{$k_ec_genes}";
							} elsif (exists $genes_and_genes{$k_ec_genes.";".$k_ec_ec}) {
								$genes_and_genes{$k_ec_genes.";".$k_ec_ec}{$mk_genes_ec.";".$k_genes_ec} = 1;
								#my %subhash = %{$genes_and_genes{$k_ec_genes}};
								#$subhash{$mk_genes_ec} = 1;
								###$subhash{$_}++ for $mk_genes_ec;
								#$genes_and_genes{$k_ec_genes} = {%subhash}; # se hace una hash de hashes,
								
							} else {
								my %subhash = ();
								$subhash{$k_ec_genes.";".$k_ec_ec} = 1;
								#$subhash{$_}++ for $k_ec_genes;
								$genes_and_genes{$mk_genes_ec.";".$k_genes_ec} = {%subhash}; # se hace una hash de hashes,
							}
						}
					}
				}
			}
			
		}
	}

	print_hash_of_hashes(\%genes_and_genes);	

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
  
  open my $OUT, "> ../Results/Gen-Gen.edgelist.txt";
  
	foreach my $masterkey (sort keys %masterhash) {
		#print $OUT "$masterkey\t"; #imprimir master key
		foreach my $key (keys %{$masterhash{$masterkey}}) {
			#$masterkey =~ s/^(.{9}).+(;.+$)/$1$2/;
			#$key =~ s/^(.{9}).+(;.+$)/$1$2/;
			print $OUT "$masterkey\t$key\n"; #use this line commenting lines 132 and 133 if you want transcripts IDs
			#print $OUT "$key\t"; #prints only subkey
		}
		#print $OUT "\n"; 
	}
	
  close $OUT;
	
  
}
