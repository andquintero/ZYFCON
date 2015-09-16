#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./4_Weight_adjacency_list.pl ../Results/Gen-Gen.nodot.list ../Results/genes_in_adjacency_matrix.list ../../GCN/S.csv ../../GCN/gene_names21122.csv 
# 0file: list of gen-gen relationships generated by 2_DataBase_network_no_duplicates.pl
# 1file: ordered list of genes in adjacency matrix
# 2file: similarity matrix
# 3file: ordered list of genes in similarity matrix

#_____________________#Load Gen-Gen relationships into a hash of hashes____________________

	print "Reading Gen-Gen network\n";
	my %genes_and_genes = create_hashofhashes_first("< $ARGV[0]");

#_____________________#Read Gene names present in adjacency matrix_________________________
	
	print "Reading Gene names in adjacency matrix\n";
	my %adjacencynames ;
	my $count = 0;
	open ADJN, "< $ARGV[1]";
	while ( my $line = <ADJN> ) {
		chomp $line;
		#$line = '"' . $line . '"';
		$adjacencynames{$line} = $count;
		$count++;
		#print "$line\n";
		#print "$adjacencynames{$line}\n";
	}
	close ADJN;
	
#_____________________#Read Gene names present in similarity matrix________________________
	
	print "Reading Gene names in similarity matrix\n";
	my @similnames ;
	open SIMN, "< $ARGV[3]";
	while ( my $line = <SIMN> ) {
		chomp $line;
		push (@similnames, $line);
		#print "$similnames[$#similnames]\n";
	}
	close SIMN;

#_____________________#Get indexes of similarity matrix to print _________________________

	print "Calculating indexes of similarity matrix to print\n";	
	my @index; #array to store indexes to print
	for my $i (0 .. $#similnames) {
		#print "$adjacencynames{$similnames[$i]}\t$similnames[$i]\n";
		
		if ( exists $adjacencynames{$similnames[$i]} ) {
			push (@index, $i);
			#print "$adjacencynames{$similnames[$i]}\t$similnames[$i]\n";
			#print "$index[$#index]\t$#index\n";
		}
		#else {print "$similnames[$i]\n";} #prints genes in similarity matrix that don't appear in adjacency matrix 
	}	

#_____________________#Similarity matrix and prints only indexes________________________

	print "Beginning to parse similarity matrix\n";
	$count = -1; #if similarity matrix doesn't have headers chage to 0
	my $pres = 0;
	my $edges = 0;
	open SIMA, "< $ARGV[2]";
	while ( my $line = <SIMA> ) {
		if ($pres<=$#index) {
		if ($count == $index[$pres]) {
			chomp $line;
			my @fields = split /,/, $line;
			#print "$#fields\n";
			
			open my $OUT, ">> ../Results/Gen-Gen_similarity.list";
			for my $j (0 .. $#index) {
				#1: rowgene:$j 2: columngene:$count
				if (exists $genes_and_genes{$similnames[$index[$j]]}{$similnames[$count]} || $genes_and_genes{$similnames[$count]}{$similnames[$index[$j]]}) {
					if ($fields[$index[$j]] > 0){
						print $OUT "$similnames[$count]\t$similnames[$index[$j]]\n";
						#print $OUT "$similnames[$count]\t$similnames[$index[$j]]\t$fields[$index[$j]]\n";
						$edges++;
						
					}
					
				}
				#print "$fields[$index[$j]]\t";
			}
			#print $OUT "\n";
			close $OUT;
			
			$pres++;
			#print "$pres\n";
		}}
		$count++;
		print "Parsing similarity matrix line: $count\n";
	}
	close SIMA; 
	
	print "\n$edges total edges\n";
	

exit 0;




#________________________________Subroutines___________________________________________________________

#Creates a hash of hashes, Firts column MasterKeys:
sub create_hashofhashes_first {
	my ($FILE) = @_;	
	my %masterhash; #hash of hashes with MasterKeys
	my $count = 0;
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	while ( my $line = <DATA> ) {
		chomp $line;
		my @columns = split (/\t/, $line);
		my %subhash;
		
		if ( exists $masterhash{$columns[0]} ) {
			$masterhash{$columns[0]}{$columns[1]} = 1;
			#%subhash = %{$masterhash{$columns[0]}};
		} else {
			%subhash = ();
			$subhash{$columns[1]} = 1;
			$masterhash{$columns[0]} = {%subhash}; # se hace una hash de hashes,
		}
		$count++;
		#print "$count\n";
		#$subhash{$_}++ for ($columns[1]);
		#$masterhash{$columns[0]} = {%subhash}; # se hace una hash de hashes,
	}
	close DATA;
	return %masterhash;
}

#________________________________________________________________________________________________________

#Reads a file of gen-gen pairs and returns only unique gene names:
sub get_only_names {
	my ($FILE) = @_;	
	my @onlynames; #array with gene names
	
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	
	print "\nGetting Genes\n";
	while ( my $line = <DATA> ) {
		
			chomp $line;
			my @columns = split (/\t/, $line);
			push (@onlynames, @columns); 
			#print "$#onlynames\n";
		
	
	}
	close DATA;
	print "Removing duplicates\n";
	@onlynames = uniq(@onlynames);
	print "Total number of genes: $#onlynames\n";
	return @onlynames;
}

#remove duplicates
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
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

#________________________________________________________________________________________________________

#Prints an array of arrays

sub print_array_of_arrays {
  my $master_reference = shift;
  my @masterarray = @$master_reference;
  
	for my $i (0 .. $#masterarray) {
		for my $j (0 .. $#{$masterarray[$i]}) {
			print "$masterarray[$i][$j]\t"; #imprimir metabolitos
		}
		print "\n"; 
	}
  
}
