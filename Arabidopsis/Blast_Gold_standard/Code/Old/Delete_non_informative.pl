#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./Delete_non_informative.pl ../Results/Gen-Gen.matrix



	my @genes_and_genes = create_arrayofarrays("< $ARGV[0]");
	#my @gene_names = get_only_names("< $ARGV[1]");

	#print_hash_of_hashes(\%genes_and_genes);
	#print "@gene_names\n";
	#print "$#gene_names\n";
	#print $genes_and_genes[0];
	my @genes_and_genes2;
	my @to_delete;
	my @gene_names;
	
	my $count = 0;
	for my $i (0 .. $#genes_and_genes) {
		my @array;
		for my $j ($count .. $#genes_and_genes) {
			push (@array, $genes_and_genes[$j][$i]); 
			#print "@array\n";
			
		}
		$genes_and_genes2[$i] = [@array]; #array of arrays
		#print "@array\n";
		$count++;
	}
	
#_________________search for non informative genes____________
	
	LINE: for my $i (0 .. $#genes_and_genes) {
		print "searching $i of $#genes_and_genes\n";
		
		my $sum1 = 0;
		my $sum2 = 0;
		#print_array_of_arrays(\@genes_and_genes2);
		#print_array_of_arrays(\@genes_and_genes);
		for my $j (0 .. $#{$genes_and_genes[$i]}) {
			 $sum1 = $genes_and_genes[$i][$j] + $sum1;
		}
		for my $j (0 .. $#{$genes_and_genes2[$i]}) {
			 $sum2 = $genes_and_genes2[$i][$j] + $sum2;
		}
		
		if ( $sum2+$sum1 > 0 ){
			#print "$i No\n"; 
			
		} else {
			#print "$i Si\n";
			push (@to_delete, $i); 
			
		}
	}
	
#_________________Delete non informative genes____________	
	for my $i (reverse 0 .. $#to_delete) {
		
		print "deleting $i of $#to_delete\n";
		
		splice @genes_and_genes, $to_delete[$i], 1;
		splice @genes_and_genes2, $to_delete[$i], 1;
		splice @gene_names, $to_delete[$i], 1;
		
	#print_array_of_arrays(\@genes_and_genes);
		
		#delete non informative indexes in 1 matrix
		for my $k (0 .. $#genes_and_genes) {
			if ( exists $genes_and_genes[$k][$to_delete[$i]] ) {
				my @array = @{$genes_and_genes[$k]};
				splice @array, $to_delete[$i], 1;
				$genes_and_genes[$k] = [@array];
			}
		}
		
		#delete non informative indexes in 2 matrix
		for my $k (0 .. $#genes_and_genes) {
			if ( exists $genes_and_genes2[$k][$to_delete[$i]] ) {
				my @array = @{$genes_and_genes2[$k]};
				splice @array, $to_delete[$i], 1;
				$genes_and_genes2[$k] = [@array];
			}
		}
		
		
	}
	
	for my $i (0 .. $#gene_names) {
		my @array = @{$genes_and_genes[$i]};
		unshift @array, $gene_names[$i];
		$genes_and_genes[$i] = [@array];
		
	}
	
	#print "@to_delete\n";
	#print "@gene_names\n";
	#print_array_of_arrays(\@genes_and_genes2);
	#print_array_of_arrays(\@genes_and_genes);

	print "printing matrix\n";
	#print_array_of_arrays(\@adjacency_matrix);

	open my $OUT, ">> ../Results/Gen-Gen_only_informative.matrix";
  	for my $i (0 .. $#genes_and_genes) {
  		for my $j (0 .. $#{$genes_and_genes[$i]}) {
  			print $OUT "$genes_and_genes[$i][$j]\t"; 
  		}
  		print $OUT "\n"; 
  	}
		

	close $OUT;
	
	


	#print_array_of_arrays(\@genes_and_genes2);

exit 0;




#________________________________Subroutines___________________________________________________________

#Creates an array of arrays:
sub create_arrayofarrays {
	my ($FILE) = @_;	
	my @masterarray; #hash of hashes with MasterKeys
	
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	my $count = 0;
	while ( my $line = <DATA> ) {
		chomp $line;
		my @columns = split (/\t/, $line);
		push (@gene_names, $columns[0]);
		shift @columns;
		print "loading $count\n";
		$masterarray[$count] = [@columns]; #array of arrays
		#print "$#columns\n";
		#print "@columns\n";
		
		$count++;
	}
	close DATA;
	return @masterarray;
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

