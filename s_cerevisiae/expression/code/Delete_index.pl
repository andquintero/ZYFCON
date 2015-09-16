#!/usr/bin/perl
    use strict;
    use warnings;

#use: ./Delete_index.pl ../Results/gene_names ../Results/mini_gene_names ../Results/index_to_delete ../Results/matrix ../Results/columns
# 1file= big list
# 2file= small list
# 3file= indexes
# 4file= matrix
# 5file= rows


	my %complete = create_hash("< $ARGV[0]");
	my %conserve = create_hash("< $ARGV[1]");
	#my @gene_names = get_only_names("< $ARGV[1]");
	my @ind ;
	#print_hash_of_hashes(\%genes_and_genes);
	
	open my $OUT, "> ../Results/index_to_delete";
	foreach my $con_keys (sort keys %conserve) {
		print $OUT "$complete{$con_keys}\n";
		push (@ind, $complete{$con_keys});
		#print "$#ind\n";
	} 
	

	close $OUT;
	
	my %index = create_hash("< $ARGV[2]");
	
#__________delete rows___________________________

open my $ROWS, "> ../Results/rows";
open my $MATRIX, "< $ARGV[3]";
my $count = 0;
while ( my $line = <$MATRIX> ) {
	if ( exists $index{$count} ) {
		print $ROWS "$line";
	}
	$count++;
}

close $ROWS;
close $MATRIX;

#__________delete columns___________________________

open my $FINAL, "> ../Results/columns";
open my $COLUMNS, "< $ARGV[4]";
print "delete columns\n";

print "$#ind\n";

while ( my $line = <$COLUMNS> ) {
	chomp $line;
	#print "$line\n";
	my @columns = split (/\t/, $line);
	#print "$#columns\n";
	#print "$columns[0]\n";
	
	for my $i ( 0 .. $#ind) {
		#if ( exists $columns[$ind[$i]-1] ) {
			print $FINAL "$columns[$ind[$i]-1]\t";
		#}
	}
	#print "$#columns\n";
	
	print $FINAL "\n";
	
}

close $COLUMNS;
	
	
	


	#print_array_of_arrays(\@genes_and_genes2);

exit 0;




#________________________________Subroutines___________________________________________________________

#Creates a hash:
sub create_hash {
	my ($FILE) = @_;	
	my %hash; #hash  with MasterKeys
	
	#open my $BLASTEDGENES, "< $ARGV[0]";
	open DATA, $FILE;
	my $count = 1;
	while ( my $line = <DATA> ) {
		chomp $line;
		$hash{$line} = $count;
				
		$count++;
	}
	close DATA;
	return %hash;
}




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
		#push (@gene_names, $columns[0]);
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


