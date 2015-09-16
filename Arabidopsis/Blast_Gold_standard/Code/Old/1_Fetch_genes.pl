#!/usr/bin/perl
    use strict;
    use warnings;
    use LWP::Simple;

#use: ./1_Fetch_genes.pl ../Data/eclist

	my $ecnumber; #ecnumber readed from eclist
	my @genename; #genes names list
	my $geneurl = 'http://rest.kegg.jp/get/'; #url to download genes associated with an EC number
	
	#read eclist
	open my $ECLIST, "< $ARGV[0]";
	while ( my $line = <$ECLIST> ) {
		chomp $line;
		$ecnumber = $line;
		print "\n$line\t";
	#___________________________#download EC number related genes___________________________________________
		my $url = "http://rest.kegg.jp/link/genes/$line";
		#my $geneurl = 'http://rest.kegg.jp/get/'; #reset url to download genes associated with an EC number
		#print "\n$url\n";
		my $ecgenes = get($url) or die 'Unable to get $url';
		@genename = split(/\n/, $ecgenes);
		#print "$genename[0] \n";
		#print "$#genename\n";
		
		for my $i (0 .. $#genename){ #keep only gene name
			my @line = split(/\t/, $genename[$i]);
			$genename[$i] = $line[1];
			#print "$genename[$i] \n";
			#$geneurl = $geneurl . '+' . $line[1];
		}

		print $#genename+1;
		print " sequences to download\n";
	#___________________________#download gene sequences______________________________________________________
		my $count = 0;
		my $step = 10;
		while ( $count <= $#genename ) {
			$geneurl = 'http://rest.kegg.jp/get/'; #reset url to download genes associated with an EC number
			$count = $count + 10;
			#print "$count\n";
			if ( $count > $#genename ) {
				$step = $#genename - ($count - 10);
				$count = $#genename;
				for my $i ($count-$step .. $count) {
					$geneurl = $geneurl . '+' . $genename[$i];		
				}
				print "retreaving $ecnumber sequence: $count\n";	
				gene_seqs_fasta($geneurl);
				last;
			}
			for my $i ($count-$step .. $count-1) {
				$geneurl = $geneurl . '+' . $genename[$i];	
			}
			print "retreaving $ecnumber sequence: $count\n";	
			gene_seqs_fasta($geneurl);
		}
	#________________________________BLAST______________________________________________________
	if ($#genename+1 > 0) {
	print "Starting Blast for $ecnumber";
	system("blastp -query ../Data/ec_genes_query.fa -db ../../../DataBases/Arabidopsis_datasets/Arabidopsis_representative_gene_prot -out ../Results/one -outfmt 7 -evalue 1e-5 -num_threads 7");
	system("rm ../Data/ec_genes_query.fa");
	my %blast_results;
	open my $BLASTR, "< ../Results/one";
	while ( my $line = <$BLASTR> ) {
		chomp $line;
		if ( substr($line, 0, 1) ne '#' ) {
			my @blastline = split(/\t/, $line);
			#$genename[$i] = $line[1];
			$blast_results{$blastline[1]} = $blastline[0];
			#print "$blastline[1]\n";		
		}

	}
	open my $OUT, ">> ../Results/Blasted_genes.list";
	foreach my $key (keys %blast_results){
		print $OUT "$blast_results{$key}\t$key\n";
	}
	close $OUT;

	}
	#________________________________END BLAST______________________________________________________
	
	} #end while
	
    exit 0;
	
	#________________________________Subroutines______________________________________________________

	#subrutin that downloads gene secuences associates to an ec number:
	sub gene_seqs_fasta {
		my ($sequrl) = @_;
		my $noparsed = get($sequrl) or die 'Unable to get $sequrl';
		my @linestoparse = split(/\n/, $noparsed);
		#print "\n$noparsed\n";
	
		open my $OUT, '>>', "../Data/ec_genes_query.fa";
	
		my $j = 0;
		ONESEQ: for my $i ($j .. $#linestoparse) {
		my @header; $header[0] = '>'; $header[1] = $ecnumber;
		#print "$i\n";
			if ( substr($linestoparse[$i], 0, 5) eq "ENTRY" ) {
				my @line = split(/\s+/, $linestoparse[$i]);
				$header[3] = $line[1];
				#print "@header \n $line[1]\n $i\n";
				#$geneurl = $geneurl . '+' . $line[1];
			
				for my $k ($i+1 .. $#linestoparse) {
					if ( substr($linestoparse[$k], 0, 8) eq "ORGANISM" ) {
						my @line = split(/\s+/, $linestoparse[$k]);
						$header[2] = $line[1];	
					
						for my $l ($k+1 .. $#linestoparse) {
							if ( substr($linestoparse[$l], 0, 5) eq "AASEQ" ) {
								print $OUT "$header[0]$header[1] $header[2]:$header[3]\n";
							
								for my $m ($l+1 .. $#linestoparse) {
									if ( substr($linestoparse[$m], 0, 4) eq '    ' ) {
									
										print $OUT split(/\s+/, $linestoparse[$m]); print $OUT "\n";
									} else {
										next ONESEQ;	
									}
								}
							}
						
						}
					}
				}
			
			}
		
		}
		close $OUT;
	}


