## input parameters ----------------------------------------
my $melgenome = $ARGV[0];
my $genetable = $ARGV[1];
my $buffer    = $ARGV[2];
my $outputfile= $ARGV[3];

## hash table of hmel assembly chr->seq ----------------------------------------
my %mel;
my $chrom;

open(fasta, $melgenome) or die "can't open $fasta $!\n";
while (my $line = <fasta>) {
	chomp($line);
	if ($line =~ /^\>/) {
		my @id = split(" ", $line);
		$chrom = $id[0];
		$chrom =~ s/\>//g;
		$mel{$chrom} = "";
	}
	else {
		$mel{$chrom} = $line;
	}
}
close $melgenome;
#foreach (sort keys %mel) { print "$_\n$mel{$_}\n"; }


## read table of genes and coordinates
print "reading genes coordinates...\n";
my %genchr;
my %gensta;
my %genend;
open(genes, $genetable) or die "can't open $genetable $!\n";
while (my $line = <genes>) {
	chomp($line);
	my @cols = split("\t", $line);
	my $genename = $cols[0];
	my $chr = $cols[1];
	my $sta = $cols[2]-1;
	my $end = $cols[3]-1;
	if (exists $genchr{$chr}) {
		$genchr{$chr} = $genchr{$chr} . "," . $genename;
		$gensta{$chr} = $gensta{$chr} . "," . $sta;
		$genend{$chr} = $genend{$chr} . "," . $end;
	}
	else {
		$genchr{$chr} = $genename;
		$gensta{$chr} = $sta;
		$genend{$chr} = $end;
	}
}

#foreach (sort keys %genchr) { print "$_\n$genchr{$_}\n"; }

## get genes sequences 
print "extracting gene sequences...\n";
open(my $fh, '>', $outputfile) or die "Could not open file '$outputfile' $!";
foreach (sort keys %genchr) {
	my $ch = $_;
	print $ch,"\n";
	# get chromosome sequence
	my $seq = $mel{$ch};
	my @sq = split("", $seq);
	# get sequences of genes in that chromosome	
	my @genes = split(",", $genchr{$ch});                       # get array of genes in chromosome
	my @sts = split(",", $gensta{$ch});
	my @ens = split(",", $genend{$ch});
	foreach (my $i = 0; $i < scalar @genes; $i++) {
		my $st = $sts[$i] - $buffer; if ($st < 0) {$st = 1;}
		my $en = $ens[$i] + $buffer;
		my @geneseq = @sq[$st .. $en];
		print $fh "\>$genes[$i]\n";
		print $fh join("", @geneseq),"\n";
	}
}




