## input parameters ==================================================
my $fasta = $ARGV[0];


## read data =========================================================
my $chrom;

open(fas, $fasta) or die "can't open $fas $!\n";
while (my $line = <fas>) {
    chomp($line);
    # get chromosome name and store it
    if ($line =~ /^\>/) {
	    my @id = split(" ", $line);
	    $chrom = $id[0];
	    $chrom =~ s/\>//g;
    }
    # split chromosomes into scaffolds
    else {
        my @scaffs = split /N+/, $line;
        my $last = scalar @scaffs;
        for (my $i=0; $i <= $last; $i++) {
    	    if (length $scaffs[$i] > 0) {
                    print ">", $chrom, "-Scaff_", $i+1, "\n";
                    print "$scaffs[$i]\n";
    	    }
        }
    }
}
