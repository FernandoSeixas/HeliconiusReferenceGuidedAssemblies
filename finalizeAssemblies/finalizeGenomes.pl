## input parameters ===================
my $species    = $ARGV[0];
my $fasta      = $ARGV[1];
my $senseTable = $ARGV[2];
my $sc2chTable = $ARGV[3];
my $ref        = $ARGV[4];

## hash table of scaff 2 chromosome =================
my %assocHash;
my $sppColumn;
my $lastCol;
open(association, $sc2chTable) or die "can't open $association $!\n";
while (my $line = <association>) {
    chomp($line);
    my @elements = split('\t', $line);
    # determine the column number of the species we're interested in
    if ($line =~ /^chromosome/) {
    	$lastCol = scalar @elements;
        for (my $i=1; $i < $lastCol; $i++) { if ($elements[$i] eq $species) { $sppColumn=$i;} }
    }
    # make association between chromosome names and scaffolds
    else {
        my $chrom = $elements[0];
        my $scaff = $elements[$sppColumn];
        $assocHash{$scaff} = $chrom;
    }
}
close $sc2chTable;
#foreach (sort keys %assocHash) { print "$_\t$assocHash{$_}\n"; }

## hash table of chromosome sense ===================
my %senseHash;
open(orientation, $senseTable) or die "can't open $orientation $!\n";
while (my $line = <orientation>) {
    chomp $line;
    # skip header line
    if ($line =~ /^chromosome/) { next; }
    # determine orientation of chromosome
    else {
    	my @elements = split('\t', $line);
        my $chrom = $elements[0];
        my $orient= $elements[$sppColumn];
        $senseHash{$chrom} = $orient;
    }
}
close $senseTable;
#foreach (sort keys %senseHash) { print "$_\t$senseHash{$_}\n"; }


## read genome file, change scaffold to chromosome names and reverse complement when needed ===================
my $seqHeader;
my $count = 100;
my $state;

my %finalGenome;

open(genome, $fasta) or die "can't open $genome $!\n";
while (my $line = <genome>) {
    chomp $line;
    # header line 
    if ($line =~ /^\>/) {
        $seqHeader = $line;
        $seqHeader =~ s/\>//g;
	# if it's a chromosome
        if (exists $assocHash{$seqHeader}) {
	    $newName = $assocHash{$seqHeader}; # get reference chromosome name
	    # get orientation
	    my $seqSense = $senseHash{$newName};
	    #print "$seqHeader\t$newName\t$seqSense\n";
	    if ($seqSense eq "+") {$state = "plus";}
	    if ($seqSense eq "-") {$state = "minus";}
            # update scaffold/chromosome name 
            $newName =~ s/$ref/00/g;
            $seqHeader = $newName;
            #print $seqHeader,"\n";
	    $finalGenome{$seqHeader} = "";

        }
	# if not a chromosome
	else {
	    $newName = "0" . $count;
            $seqHeader = $newName;
            #print $seqHeader,"\n";
	    $finalGenome{$seqHeader} = "";
	    $state="unknown";
	    $count++;
	}
    }
    else {
	if ($state eq "minus") {
	    #print "reverse\n"; 
	    my $revcomp = reverse $line;
	    $revcomp =~ tr/ATGCatgc/TACGtacg/;
	    #print $revcomp,"\n";
	    $finalGenome{$seqHeader} = $finalGenome{$seqHeader} . $revcomp;

	}
	else { 
	    #print $line,"\n";
	    $finalGenome{$seqHeader} = $finalGenome{$seqHeader} . $line;
	}
    }
}

foreach (sort {$a <=> $b} (keys %finalGenome) ) { print ">$species$_\n", "$finalGenome{$_}\n"; }

