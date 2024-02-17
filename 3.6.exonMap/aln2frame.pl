## read/determine variables
my $aln = $ARGV[0];
my $sta = $ARGV[1];
my $exc = $ARGV[2];
my $frame = $ARGV[3];

my $cnt = 1;
my $stabase = 1;
my $endbase = 1;

# first pass alignment to determine firsr valid base (using melpomene as reference)
open(fas, $aln) or die "can't open $aln $!\n";
while (my $line = <fas>) {
    chomp($line);
    if ($cnt <= 1) {
        if ($line =~/^>/) { next;}
        else {
            my @bases = split("", $line); # split sequence by nucleotides
            for (my $i=0; $i < scalar @bases; $i++) {
		if ($bases[$i] eq "-") {next;}
		else {
		    $stabase = $i;
		    last;
		}
            }
            for (my $i=0; $i < scalar @bases; $i++) {
		if ($bases[$i] eq "-") {next;}
		else {$endbase = $i;}
            }
        }
    }
    else {next;}
    $cnt++;
}
close $fas;

open(fas, $aln) or die "can't open $aln $!\n";
while (my $line = <fas>) {
    chomp($line);
    if ($line =~/^>/) { 
	    print "$line\n"; 
	    }
    else {
	    my @bases = split("", $line);
	    my @seqbnd = @bases[ ($stabase+$frame-1) .. ($endbase-$exc) ];
	    #my @seqbnd = @bases[ ($stabase+$frame-1) .. ($stabase+$end-1) ];
	    print join("", @seqbnd),"\n";
    }
}

#print "$stabase\n";
#print "$endbase\n";
