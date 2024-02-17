my $fasta = $ARGV[0];
my $chr = $ARGV[1];
my $sta = $ARGV[2];
my $end = $ARGV[3];
my $str = $ARGV[4];

if ($str eq '') { $str = "+";} # assumes positive strand if no argument given
if ($str eq "+") {
	$sta = $sta - 1; 
	$end = $end - 1;
	if ($sta < 0) {$sta = 0;} 
}
if ($str eq "-") {
	$sta = $sta - 1;
	$end = $end - 1;
	if ($sta < 0) {$sta = 0;} 
}

my $currentChrom;

open(fas, $fasta) or die "can't open $fasta $!\n";
while (my $line = <fas>) {
	chomp($line);
	if ($line =~/^>/) {
		$line =~ s/\>//g;
		my @cols = split(" ", $line);
		$currentChrom = $cols[0];
		}
	else {
		if ($currentChrom eq $chr) {
			my @cols = split("", $line);
			my @seq = @cols[ $sta .. $end];
			my $origin_seq = join("", @seq);
			if ($str eq "-") {
				my $revcomp = reverse $origin_seq;
				$revcomp =~ tr/ATGCatgc/TACGtacg/;
				$origin_seq = $revcomp;
			}
			#print ">$chr.$sta.$end\n";i
			print $origin_seq, "\n";
		}
		else { next; }
	}
}
