## load modules
use POSIX qw/floor/;


## read/determine variables
my $aln = $ARGV[0];
my @stopcodon = qw(TAA TAG TGA);
my $seqname;

my $melcnv=0;$melstp=0; my $melfrm = 0; my $melsub = 0;
my $heccnv=0;$hecstp=0; my $hecfrm = 0; my $hecsub = 0;
my $elecnv=0;$elestp=0; my $elefrm = 0; my $elesub = 0;
my $parcnv=0;$parstp=0; my $parfrm = 0; my $parsub = 0;


# look for stop codons
open(fas, $aln) or die "can't open $aln $!\n";
while (my $line = <fas>) {
    chomp($line);
    if ($line =~/^>/) { $seqname = $line}
    else {
	$line =~ tr/atgc/ATGC/;
	$line =~ s/-//gi;
	my $seqlen = length $line; my $frame = ($seqlen/3) - floor($seqlen/3);
	my $stp = 0;
        my @codons;
        push @codons, substr($line, 0, 3, "") while length($line);
	for (my $i=0; $i < ((scalar @codons)-1); $i++) {
	    if ($codons[$i] ~~ @stopcodon) { $stp++; }
	}
	if ($seqname =~ /hmel/) { $melcnv++; if ($stp > 0) {$melstp++; if ($frame != 1) {$melfrm++}; if ($frame == 1) {$melsub++} } }
	if ($seqname =~ /hhec/) { $heccnv++; if ($stp > 0) {$hecstp++; if ($frame != 1) {$hecfrm++}; if ($frame == 1) {$hecsub++} } }
	if ($seqname =~ /hele/) { $elecnv++; if ($stp > 0) {$elestp++; if ($frame != 1) {$elefrm++}; if ($frame == 1) {$elesub++} } }
	if ($seqname =~ /hpar/) { $parcnv++; if ($stp > 0) {$parstp++; if ($frame != 1) {$parfrm++}; if ($frame == 1) {$parsub++} } }
	#print "$seqname\t$stp\t$frame\n"; }
    }
}

#print "hmel.cnv\thhec.cnv\thele.cnv\thpar.cnv\t";
#print "hmel.stp\thhec.stp\thele.stp\thpar.stp\t";
#print "hmel.frm\thhec.frm\thele.frm\thpar.frm\t";
#print "hmel.sub\thhec.sub\thele.sub\thpar.sub\n";
print "$melcnv\t$heccnv\t$elecnv\t$parcnv\t";
print "$melstp\t$hecstp\t$elestp\t$parstp\t";
print "$melfrm\t$hecfrm\t$elefrm\t$parfrm\t";
print "$melsub\t$hecsub\t$elesub\t$parsub\n";

