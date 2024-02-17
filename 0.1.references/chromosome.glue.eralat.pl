## libraries
use strict; # help
use warnings;

## input parameters 
my $fasta = $ARGV[0];
my $prefix=$ARGV[1];

## new variables
my %chromosomes;
my $chromID; 
my $mergeDE=0;    # 0 for no merge and 1 to merge;

## merge scaffolds in same chromosome
open('fa', '<', $fasta) or die "can't open $fasta $!\n";
while (my $line = <fa>) {
    chomp($line);
    if ($line =~ /^\>/) {
        ## determine chromosome
        my @fields = split(" ", $line);
        my $chrname = $fields[0];
        $chrname =~ s/>$prefix//g ;
	my @newname = split("_", $chrname);
        # unique chromosome identifier
	$chromID = $newname[0];
        # decide if in the same chromosome
        if ($chromID ne "00") {
            if (exists $chromosomes{$chromID}) {
                $mergeDE=1;
                } # MERGE next when reading the sequence
            else {
                $mergeDE=0; 
                $chromosomes{$chromID}="";
                } # do NOT merge
        }
        else {
            $mergeDE=0;
            $chromID=$chrname;
            $chromosomes{$chromID}="";
        }        
    }
    else {
        if ($mergeDE == 0) {
            $chromosomes{$chromID} .= $line; ## NEW sequence
        }
        if ($mergeDE == 1) {
            $chromosomes{$chromID} .=  "N"x100; ## add 100 Ns
            $chromosomes{$chromID} .= $line; ## GLUE sequence
        } 
    }
}    

for my $k (sort keys %chromosomes) {
    print ">$prefix$k\n";
    print "$chromosomes{$k}\n";
}
