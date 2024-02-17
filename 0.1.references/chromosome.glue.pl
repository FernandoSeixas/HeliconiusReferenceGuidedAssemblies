## libraries
use strict; # help
use warnings;
use  Scalar::Util qw(looks_like_number);


## input parameters 
my $species = $ARGV[0];
my $fasta   = $ARGV[1];
my $prefix  = $ARGV[2];
my $newprf  = $ARGV[3];

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
        my $chrname; 
        # unique chromosome identifier
        if ($species eq "heradem") {
            $chrname = $fields[0];
            $chrname =~ s/>$prefix//g ;                 # remove leading prefix
            my @chars = split("", $chrname);
            $chromID = $chars[0] . $chars[1];
        }
        if ($species eq "heralat") {
            $chrname = $fields[0];
            $chrname =~ s/>$prefix//g ;                 # remove leading prefix
            my @chars = split("_", $chrname);
            $chromID = $chars[0];
        }
        if ($species eq "sara") {
            $chrname = $fields[1];
            $chrname =~ s/$prefix//g ;                  # remove leading prefix
            my @chars = split("_", $chrname);
            $chromID = $chars[1];
        }
        if ($species eq "hcha_pr") {
            $chrname = $fields[0];
            $chrname =~ s/>$prefix//g ;                 # remove leading prefix
            my @chars = split("", $chrname);
            $chromID = $chars[0] . $chars[1]; 
        }
        if ($species eq "hcha_tx") {
            $chrname = $fields[0];
            $chrname =~ s/>$prefix//g ;                 # remove leading prefix
            my @chars = split("_", $chrname);
            $chromID = $chars[0];
        }
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
            $chromosomes{$chromID} .=  "N"x500; ## add 500 Ns
            $chromosomes{$chromID} .= $line; ## GLUE sequence
        } 
    }
}    


for my $k (sort { $a <=> $b } keys %chromosomes) {
    if (looks_like_number($k)) {
        print ">$newprf$k\n";
        print "$chromosomes{$k}\n";        
    }
}
