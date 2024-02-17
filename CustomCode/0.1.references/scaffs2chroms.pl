## input parameters ========================================
my $scaffSizes = $ARGV[0];
my $ref = $ARGV[1];

## read data ===============================================
open(input, $scaffSizes) or die "can't open $scaffSizes $!\n";
while (my $line = <input>) {
  # read data
  chomp($line);
  my @cols = split("\t", $line);
  my $scaffold = $cols[0];
  my $scaSize  = $cols[1];
  # define chromosome group
  if ($ref eq "heradem") { 
	  my @scaNb = split("Herato", $scaffold);
	  my @scaDigits = split("", $scaNb[1]);
	  my $scaGroup = $scaDigits[0] . $scaDigits[1];
	  #my $scaOrder = $scaDigits[2] . $scaDigits[3] . $scaDigits[4] . $scaDigits[5];
	  print "Herato$scaGroup\t$scaffold\n";
  }
  if ($ref eq "hmelv25") { 
	  my @scaNb = split("Hmel2", $scaffold);
	  my @scaDigits = split("", $scaNb[1]);
	  my $scaGroup = $scaDigits[0] . $scaDigits[1];
	  if ($scaGroup == 00) { print "$scaffold\t$scaffold\n";}
	  if ($scaGroup != 00) { print "Hmel2$scaGroup\t$scaffold\n";}
	  #my $scaOrder = $scaDigits[2] . $scaDigits[3] . $scaDigits[4] . $scaDigits[5];
	  }
}
close $scaffSizes;

