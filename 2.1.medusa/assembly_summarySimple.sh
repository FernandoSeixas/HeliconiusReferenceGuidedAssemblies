while getopts i:n: option; do
  case "${option}" in
    i) FILE=${OPTARG};;
    n) SPPN=${OPTARG};;
  esac
done

NAME=${FILE/.fasta/}
NAME=${NAME/.fa/}

## ------------------------------------------------------------
## get scaffold lengths and missing data (Ns) per scaffold ----
#grep -v ">" $NAME.fasta | sed 's/A//g' | sed 's/T//g' | sed 's/G//g' | sed 's/C//g' | awk '{print length($1)}' > $NAME.scaff_missN.txt
#grep -v ">" $NAME.fasta | awk '{print length($1)}' > $NAME.scaff_sizes.txt
grep -v ">" $NAME.fa | sed 's/A//g' | sed 's/T//g' | sed 's/G//g' | sed 's/C//g' | awk '{print length($1)}' > $NAME.scaff_missN.txt
grep -v ">" $NAME.fa | awk '{print length($1)}' > $NAME.scaff_sizes.txt

# print basic stats
nscaffs=`wc -l $NAME.scaff_sizes.txt | cut -d' ' -f1`
tsize=`awk '{sum=sum+$1} END { printf "%s\t", sum }' $NAME.scaff_sizes.txt`
lsize=`sort -n $NAME.scaff_sizes.txt | tail -n 21 | awk '{sum=sum+$1} END {printf "%s\n", sum}'`
printf "%s\t%s\t%s\t%s\n" $NAME $nscaffs $tsize $lsize

## ------------------------------------------------------------
## in R 
module load R/3.5.1-fasrc01
export R_LIBS_USER=$HOME/apps/R_v3.5.1:$R_LIBS_USER
Rscript ~/code/heliconius_seixas/1.pseudo_references/2.1.medusa/fasta2percmissingSimple.R -g $NAME.scaff_missN.txt -s $NAME.scaff_sizes.txt -n $SPPN
