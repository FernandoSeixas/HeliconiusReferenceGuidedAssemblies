while getopts i: option; do
  case "${option}" in
    i) NAME=${OPTARG};;
  esac
done

grep -v ">" $NAME.fasta | \
sed 's/A/0/g' | \
sed 's/T/0/g' | \
sed 's/G/0/g' | \
sed 's/C/0/g' | \
sed 's/N/1/g' | \
tr -s '0' | \
tr 0 \\n | awk 'BEGIN{RS="\n\n" ; ORS="\n";}{ print }' | \
awk '{print length($1)}' | sort -n | uniq -c | awk '{print $2,$1}' > $NAME.gapsize


