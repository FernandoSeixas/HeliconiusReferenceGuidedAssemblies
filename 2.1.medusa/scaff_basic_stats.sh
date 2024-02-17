while getopts f:o: option; do
  case "${option}" in
    f) FILE=${OPTARG};;
    o) OUTP=${OPTARG};;
  esac
done

## calculate total length of assembly
echo "assembly total length: " > $OUTP
grep -v ">" $FILE | awk '{len=len+length($1)} END {print(len)}' >> $OUTP
echo "" >> $OUTP
echo "" >> $OUTP

## calculate the length of each scaffold and order
echo "scaffolds lengths (smallest > largerst):" >> $OUTP
grep -v ">" $FILE | awk '{print length($1)}'  | sort -n >> $OUTP

