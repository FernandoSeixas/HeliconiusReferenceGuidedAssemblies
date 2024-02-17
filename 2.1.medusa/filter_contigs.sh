while getopts a:m:b:s: option; do
  case "${option}" in
    a) ASSEMBLY=${OPTARG};;
    m) MINLEN=${OPTARG};;
    b) OUTbig=${OPTARG};;
    s) OUTsma=${OPTARG};;
  esac
done


# subset contigs based on size 
awk -v min="$MINLEN" 'length($1) > min {printf "%s\n%s\n",f,$1} {f=$1}' $ASSEMBLY > $OUTbig
awk -v min="$MINLEN" 'length($1) < min && length($1) > 200 {printf "%s\n%s\n",f,$1} {f=$1}' $ASSEMBLY > $OUTsma

