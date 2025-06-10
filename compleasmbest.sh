# This helper script converts Compleasm results into the best-hit table for homeologue assignment
GENBASE=$1
RESULTS=$2
N=$3
BESTN=1
if [ ! -z "$4" ]; then
  BESTN=$4
fi
echo "Extracting best $BESTN hits from $RESULTS for $GENBASE..."

TXT=$GENBASE.best.txt
if [ -f "$TXT" ]; then
  cp -v $TXT $TXT.bak
fi
for i in $(seq 1 $N); do
  CHR=$(printf "chr%02d" "$i")
  #echo chr$i "->" $CHR
  GENES=$(awk -v chr="$CHR" '$3 == chr {print $1;}' $RESULTS)
  
  
  awk -v genes="$GENES" '
    BEGIN {
      split(genes, glist, " ")
      for (i in glist) gene[glist[i]] = 1
    }
    gene[$1]
    ' "$RESULTS" | awk '{print $3;}' | sort | uniq -c | sort -n -r | awk -v chr=$CHR '$2 != chr {print $1, chr, $2;}' | head -n $BESTN | tee -a $TXT
done