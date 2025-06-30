SHAMROCK=$0
SHAMDIR=$(dirname $SHAMROCK)
TELOCIRAPTOR="python3 $SHAMDIR/../telociraptor/code/telociraptor.py"

GENBASE=$1
SEQIN=$GENBASE.fasta
#if [ ! -f "$SEQIN" ]; then
#  echo "Cannot find $SEQIN"; exit 1
#fi


# [ ] : Add MINCHROM setting?

# ~~~~~~ Shamrock example ~~~~~~~~
# Download from ENA
if [ ! -z "$2" ] && [ ! -f "$SEQIN" ]; then
  if [[ "$2" == G* ]]; then
    echo "Downloading $2 as assembly..."
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/$2?download=true&gzip=true" -O $SEQIN.gz
    unpigz $SEQIN.gz
  else
    echo "Downloading $2 as bioproject..."
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/links/study?accession=$2&result=sequence" -O $SEQIN
  fi
fi

if [ ! -f "$SEQIN" ]; then
  echo "Cannot find $SEQIN"; exit 1
fi


# Reformat names without pipes first for Telociraptor:
sed -E 's/ENA\|([A-Z0-9]+)/\1 ENA|\1/g' $SEQIN | sed -E 's/>([^:]+)chromosome: ([0-9]+)/>chr\2 \1chromosome: \2/g' | tee ${SEQIN/.fasta/.fna} | grep '>'

# Use Telociraptor to rename and reorder the sequences. 
$TELOCIRAPTOR -seqin ${SEQIN/.fasta/.fna} -chromsort -basefile ${SEQIN/.fasta/}


SEQOUT=${SEQIN/.fasta/}.tweak.fasta
ls -lrth $SEQOUT
if [ -f "$SEQOUT" ]; then
  echo "Ready for use as input: $SEQOUT"
fi

