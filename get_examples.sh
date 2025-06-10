SHAMROCK=$0
SHAMDIR=$(dirname $SHAMROCK)
TELOCIRAPTOR="python3 ../$SHAMDIR/../telociraptor/code/telociraptor.py"

if [ ! -d "example" ]; then
  mkdir -v example
fi
cd example

TESTV=$($TELOCIRAPTOR --version) && echo "[$(date)] Telociraptor: $TESTV"
if [ -z "$TESTV" ]; then
	echo "ERROR! Cannot find Telociraptor. Kill script or manually reformat downloaded sequences."
fi


# ~~~~~~ Shamrock example ~~~~~~~~
# Download from ENA
echo "[$(date)] Downloading Shamrock genome from ENA..."
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/links/study?accession=PRJEB62713&result=sequence" -O drTriDubi3.fasta

# Reformat names without pipes first for Telociraptor:
echo "[$(date)] Reformatting names for Telociraptor compatibility..."
sed -E 's/ENA\|([A-Z0-9]+)/\1 ENA|\1/g' drTriDubi3.fasta | sed -E 's/>([^:]+)chromosome: ([0-9]+)/>chr\2 \1chromosome: \2/g' > drTriDubi3.fna

# Use Telociraptor to rename and reorder the sequences. 
echo "[$(date)] Telociraptor sorting and reformatting names..."
$TELOCIRAPTOR -seqin drTriDubi3.fna -chromsort -basefile drTriDubi3

cd ..

SEQIN=example/drTriDubi3.tweak.fasta
echo "[$(date)] Ready for use as input: $SEQIN"


# ~~~~~~ Hemp-nettle example ~~~~~~~~
$SHAMDIR/chromformat.sh daGalTetr1 PRJEB69540
