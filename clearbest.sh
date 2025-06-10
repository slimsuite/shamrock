# Delete the files needed to redefine the best pairs of chromosomes and re-run
GENBASE=$1

rm -v kmc_$GENBASE/*.allo*
rm -v kmc_$GENBASE/*.alt.*

rm -v $GENBASE.0[34567]*.done
rm -v $GENBASE.parent.[12].fasta
rm -v $GENBASE.*.best.txt
rm -v $GENBASE.allo*.csv

mv -v $GENBASE.rawk.pdf $GENBASE.badbest.rawk.pdf 
mv -v $GENBASE.shamrock.pdf $GENBASE.badbest.shamrock.pdf
mv -v $GENBASE.parent.1.txt $GENBASE.badbest.parent.1.txt
mv -v $GENBASE.parent.2.txt $GENBASE.badbest.parent.2.txt

