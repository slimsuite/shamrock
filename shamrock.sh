##################################################################
### SHAMROCK: Splitting Homeologue Ancestors by Mapping  ~~~~~ ###
###             Repeats with Overlapping Common Kmers    ~~~~~ ###
### MAIN WORKFLOW EXECUTION SCRIPT                       ~~~~~ ###
### VERSION: 0.0.1                                       ~~~~~ ###
### LAST EDIT: 07/06/25                                  ~~~~~ ###
### AUTHORS: Richard Edwards 2025                        ~~~~~ ###
### CONTACT: https://github.com/slimsuite/shamrock      ~~~~~ ###
##################################################################

# This is the primary script for the SHAMROCK allotetraploid partitioning workflow.
# Usage: ./shamrock.sh <INPUT FASTA>

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial working version.

####################################### ::: TO DO ::: ##############################################
# [ ] : Separate out the preflight checks into preflight.exec so it can be run separately.
# [ ] : Concatenate the parents into two files.
# [ ] : Consider running Compleasm and ChromSyn on the two sets of parents.
# [ ] : Add a check that the sequences are actually found, i.e. the formatting is correct.

####################################### ::: SETUP ::: ##############################################

### ~~~~~~~ SHAMROCK code paths ~~~~~~~~ ##
SHAMROCK=$0
SHAMDIR=$(dirname $SHAMROCK)

### ~~~~~~~ Input fasta file ~~~~~~~~ ##
SEQIN=$1
if [ -z "$SEQIN" ]; then
	echo "Usage: ./shamrock.sh <INPUT FASTA> [<k>] [debug]"; exit 1
fi

if [ ! -f "$SEQIN" ]; then
	echo "Input sequence file not found '$SEQIN'"
	echo "Usage: ./shamrock.sh <INPUT FASTA> [<k>] [debug]"; exit 1
fi
#i# Output prefix
GENBASE=$(basename $SEQIN | awk -F '.' '{print $1;}')
LOG=$GENBASE.log
echo "[$(date)] Running SHAMROCK: $SHAMROCK" | tee -a $LOG
echo "[$(date)] Log output: $LOG" | tee -a $LOG

### ~~~~~~~~~~ Set variables ~~~~~~~~~~ ###
#i# Input fasta file
echo "[$(date)] Input sequences: $SEQIN" | tee -a $LOG
echo "[$(date)] Output prefix: $GENBASE" | tee -a $LOG
#i# Number of chromosomes
N=$(grep -c "^>chr" $SEQIN)
echo "[$(date)] Number of chromosomes: $N" | tee -a $LOG
if [ "$N" == "0" ]; then
  echo "No chromosomes (names chrXX) found in $SEQIN"; exit 1
fi
TEMPDIR=tmp_$GENBASE
KMCDIR=kmc_$GENBASE
#i# Set K and Debug
K=31
if [ ! -z "$2" ]; then
	if [ "$2" != "debug" ]; then
	  K=$2
	fi
fi
echo "[$(date)] kmer length: $K" | tee -a $LOG
#i# Debugging
DEBUG=FALSE
if [ "$2" == "debug" ] || [ "$3" == "debug" ]; then
  DEBUG=TRUE
fi
echo "[$(date)] Debug mode: $DEBUG" | tee -a $LOG
#i# Run Directory
RUNDIR=$(pwd)
echo "[$(date)] Run from directory: $RUNDIR" | tee -a $LOG

##################################### ::: PREFLIGHT CHECKS ::: #######################################
echo "[$(date)] Preflight checks of dependencies..."
#i# KMC
TESTV=$(kmc | head -n 1) && echo "[$(date)] KMC: $TESTV" | tee -a $LOG
if [ -z "$TESTV" ]; then
	echo "ERROR! Cannot run: kmc"; exit 1
fi

#i# Directory setups
if [ ! -d "$TEMPDIR" ]; then
  mkdir -v $TEMPDIR
fi
if [ ! -d "$KMCDIR" ]; then
  mkdir -v $KMCDIR
fi

####################################### ::: SHAMROCK ::: ##############################################

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 1: KMC kmer profiles per chromosome ~~~~ ###
## Goal: Pull out each chromosome into a file and generate kmer profiles.
## Rationale: Chromosome kmer profiles will be used to identify homeologues and then allokmers.
## Method: Pull out each chromosome into a file and run KMC.
## Key inputs: Genome Assembly ($SEQIN)
## Key outputs: 
## - Individual chromosome fasta files ($KMCDIR/$GENBASE.$CHR.fasta)
## - Individual chromosome KMC profiles.
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.01-ChromKMC
START="[$(date)] Step $STEP"
KEYOUT=$STEP.done
if [ ! -f "$KEYOUT" ]; then

	#i# Run Telociraptor on Hap1
	echo "[$(date)] Pull out each chromosome into a file and generate kmer profiles with KMC..."

  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    grep -A 1 ">$CHR" $SEQIN | tee $KMCDIR/$GENBASE.$CHR.fasta | grep ">"
    kmc -k$K -ci1 -fm $KMCDIR/$GENBASE.$CHR.fasta $KMCDIR/$GENBASE.$CHR.k$K $TEMPDIR
  done && echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output file
if [ ! -f "$KEYOUT" ]; then
	echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
	exit 1
fi

# Outputs:
# $KMCDIR/
# ├── $GENBASE.$CHR.k$K.kmc_pre
# ├── $GENBASE.$CHR.k$K.kmc_suf
if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 2: Generate pairwise intersects of all kmers. ~~~~ ###
## Goal: Generate pairwise intersects of all kmers.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.02-KmerIntersect
START="[$(date)] Step $STEP"
CSV=$GENBASE.k$K.csv
KEYOUT=$CSV
KEYOUTS="$KEYOUT $STEP.done"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Running $STEP ..." | tee -a $LOG
	echo chri,chrj,ktype,knum | tee $CSV
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    for j in $(seq 1 $N); do
      CHRJ=$(printf "chr%02d" "$j")
      kmc_tools simple $KMCDIR/$GENBASE.$CHR.k$K $KMCDIR/$GENBASE.$CHRJ.k$K intersect $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K
      kmc_dump $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K.txt
      KNUM=$(wc -l $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K.txt | awk '{print $1;}')
      echo "$CHR,$CHRJ,kmer,$KNUM" | tee -a $CSV
    done
  done && echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG


else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 3: Pair up chromosomes based on all chromosome kmers. ~~~~ ###
## Goal: Pair up chromosomes based on all chromosome kmers.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.03-PairChrom
START="[$(date)] Step $STEP"
KEYOUT=$GENBASE.k$K.best.txt
KEYOUTS="$KEYOUT $STEP.done"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Running $STEP ..." | tee -a $LOG
  # Pair up chromosomes based on all chromosome kmers.
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
  done
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    grep "^$CHR," $CSV | sed 's/,/ /g' | awk '$1 != $2 {print $4, $1, $2;}' | sort -n -r | head -n1 | tee -a $GENBASE.k$K.best.txt
  done
  awk '{print $2, $3;}' $GENBASE.k$K.best.txt | tee $GENBASE.k$K.best.tmp
  awk '{print $3, $2;}' $GENBASE.k$K.best.txt | tee -a $GENBASE.k$K.best.tmp
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    rm $KMCDIR/$GENBASE.$CHR.alt.fasta
    for CHRJ in $(grep "^$CHR" $GENBASE.k$K.best.tmp | sort | uniq | awk '{print $2;}'); do
      echo "$CHR alt -> $CHRJ"
      cat $KMCDIR/$GENBASE.$CHRJ.fasta >> $KMCDIR/$GENBASE.$CHR.alt.fasta
    done
  done
  wc $KMCDIR/$GENBASE.*.alt.fasta
  rm -v $GENBASE.k$K.best.tmp && echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 4: Generate allokmers. ~~~~ ###
## Goal: Generate allokmers by subtraction of chromosome pairs.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.04-AlloKmers
START="[$(date)] Step $STEP"
KEYOUT=$STEP.done
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Running $STEP ..." | tee -a $LOG
	# Generate allokmers by subtraction of chromosome pairs.
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    kmc -k$K -ci1 -fm $KMCDIR/$GENBASE.$CHR.alt.fasta $KMCDIR/$GENBASE.$CHR.alt.k$K $TEMPDIR
    kmc_tools simple $KMCDIR/$GENBASE.$CHR.k$K $KMCDIR/$GENBASE.$CHR.alt.k$K kmers_subtract $KMCDIR/$GENBASE.${CHR}.allok$K
  done &&	echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 5: Generate pairwise intersects of all allokmers. ~~~~ ###
## Goal: Generate pairwise intersects of all allokmers.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
CSV=$GENBASE.allok$K.csv
START="[$(date)] Step $STEP"
STEP=$GENBASE.05-AlloKmerIntersect
KEYOUT=$CSV
KEYOUTS="$KEYOUT $STEP.done"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Running $STEP ..." | tee -a $LOG
  # Generate pairwise intersects of all allokmers.
  echo chri,chrj,ktype,knum | tee $CSV
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    for j in $(seq 1 $N); do
      CHRJ=$(printf "chr%02d" "$j")
      kmc_tools simple $KMCDIR/$GENBASE.$CHR.allok$K $KMCDIR/$GENBASE.$CHRJ.allok$K intersect $KMCDIR/$GENBASE.${CHR}-${CHRJ}.allok$K
      kmc_dump $KMCDIR/$GENBASE.${CHR}-${CHRJ}.allok$K $KMCDIR/$GENBASE.${CHR}-${CHRJ}.allok$K.txt
      KNUM=$(wc -l $KMCDIR/$GENBASE.${CHR}-${CHRJ}.allok$K.txt | awk '{print $1;}')
      echo "$CHR,$CHRJ,allo,$KNUM" | tee -a $CSV
    done
  done && echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 6: Run R Script. ~~~~ ###
## Goal: Cluster by allokmers, generate graphics and split into parents.
## Approach: Run accompanying shamrock.R Rscript.
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.06-Rscript
START="[$(date)] Step $STEP"
KEYOUT=$GENBASE.shamrock.pdf
KEYOUTS="$KEYOUT $GENBASE.parent.1.txt $GENBASE.parent.2.txt $STEP.done"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Running Rscript ..." | tee -a $LOG
	Rscript $SHAMDIR/shamrock.R basefile=$GENBASE k=$K log=$LOG
	echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 7: Make parents. ~~~~ ###
## Goal: Compile the two parent fasta files from the Shamrock clustering.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=$GENBASE.07-Parents
START="[$(date)] Step $STEP"
KEYOUT=$STEP.done
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ]; then

	# Execute code for step
	echo "[$(date)] Generating parent fasta files..." | tee -a $LOG
	for P in 1 2; do
	  PFILE=$GENBASE.parent.$P.fasta
	  if [ -f "$PFILE" ]; then
	    rm $PFILE
	  fi
	  for CHR in $(cat $GENBASE.parent.$P.txt); do
	    cat $KMCDIR/$GENBASE.$CHR.fasta >> $PFILE
	  done
	  wc $PFILE
	done && echo -e "$START\n[$(date)] Step $STEP complete" | tee $STEP.done | tee -a $LOG

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi
# Check for key output files
for KEYOUT in $KEYOUTS; do
	if [ ! -f "$KEYOUT" ]; then
		echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
		exit 1
	fi
done

if [ "$RUNMODE" = "$STEP" ]; then
	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
	exit 0;
fi

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###




####################################### ::: FINISH ::: #############################################
echo "#---------------" | tee -a $LOG
echo "[$(date)]: Run complete" | tee -a $LOG
echo "#---------------" | tee -a $LOG
exit 0
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
