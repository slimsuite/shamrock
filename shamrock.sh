##################################################################
### SHAMROCK: Separating Homeologue Ancestors by Mapping ~~~~~ ###
###             Repeats with Overlapping Common Kmers    ~~~~~ ###
### MAIN WORKFLOW EXECUTION SCRIPT                       ~~~~~ ###
### VERSION: 0.4.0                                       ~~~~~ ###
### LAST EDIT: 11/06/25                                  ~~~~~ ###
### AUTHORS: Richard Edwards 2025                        ~~~~~ ###
### CONTACT: https://github.com/slimsuite/shamrock      ~~~~~ ###
##################################################################

# This is the primary script for the SHAMROCK allotetraploid partitioning workflow.
# Usage: ./shamrock.sh <INPUT FASTA> [<RUNMODE>]

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial working version.
# v0.2.0 : Added capacity to generate and use $GENBASE.best.txt.
# v0.3.0 : Added simple shamrock.config file and compleasm homeologue assignment.
# v0.3.1 : Modified the R output and added clustering of higher ploidies.
# v0.4.0 : Improved control and reporting of running steps. Renamed "parents" as "subgenomes".
VERSION=v0.4.0

####################################### ::: TO DO ::: ##############################################
# [?] : Separate out the preflight checks into preflight.exec so it can be run separately.
# [Y] : Concatenate the parents into two files.
# [ ] : Consider running Compleasm and ChromSyn on the two sets of parents.
# [ ] : Add a check that the sequences are actually found, i.e. the formatting is correct.
# [Y] : Make it easier to use an existing *.best.txt file.
# [Y] : Add a compleasm-based pairing of chromosomes to replace the kmer-based approach.
# [Y] : Add a simple shamrock.config file to easily over-ride some of the default settings without argv handling.
# [ ] : Add option to look at higher ploidy. (Can do this for homeologues)
# [ ] : Add option to pre-screen for polyploidy using compleasm and a duplication threshold. (Move compleasm step higher.)

####################################### ::: SETUP ::: ##############################################

### ~~~~~~~ SHAMROCK code paths ~~~~~~~~ ##
SHAMROCK=$0
SHAMDIR=$(dirname $SHAMROCK)

### ~~~~~~~ Input fasta file ~~~~~~~~ ##
SEQIN=$1
if [ -z "$SEQIN" ]; then
	echo "Usage: ./shamrock.sh <INPUT FASTA> [<RUNMODE>]"; exit 1
fi
if [ "$SEQIN" == "--version" ]; then
  echo "shamrock $VERSION"; exit 0
fi
if [ "$SEQIN" == "--help" ]; then
	echo "Usage: ./shamrock.sh <INPUT FASTA> [<RUNMODE>]"
	echo "Settings can be overridden using shamrock.config"
  echo "Please see https://github.com/slimsuite/shamrock for help with usage."; exit 0
fi
if [ ! -f "$SEQIN" ]; then
	echo "Input sequence file not found '$SEQIN'"
	echo "Usage: ./shamrock.sh <INPUT FASTA>"; exit 1
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
  echo "[$(date)] No chromosomes (names chrXX) found in $SEQIN"; exit 1
fi
TEMPDIR=tmp_$GENBASE
KMCDIR=kmc_$GENBASE
#i# Set Run Mode
RUNMODE=completion
if [ ! -z "$2" ]; then
  RUNMODE=$2
fi
echo "[$(date)] Running until step: $RUNMODE" | tee -a $LOG

### ~~~~~~~~~~ Run settings ~~~~~~~~~~ ###
#i# Set K
K=31
#i# Best homeologue matching strategy (compleasm/kmer/<FILE>)
BEST=compleasm
if [ -f $GENBASE.best.txt ]; then
  echo "[$(date)] Using $GENBASE.best.txt as paired chromosomes ..." | tee -a $LOG
  BEST=$GENBASE.best.txt
fi
#i# Compleasm lineage
LINEAGE=embryophyta
#i# Number of Best Homeologues
BESTN=1
#i# Debugging
DEBUG=FALSE
#i# Dev settings
DEV=FALSE
#i# Number of threads for compleasm etc.
THREADS=$(nproc)
#i# Config override
if [ -f "shamrock.config" ]; then
	echo "[$(date)] Loading settings from shamrock.config" | tee -a $LOG
	source shamrock.config
fi
if [[ ! "$THREADS" =~ ^-?[0-9]+$ ]]; then
  echo "[$(date)] Thread count not recognised - setting to 16: $THREADS" | tee -a $LOG
  THREADS=16
fi
echo "[$(date)] kmer length: $K" | tee -a $LOG
echo "[$(date)] Best Homeologue identification strategy: $BEST" | tee -a $LOG
if [ $BEST == "compleasm" ]; then
  echo "[$(date)] Compleasm lineage: $LINEAGE" | tee -a $LOG
fi
echo "[$(date)] Number of threads: $THREADS" | tee -a $LOG
echo "[$(date)] Number of best hits for homeologues: $BESTN" | tee -a $LOG
echo "[$(date)] Debug mode: $DEBUG" | tee -a $LOG
echo "[$(date)] Dev mode: $DEV" | tee -a $LOG

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
STEPN=0

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
STEP=ChromKMC
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
KEYOUT=$KMCDIR/$GENBASE.$(printf "chr%02d" "$N").k${K}.kmc_suf
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Pull out each chromosome into a file and generate kmer profiles with KMC..."

  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    grep -A 1 ">$CHR" $SEQIN | tee $KMCDIR/$GENBASE.$CHR.fasta | grep ">"
    kmc -k$K -ci1 -fm $KMCDIR/$GENBASE.$CHR.fasta $KMCDIR/$GENBASE.$CHR.k$K $TEMPDIR
  done && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi

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
## Approach: Calculate kmers shared between chromsomes, or Compleasm gene predictions.
## Rationale: Unless there is a lot of recombination, homeologues will have maximal shared Duplicated BUSCO genes or Kmers.
## Method: KMC intersect of individual chromosome KMC dumps, or Compleasm run.
##    Will bypass if $BEST is pointing to a manually constructed file.
## Key inputs: KMC dumps (for KMC) or assembly fasta (for Compleasm)
## Key outputs: Pairwise KMC intersect CSV, or compleasm full table TSV.
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=ChromIntersect
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
CSV=$GENBASE.k$K.csv
TSV=$GENBASE.${LINEAGE}_odb12.tsv
#i# Set $KEYOUT depending on $BEST method.
KEYOUT=$CSV
if [ "$BEST" == "compleasm" ]; then
  KEYOUT=$TSV
fi
if [ -f "$BEST" ]; then
  KEYOUT=$BEST
fi
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

  # Use compleasm to make the best file
  if [ "$BEST" == "compleasm" ]; then
    echo "[$(date)] Running Compleasm ..." | tee -a $LOG
    COMPRES=$GENBASE.compleasm/${LINEAGE}_odb12/full_table.tsv
    if [ ! -f "$COMPRES" ]; then
      compleasm run -a $SEQIN -o $GENBASE.compleasm -t $THREADS -l $LINEAGE
    fi
    $SHAMDIR/compleasmbest.sh $GENBASE $COMPRES $N $BESTN
    cp -v $GENBASE.best.txt $GENBASE.k$K.best.txt
    cp -v $COMPRES $TSV && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG
  fi

  # Else, will need kmer intersect:
  if [ "$BEST" != "compleasm" ] && [ ! -f "$BEST" ]; then
  	echo "[$(date)] Generating kmer intercepts ..." | tee -a $LOG
  	echo chri,chrj,ktype,knum | tee tmp.$CSV
    for i in $(seq 1 $N); do
      CHR=$(printf "chr%02d" "$i")
      echo chr$i "->" $CHR
      for j in $(seq 1 $N); do
        CHRJ=$(printf "chr%02d" "$j")
        kmc_tools simple $KMCDIR/$GENBASE.$CHR.k$K $KMCDIR/$GENBASE.$CHRJ.k$K intersect $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K
        kmc_dump $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K.txt
        KNUM=$(wc -l $KMCDIR/$GENBASE.${CHR}-${CHRJ}.k$K.txt | awk '{print $1;}')
        echo "$CHR,$CHRJ,kmer,$KNUM" | tee -a tmp.$CSV
      done
    done && mv tmp.$CSV $CSV && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG
  fi

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi

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
### ~~~~ STEP 3: Pair up chromosomes based on chromosome intersects. ~~~~ ###
## Goal: Pair up chromosomes based on all shared kmers or BUSCO genes from previous step.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=PairChrom
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
KEYOUT=$KMCDIR/$GENBASE.$(printf "chr%02d" "$N").alt.fasta
KEYOUTS="$KEYOUT $GENBASE.k$K.best.txt"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Establishing homeologous chromosomes ..." | tee -a $LOG
  if [ -f $BEST ]; then
    echo "[$(date)] Using $BEST as paired chromosomes ..." | tee -a $LOG
    cp -v $BEST $GENBASE.k$K.best.txt
  fi
  if [ ! -f $GENBASE.k$K.best.txt ]; then
    for i in $(seq 1 $N); do
      CHR=$(printf "chr%02d" "$i")
      grep "^$CHR," $CSV | sed 's/,/ /g' | awk '$1 != $2 {print $4, $1, $2;}' | sort -n -r | head -n $BESTN | tee -a $GENBASE.k$K.best.txt
    done
  fi
  awk '{print $2, $3;}' $GENBASE.k$K.best.txt | tee $GENBASE.k$K.best.tmp
  awk '{print $3, $2;}' $GENBASE.k$K.best.txt | tee -a $GENBASE.k$K.best.tmp
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    if [ -f "$KMCDIR/$GENBASE.$CHR.alt.fasta" ]; then
      rm $KMCDIR/$GENBASE.$CHR.alt.fasta
    fi
    for CHRJ in $(grep "^$CHR" $GENBASE.k$K.best.tmp | sort | uniq | awk '{print $2;}'); do
      echo "$CHR alt -> $CHRJ"
      cat $KMCDIR/$GENBASE.$CHRJ.fasta >> $KMCDIR/$GENBASE.$CHR.alt.fasta
    done
  done
  wc $KMCDIR/$GENBASE.*.alt.fasta
  rm -v $GENBASE.k$K.best.tmp && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi
  
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
STEP=AlloKmers
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
KEYOUT=$KMCDIR/$GENBASE.$(printf "chr%02d" "$N").allok${K}.kmc_suf
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Generate allokmers by subtraction of chromosome pairs ..." | tee -a $LOG
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    kmc -k$K -ci1 -fm $KMCDIR/$GENBASE.$CHR.alt.fasta $KMCDIR/$GENBASE.$CHR.alt.k$K $TEMPDIR
    kmc_tools simple $KMCDIR/$GENBASE.$CHR.k$K $KMCDIR/$GENBASE.$CHR.alt.k$K kmers_subtract $KMCDIR/$GENBASE.${CHR}.allok$K
  done && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi
  
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
STEP=AlloKmerIntersect
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
CSV=$GENBASE.allok$K.csv
KEYOUT=$CSV
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Generate pairwise intersects of all allokmers. ..." | tee -a $LOG
  echo chri,chrj,ktype,knum | tee tmp.$CSV
  for i in $(seq 1 $N); do
    CHR=$(printf "chr%02d" "$i")
    echo chr$i "->" $CHR
    for j in $(seq 1 $N); do
      CHRJ=$(printf "chr%02d" "$j")
      KBASE=$KMCDIR/$GENBASE.${CHR}-${CHRJ}.allok$K
      if [ ! -f "$KBASE.txt" ] || [ "$RUNMODE" == "force" ]; then
        kmc_tools simple $KMCDIR/$GENBASE.$CHR.allok$K $KMCDIR/$GENBASE.$CHRJ.allok$K intersect $KBASE
        kmc_dump $KBASE $KBASE.txt
      fi
      KNUM=$(wc -l $KBASE.txt | awk '{print $1;}')
      echo "$CHR,$CHRJ,allo,$KNUM" | tee -a tmp.$CSV
    done
  done && mv tmp.$CSV $CSV && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi
  
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
STEP=Rscript
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
KEYOUT=$GENBASE.shamrock.pdf
KEYOUTS="$KEYOUT $GENBASE.subgenome.1.txt $GENBASE.subgenome.2.txt"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Running Rscript ..." | tee -a $LOG
	Rscript $SHAMDIR/shamrock.R basefile=$GENBASE k=$K partition=$((BESTN + 1)) | tee -a $LOG \
	  && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi
  
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
## Goal: Compile the two subgenome fasta files from the Shamrock clustering.
## Approach: 
## Rationale: 
## Method: 
## Key inputs: 
## Key outputs: 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=SubGenomes
STEPN=$((STEPN + 1))
DONE=$GENBASE.$(printf "%02d" "$STEPN").$STEP.done
START="[$(date)] Step $STEPN $STEP started"
KEYOUT=$GENBASE.subgenome.$((BESTN + 1)).fasta
KEYOUTS="$KEYOUT"
if [ ! -f "$KEYOUT" ] || [ "$RUNMODE" == "force" ] || [ "$RUNMODE" == "$STEP" ]; then

	echo "[$(date)] Generating subgenome fasta files..." | tee -a $LOG
	for P in $(seq 1 $((BESTN + 1))); do
	  PFILE=$GENBASE.subgenome.$P.fasta
	  if [ -f "$PFILE" ]; then
	    rm $PFILE
	  fi
	  for CHR in $(cat $GENBASE.subgenome.$P.txt); do
	    cat $KMCDIR/$GENBASE.$CHR.fasta >> $PFILE
	  done
	  wc $PFILE
	done && echo -e "$START\n[$(date)] Step $STEPN $STEP complete" | tee $DONE | tee -a $LOG

  if [ "$RUNMODE" == "next" ]; then
  	echo "[$(date)] RUNMODE=$RUNMODE -> Exiting run." | tee -a $LOG
  	exit 0;
  fi
  
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
