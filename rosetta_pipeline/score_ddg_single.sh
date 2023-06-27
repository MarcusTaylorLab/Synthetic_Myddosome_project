#!/bin/bash -l
#
# Name of SLURM-job
#SBATCH -J RELAX_ARRAY
#
# Standard output and error:
#SBATCH -o job_%A_%a.out        # Standard output, %A = job ID, %a = job array index
#
# Number of nodes and MPI tasks per node:
#SBATCH --ntasks=1
#SBATCH --constraint="gpu"
#
# --- default case: use a single GPU on a shared node ---
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000
#
# Define the SLURM-array to analyse multiple pdb structures in parallel
#SBATCH --array=1-1
#
#SBATCH --mail-type=NONE
#SBATCH --mail-user=lichtenstein@mpiib-berlin.mpg.de
#SBATCH --time=23:00:00

# Define the directory variables
ROSETTA3='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/main/source'
MYROSETTA='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline'

# Define the path to the structures to be processed
#INPUT_LIST=$MYROSETTA/flags/score_list.txt

# Compute number of structures in list
#NUM=$(wc -l < $INPUT_LIST)

#NUM=$((NUM+1))

#echo "Number of structures:" $NUM

# Define the variable S which is equal to the nth line in pdb_list - pdb_list contains the structures and their paths to be analysed.
#S=$( head -${SLURM_ARRAY_TASK_ID} $INPUT_LIST | tail -1 )

#echo "Array-ID is:" ${SLURM_ARRAY_TASK_ID}

S="from_rcsb/HA13_trimmed.pdb"

echo "path to structure is:" $S

# extract structure name from S
filename=$S
title=$(basename "$filename" .pdb)

echo "name of structure is:" $S

# create a new variable DIR for this new folder
DIR=$MYROSETTA/output_files/$title

cd $MYROSETTA

# determine the interfaces of chain A in r-script
#output=$(R --vanilla --slave -f find_interface.R --args HOME_DIRECTORY=$MYROSETTA STRUCTURE_PATH=$S)
#chains=$(echo $output | cut -d '"' -f 2)

chains="Q V W"

echo "chains interfacing with A are:" $chains

cd $DIR

# compute the total_score all-atom-rmsd and ca-rmsd. Output is named default.sc
$ROSETTA3/bin/score.default.linuxgccrelease -in:file:s $DIR/*.pdb -in:file:native $MYROSETTA/$S -ignore_waters
awk '{print $NF, $2, $(NF-11), $(NF-1)}' default.sc | tail -n +1 > total_scores_$title.sc

# compute the different interface energies using InterfaceAnalyzer
for i in $chains
do
$ROSETTA3/bin/InterfaceAnalyzer.mpi.linuxgccrelease -in:file:s $DIR/*.pdb -interface A_$i -pack_input true -pack_separated true -out:file:score_only interface_A.sc -compute_packstat
mv interface_A.sc interface_A$i.sc
echo dG_A$i > dG_A$i.sc
awk '{print $6}' interface_A$i.sc | tail -n +3 >> dG_A$i.sc
done

paste total_scores_$title.sc dG* > scores_$title.sc
rm dG* total_scores_$title.sc
