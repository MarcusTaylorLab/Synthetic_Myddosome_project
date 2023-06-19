#!/bin/bash -l
#
# Name of SLURM-job
#SBATCH -J RELAX_ARRAY
#
# Standard output and error:
#SBATCH -o ./job.out.%j
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=120000
#
#SBATCH --mail-type=NONE
#SBATCH --mail-user=lichtenstein@mpiib-berlin.mpg.de
#SBATCH --time=23:00:00

# Define the directory variables
source ../user_parameters.sh

# Define the variable S which is equal to the nth line in pdb_list - pdb_list contains the structures and their paths to be analysed.
S=$( head -${SLURM_ARRAY_TASK_ID} $INPUT_LIST | tail -1 )

echo "Array-ID is:" ${SLURM_ARRAY_TASK_ID}
echo "path to structure is:" $S

# extract structure name from S
filename=$S
title=$(basename "$filename" .pdb)

echo "name of structure is:" $S

# create a folder for each structure to be analysed
mkdir $MYROSETTA/output_files/$title

# create a new variable DIR for this new folder
DIR=$MYROSETTA/output_files/$title

cd $MYROSETTA

# determine the interfaces of chain A in r-script
output=$(R --vanilla --slave -f find_interface.R --args HOME_DIRECTORY=$MYROSETTA STRUCTURE_PATH=$S)
chains=$(echo $output | cut -d '"' -f 2)

echo "chains interfacing with A are:" $chains

# compute 5 relaxed structures for further ddG calculations
$ROSETTA3/bin/relax.mpi.linuxgccrelease -in:file:s $S -relax:constrain_relax_to_start_coords -out:suffix _startconstraints -out:path:all $DIR -nstruct 5

echo "finished relaxation"

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
