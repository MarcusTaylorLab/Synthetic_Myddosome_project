#!/bin/bash -l
#
# Name of SLURM-job
#SBATCH -J RELAX_SINGLE
#
# Initial working directory:
#SBATCH -D '/u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline/slurm_logs'
#
# Standard output and error:
#SBATCH -o ./job.out.%j
#
# Run a job on 2 cores of a shared node using CPU
#SBATCH --ntasks=1
#
#SBATCH --cpus-per-task=2
#SBATCH --mem=5000
#
#SBATCH --mail-type=NONE
#SBATCH --mail-user=lichtenstein@mpiib-berlin.mpg.de
#SBATCH --time=23:00:00

# Define the directory variables
ROSETTA3='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/main/source'
MYROSETTA='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline'

# Provide the relative path of structure to be analysed
S='from_rcsb/6I3N.pdb'

echo "path to structure is:" $S

# extract structure name from S
filename=$S
title=$(basename "$filename" .pdb)

echo "name of structure is:" $title

# create a folder for each structure to be analysed
mkdir $MYROSETTA/output_files/$title/cpu

# create a new variable DIR for this new folder
DIR=$MYROSETTA/output_files/$title/cpu

cd $MYROSETTA

# determine the interfaces of chain A in r-script
output=$(R --vanilla --slave -f find_interface.R --args HOME_DIRECTORY=$MYROSETTA STRUCTURE_PATH=$S)
chains=$(echo $output | cut -d '"' -f 2)

echo "chains interfacing with A are:" $chains

# compute 5 relaxed structures for further ddG calculations
$ROSETTA3/bin/relax.mpi.linuxgccrelease -in:file:s $S -relax:constrain_relax_to_start_coords -out:suffix _startconstraints -out:path:pdb $DIR -nstruct 5

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
