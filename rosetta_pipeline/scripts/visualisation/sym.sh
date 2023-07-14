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
# Number of nodes and MPI tasks per node:
#SBATCH --ntasks=1
#SBATCH --constraint="gpu"
#
# We will use 4 GPUs:
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=72
#SBATCH --mem=500000
#
#SBATCH --mail-type=NONE
#SBATCH --mail-user=lichtenstein@mpiib-berlin.mpg.de
#SBATCH --time=24:00:00

conda activate rosetta

# Load compiler and MPI modules with explicit version specifications,
# consistently with the versions used to build the executable.
module purge
#module load intel/19.1.3 impi/2019.9
module load intel/21.7.1 impi/2021.7
module load cuda/11.2

# Define the directory variables
ROSETTA3='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/main/source'
MYROSETTA='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline'

# put temporary files into a ramdisk
export TMPDIR=${JOB_SHMTMPDIR}

### ENABLE CUDA UNIFIED MEMORY ###############################################
export TF_XLA_FLAGS=--tf_xla_enable_xla_devices
export TF_FORCE_UNIFIED_MEMORY=1
     # Enable jax allocation tweak to allow for larger models, note that
     # with unified memory the fraction can be larger than 1.0 (=100% of single GPU memory):
     # https://jax.readthedocs.io/en/latest/gpu_memory_allocation.html
     # When using 4 GPUs:
export XLA_PYTHON_CLIENT_MEM_FRACTION=4

     # run threaded tools with the correct number of threads (MPCDF customization)
export NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd /u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline/output_files/sym

$ROSETTA3/src/apps/public/symmetry/make_symmdef_file.pl -m HELIX -a U -b T -t 10 -p 6E9X_startconstraints_0001_renumbered.pdb > 6E9X_startconstraints_0001_renumbered.symm

echo 'finished symmetry'
