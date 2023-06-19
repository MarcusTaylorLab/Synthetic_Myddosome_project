# Installing and setting up ddG pipeline on Raven
## Installation of ROSETTA3

- Get an academic license for Rosetta: http://c4c.uwc4c.com/express_license_technologies/rosetta

- Downloading Rosetta on raven: 
```
curl -u user_name:password https://www.rosettacommons.org/downloads/academic/2023/wk6/rosetta.source.release-340.tar.bz2 -o rosetta.source.release-340.tar.bz2
```

- Downloading Rosetta on local computer download the source code here: https://www.rosettacommons.org/software/academic

- Uncompress downloaded code, you can do this directly in your home directory:
```
tar xjf rosetta.source.release-340.tar.bz2
```
- Next, navigate to the source directory: ```cd path/to/rosetta/main/source```
- Create a virtual python environment using conda:
```
conda create --name rosetta python=3.10
```
- activate environment with conda ``` activate rosetta ``` and deactivate with ```conda deactivate```
- install the compiler SCons in environment using
```
pip install SCons
```
- The compilation command is
```
./scons.py -j 36 mode=release bin extras=mpi
```
- If you want to install it on your local machine (works only on UNIX operating system) use ```./scons.py -j<number_of_processors_to_use> mode=release bin```. Replace with a number one processor fewer than your computer has. Expect compilation to take a while (hours on one processor).

## Setting up ddG pipeline

- Clone this directory onto raven
```
git clone https://github.com/m4uriz/Synthetic_Myddosome_project
```
- Navigate to the folder rosetta_pipeline ```cd rosetta_pipeline```
- prepare the environment by creating necessary folder structure:
```
mkdir output_files slurm_logs from_af2
```
- open R by typing ```R```
- install the package bio3d with the command ```install.packages("bio3d")```

## Submitting relaxation of a list of structures followed by interface scoring

- create a list containing in each line the relative paths of the structures you want to relax and score. It is important to provide the relative paths starting at your home directory (rosetta_pipeline), as ROSETTA3 has had problems with absolute file paths. A template for this can be found in ```flags/test.txt```
```
from_rcsb/1qys.pdb
from_rcsb/3R2X.pdb
```
- mdodify ```user_parameters``` so it contains the filepath to your ROSETTA3 installation (ROSETTA3) and to your rosetta_pipeline folder (MYROSETTA)
- provide the file path for the list of structures to be analysed to the INPUT variable
- submit the submit_relax_ddg.sh script
```
./submit_relax_ddg.sh
```

## Downloading pdb structures directly from rcsb

from: https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script

- navigate to the folder from_rcsb ```cd from_rcsb```
- the file ```list_file.txt``` is a plain text file that contains a comma separated list of PDB ids. You can modify this list and insert the pdb ids you want to download
- Download the structures in ```list_file.txt``` in the .pdb format:
```
./batch_download.sh -f list_file.txt -p
```
- Obtain full help on the batch download shell script at the command line with: ```./batch_download.sh -h```
