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
install the compiler SCons in environment using
```
pip install SCons
```
- The compilation command is
```
./scons.py -j 36 mode=release bin extras=mpi
```
- If you want to install it on your local machine (works only on UNIX operating system) use ```./scons.py -j<number_of_processors_to_use> mode=release bin```. Replace with a number one processor fewer than your computer has. Expect compilation to take a while (hours on one processor).

## Setting up ddG pipeline
