#!/bin/bash -l

# Define environment variables
ROSETTA3='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/main/source'
MYROSETTA='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/data'

# Define the path to the structures to be processed
INPUT_LIST=$MYROSETTA/flags/pdb_list

# Compute number of structures in list
NUM=$(wc -l INPUT_LIST)