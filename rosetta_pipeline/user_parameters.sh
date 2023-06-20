#!/bin/bash -l

# Define environment variables
ROSETTA3='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/main/source'
MYROSETTA='/u/mlicht/rosetta_2023.6/rosetta.source.release-340/Synthetic_Myddosome_project/rosetta_pipeline'

# Define the path to the structures to be processed
INPUT_LIST=$MYROSETTA/flags/test.txt

# Compute number of structures in list
NUM=$(wc -l < $INPUT_LIST)

NUM=$((NUM+1))

echo "Number of structures:" $NUM
