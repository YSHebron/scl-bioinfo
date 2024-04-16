#!/bin/bash
# This program can both perform batch clustering and batch evaluation of the results of the clustering.

codepath="code/PC2P/PC2P.py"
inputdir="data/Yeast/FilteredPPINs"
outputdir="data/Results/PC2P/FilteredRay"
pool_thresh=100

for inputfile in "$inputdir"/*
do
    python $codepath $inputfile $outputdir -p --pool_thresh $pool_thresh
    echo ================
done