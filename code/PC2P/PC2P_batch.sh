#!/bin/bash

codepath="code/PC2P/PC2P.py"
inputdir="code/PC2P/Yeast/FilteredPPINs"
outputdir="code/PC2P/Results/FilteredRay"
pool_thresh=100

for inputfile in "$inputdir"/*
do
    python $codepath $inputfile $outputdir -p --pool_thresh $pool_thresh
done