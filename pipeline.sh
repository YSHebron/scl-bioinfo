#!/bin/bash

# P5COMP: Parameter-free Pipeline for Predicting Problematic Protein Complexes
# Can only currently be run from repository root.
# Disclaimer: Parameter-free in a very specific, clustering sense.

set -e

help() {
echo "usage: ./pipeline.sh [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W).

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v s) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)
    -h                  show this help information"
}

log() {
    echo "$@" >> ./debug.txt
}

validate_file() {
    if [[ ! -f $1 ]]; then 
        echo "$1 is not a valid file."
        exit 1
    fi
}

validate_dir() {
    if [[ ! -d $1 ]]; then 
        echo "$1 is not a valid directory."
        exit 1
    fi
}

if [ $# -eq 0 ]; then
    echo "No options supplied."
    help
    exit 1
fi

# Developer parameters
filteredfile="data/Interm/filtered_ppin.txt"
decompfile="data/Interm/decomp_ppin.txt"
hubfile="data/Interm/hub_proteins.txt"
iAdjustCD_outfile="data/Interm/ppin_adjusted.txt"

# Parse user-defined parameters
ppinfile=
reffile=
outputdir=
while getopts ":hp:r:o:n:" opt; do
    case ${opt} in
        h)
            help
            exit 0
            ;;
        p)
            ppinfile=$OPTARG
            ;;
        r)
            reffile=$OPTARG
            ;;
        o)
            outputdir=$OPTARG
            ;;
        n)
            negfile=$OPTARG
            ;;
        \?)
            echo "Error: Invalid option: -$OPTARG."
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done

if [[ -z $ppinfile || -z $reffile || -z $outputdir ]]; then
    echo "Error: Missing -p, -r, and/or -o arguments. See -h."
    exit 1
fi
validate_file $ppinfile && validate_file $reffile && validate_dir $outputdir

echo "Running P5COMP..."
printf "PPIN:\t%s\n" $(realpath "$ppinfile" -q)
printf "Ref:\t%s\n" $(realpath "$reffile" -q)
printf "Output:\t%s\n" $(realpath "$outputdir" -q)

# Denoising -> data/Interm/filtered_ppin.txt
## Filtering: Negatome and PerProteinPair filtering.
## Note: This pipeline is packaged with Negatome 2.0 datasets.
python code/filtering.py $ppinfile $reffile $filteredfile --negfile $negfile --confidence 0.33

### DECOMP 1: Hub Removal
### -> data/Interm/ppin_adjusted.txt
### -> data/Interm/decomp_ppin.txt, data/Interm/hub_proteins.txt 
python code/iAdjustCD.py $filteredfile $iAdjustCD_outfile
python code/hub_remove.py $filteredfile $hubfile $decompfile

# Parallel Clustering
## 1. PC2P
predictsfile_PC2P="${outputdir}/PC2P_predicted.txt"
postprocessed_PC2P="${outputdir}/PC2P_postprocessed.txt"
python code/PC2P/PC2P.py $decompfile $predictsfile_PC2P -p mp

### DECOMP 2: Hub Return -> data/Results/Dummy/{method}_predicted.txt
python code/hub_return.py $predictsfile_PC2P $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_PC2P

## 2. CUBCO+


## 3. FINCH


## 4. ClusterOne


# Ensemble Clustering

# Evaluation
python code/PC2P/PC2P_eval.py $postprocessed_PC2P $reffile $outputdir