#!/bin/bash

# P5COMP: Parameter-free Pipeline for Predicting Problematic Protein Complexes

set -e

help() {
echo "usage: ./pipeline.sh [-i [ppinfile]] [-r [reffile]] [-o [outputdir]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.

options:
    -i [ppinfile]       path to PPIN file (.txt, .csv, .tsv) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -h                  show this help information"
}

log() {
    echo "$@" >> ./debug.txt
}

if [ $# -eq 0 ]
then
echo "No options supplied."
help
exit 1
fi

# Parse command line
ppinfile=
reffile=
outputdir=
while getopts ":hi:r:o:" opt; do
    case ${opt} in
        h)
            help
            exit 0
            ;;
        i)
            ppinfile=$OPTARG
            if [[ ! -f $ppinfile ]]; then 
                echo "$ppinfile is not a valid file."
                exit 1
            fi
            ;;
        r)
            reffile=$OPTARG
            if [[ ! -f $reffile ]]; then 
                echo "$reffile is not a valid file."
                exit 1
            fi
            ;;
        o)
            outputdir=$OPTARG
            if [[ ! -d $outputdir ]]; then 
                echo "$outputdir is not a valid directory."
                exit 1
            fi
            ;;
        \?)
            echo "Invalid option: -$OPTARG."
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires a path argument."
            exit 1
            ;;
    esac
done

echo "Running P5COMP..."
printf "CWD:\t%s\n" "$(pwd)"
printf "PPIN:\t%s\n" "$(realpath "$ppinfile")"
printf "Ref:\t%s\n" "$(realpath "$reffile")"
printf "Output:\t%s\n" "$(realpath "$outputdir")"

# Denoising
## Filtering

# Parallel Clustering

## 1. PC2P


## 2. CUBCO+


## 3. FINCH


## 4. ClusterOne


# Ensemble Clustering

# Evaluation