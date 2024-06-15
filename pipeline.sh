#!/bin/bash

# P5COMP: Parameter-free Pipeline for Predicting Problematic Protein Complexes
# Can only currently be run from repository root.
# Disclaimer: This still contain necessary usage parameters such as paths to files and directories
# Note that usage parameters are not clustering parameters.

set -e

help() {
echo "usage: ./pipeline.sh [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-n [negfile]] [-f [filter]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
    -n [negfile]        path to negatome (.txt) where each row is (u v) (optional)
    -f [filter]         filtering type (perpair or perprotein)
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
    echo "Error: No options supplied. See help below."
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
filtering=
while getopts ":hp:r:o:n:f:" opt; do
    case ${opt} in
        h)
            help
            exit 0
            ;;
        f)
            filtering=$OPTARG
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
    echo "Error: Missing -p, -r, and/or -o arguments. See help (-h)."
    exit 1
fi
mkdir -p $outputdir
validate_file $ppinfile && validate_file $reffile && validate_dir $outputdir

echo "Running P5COMP..."
printf "PPIN:\t%s\n" $(realpath "$ppinfile" -q)
printf "Ref:\t%s\n" $(realpath "$reffile" -q)
printf "Output:\t%s\n" $(realpath "$outputdir" -q)

# Denoising -> data/Interm/filtered_ppin.txt
## Filtering: Negatome and either PerProteinPair / PerProtein filtering.
## Note: This pipeline is packaged with Negatome 2.0 datasets.
echo "1. Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --negfile $negfile --filtering $filtering

### DECOMP 1: Hub Removal
### -> data/Interm/ppin_adjusted.txt
### -> data/Interm/decomp_ppin.txt, data/Interm/hub_proteins.txt 
echo "2. Hub Removal..."
python code/iAdjustCD.py $filteredfile $iAdjustCD_outfile
python code/hub_remove.py $filteredfile $hubfile $decompfile
# At this point, hubfile contains the hub proteins while decompfile contains the decomposed PPIN.
# The iAdjustCD_outfile contains the rescored PPIN, for use in hub_return.

echo "Clustering Algorithms..."
# Parallel Clustering (Ensemble Clustering)
## 1. PC2P with hub return -> data/Results/Dummy/PC2P_predicted.txt, data/Results/Dummy/PC2P_postprocessed.txt
echo "Running PC2P..."
predictsfile_PC2P="${outputdir}/PC2P_predicted.txt"
postprocessed_PC2P="${outputdir}/PC2P_postprocessed.txt"
python code/PC2P/PC2P.py $decompfile $predictsfile_PC2P -p mp
python code/hub_return.py $predictsfile_PC2P $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_PC2P

## NOTE: UNCOMMENT ONLY THE FOLLOWING python LINES ONCE CUBCO+ and ClusterOne CAN SUPPORT DECOMP
## For parallelization, might separate clustering algorithms from DECOMP (but seq ok for now)
## TODO: Simplify variables (predictsfile_method to just predictsfile, etc.)

## 2. CUBCO+ with hub return -> data/Results/Dummy/CUBCO+_predicted.txt, data/Results/Dummy/CUBCO+_postprocessed.txt
## Note: omit '+' character from varnames
echo "Running CUBCO+..."
predictsfile_CUBCO="${outputdir}/CUBCO+_predicted.txt"
postprocessed_CUBCO="${outputdir}/CUBCO+_postprocessed.txt"
python code/CUBCO+/CUBCO.py $decompfile $outputdir $predictsfile_CUBCO
python code/hub_return.py $predictsfile_CUBCO $iAdjustCD_outfile $hubfile $filteredfile $postprocessed_CUBCO

## 3. ClusterOne with hub return -> data/Results/Dummy/ClusterOne_predicted.txt, data/Results/Dummy/ClusterOne_postprocessed.txt
# Insert ClusterOne code
echo "Running ClusterOne..."
predictsfile_ClusterOne="${outputdir}/ClusterOne_predicted.txt"
postprocessed_ClusterOne="${outputdir}/ClusterOne_postprocessed.txt"
jarPath="code/ClusterOne/cluster_one-1.0.jar"
java -jar $jarPath $filteredfile > $predictsfile_ClusterOne

## Score clusters
python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile_ClusterOne $postprocessed_ClusterOne

# Ensemble Clustering
echo "Finale: Running Ensemble Clustering..."
P5COMP_clusters="${outputdir}/P5COMP_clusters.txt"
python code/ensemble.py $postprocessed_ClusterOne $postprocessed_CUBCO $postprocessed_PC2P $P5COMP_clusters

# Evaluation (currently assumes running from root)
python code/eval.py $postprocessed_ClusterOne $postprocessed_CUBCO $postprocessed_PC2P $P5COMP_clusters $reffile $outputdir
