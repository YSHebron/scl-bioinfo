#!/bin/bash

# P5COMP: Parameter-free Pipeline for Predicting Problematic Protein Complexes
# Can only currently be run from repository root.
# Disclaimer: This still contain necessary usage parameters such as paths to files and directories
# Note that usage parameters are not clustering parameters.

# Sample: ./pipeline.sh -p data/Yeast/Collins/collins2007.txt -r data/Yeast/CYC_complexes.txt -o data/Results/P5COMP -n data/Negatome/negatome_2_mix_mapped.txt -f perpair

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
    -a [attribs]        attributes for evaluation file name of format 'algo-goldstd-ppin', ex: P5COMP-CYC-Collins
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
attribs=
while getopts ":hp:r:o:n:f:a:" opt; do
    case ${opt} in
        h)
            help
            exit 0
            ;;
        f)
            filtering=$OPTARG
            ;;
        a)
            attribs=$OPTARG
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

echo "RUNNING INDEPENDENT CLUSTERING ALGOS:"
echo "Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --filtering $filtering

echo "IND 1: Running PC2P..."
predictsfile_PC2P="${outputdir}/PC2P_predicted.txt"
python code/PC2P/PC2P.py $filteredfile $predictsfile_PC2P -p mp

echo "IND 2: Running CUBCO+..."
predictsfile_CUBCO="${outputdir}/CUBCO+_predicted.txt"
python code/CUBCO+/CUBCO.py $filteredfile $outputdir $predictsfile_CUBCO

echo "IND 3: Running ClusterOne..."
predictsfile_ClusterOne="${outputdir}/ClusterOne_predicted.txt"
postprocessed_ClusterOne="${outputdir}/ClusterOne_postprocessed.txt"
jarPath="code/ClusterOne/cluster_one-1.0.jar"
java -jar $jarPath $filteredfile > $predictsfile_ClusterOne

## Score ClusterOne clusters
python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile_ClusterOne $postprocessed_ClusterOne

# Evaluation (currently assumes running from root)
python code/eval2.py $predictsfile_PC2P $reffile results.csv auc_pts.csv --attribs $attribs
python code/eval2.py $predictsfile_CUBCO $reffile results.csv auc_pts.csv --attribs $attribs
python code/eval2.py $postprocessed_ClusterOne $reffile results.csv auc_pts.csv --attribs $attribs