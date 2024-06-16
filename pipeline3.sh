#!/bin/bash

set -e

help() {
echo "usage: ./pipeline.sh [-p [ppinfile]] [-r [reffile]] [-o [outputdir]] [-f [filter]] [-h]
    
Runs P5COMP on the given PPIN file (ppinfile) and evaluates against the given gold standard (reffile).
Final predicted clusters will be written in outputdir.
Important: Protein names (PID) should be in gene name (ordered locus) or KEGG format (ex. YLR075W) to match gold standards.

options:
    -p [ppinfile]       path to PPIN file (.txt) where each row is (u v w) (required)
    -r [reffile]        path to gold standard or reference complexes file (.txt) (required)
    -o [outputdir]      path to output directory (required)
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

# Parse user-defined parameters
ppinfile=
reffile=
outputdir=
filtering=
attribs=
resfile=
while getopts ":hp:r:o:f:a:R:" opt; do
    case ${opt} in
        h)
            help
            exit 0
            ;;
        f)
            filtering=$OPTARG
            ;;
        R)
            resfile=$OPTARG
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

echo "RUNNING INDEPENDENT CLUSTERING ALGO: ${attribs}"
echo "Denoising..."
python code/filtering.py $ppinfile $reffile $filteredfile --filtering $filtering

method=${attribs%%-*}
case "$method" in
   "PC2P")
        echo "IND 1: Running PC2P..."
        predictsfile="${outputdir}/PC2P_predicted.txt"
        postprocessed="${outputdir}/PC2P_postprocessed.txt"
        python code/PC2P/PC2P.py $filteredfile $predictsfile -p mp
        python code/PC2P/PC2P_scoring.py $filteredfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile $resfile auc_pts.csv --attribs $attribs
   ;;
   "CUBCO+")
        echo "IND 2: Running CUBCO+..."
        predictsfile="${outputdir}/CUBCO+_predicted.txt"
        postprocessed="${outputdir}/CUBCO+_postprocessed.txt"
        python code/CUBCO+/CUBCO.py $filteredfile $outputdir $predictsfile
        python code/CUBCO+/cubco_scoring.py $filteredfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile $resfile auc_pts.csv --attribs $attribs
   ;;
   "ClusterOne")
        echo "IND 3: Running ClusterOne..."
        predictsfile="${outputdir}/ClusterOne_predicted.txt"
        postprocessed="${outputdir}/ClusterOne_postprocessed.txt"
        jarPath="code/ClusterOne/cluster_one-1.0.jar"
        java -jar $jarPath $filteredfile > $predictsfile

        ## Score ClusterOne clusters
        python code/ClusterOne/cluster_one_scoring.py $ppinfile $predictsfile $postprocessed
        python code/eval2.py $postprocessed $reffile $resfile auc_pts.csv --attribs $attribs
   ;;
esac
