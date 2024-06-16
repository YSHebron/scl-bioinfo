#!/bin/bash

declare -a methods=(
    [1]=PC2P
    [2]=CUBCO+
    [3]=ClusterOne
)

declare -a gldstds=(
    [0]=CYC
    [1]=SGD
)

declare -a ppins=(
    [0]=Collins
    [1]=Gavin
    [2]=KroganCore
    [3]=KroganExt
    [4]=BIM
)

# Results File
resultsfile=$1
echo "Method, GldStd, PPIN, Predicts, Refs, Precision, Recall, F1-score, F2-score, AUC-PR, MMR, Sensitivity, PPP, Accuracy, F-Match, Separation" > $resultsfile

p=
r=
o="data/Results/P5COMP"
for gldstd in "${gldstds[@]}"; do
    case "$gldstd" in
        "CYC") r=eval/CYC2008.txt
        ;;
        "SGD") r=eval/SGD.txt
        ;;
    esac
    for ppin in "${ppins[@]}"; do
        case "$ppin" in
            "Collins") p=eval/Collins.txt
            ;;
            "Gavin") p=eval/Gavin.txt
            ;;
            "KroganCore") p=eval/KroganCore.txt
            ;;
            "KroganExt") p=eval/KroganExt.txt
            ;;
            "BIM") p=eval/BIM.txt
            ;;
        esac
        # P5COMP
        ./pipeline2.sh -p $p -r $r -o $o \
            -n data/Negatome/negatome_2_mix_mapped.txt -f perpair -a "P5COMP-${gldstd}-${ppin}" -R $resultsfile
        # PC2P, CUBCO+, and ClusterOne
        for method in "${methods[@]}"; do
            ./pipeline3.sh -p $p -r $r -o $o -f perpair -a "${method}-${gldstd}-${ppin}" -R $resultsfile
        done
    done
done