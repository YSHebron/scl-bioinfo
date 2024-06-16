#!/bin/bash

echo "Method, GldStd, PPIN, Predicts, Refs, Precision, Recall, F1-score, F2-score, AUC-PR, MMR, Sensitivity, PPP, Accuracy, F-Match, Separation" > results.csv

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

p=
r=
o="data/Results/P5COMP"
for gldstd in "${gldstds[@]}"; do
    if [ "$gldstd" = "CYC" ]; then
        r=eval/CYC2008.txt
    elif [ "$gldstd" = "SGD" ]; then
        r=eval/SGD.txt
    fi
    for ppin in "${ppins[@]}"; do
        if [ "$ppin" = "Collins" ]; then
            p=eval/Collins.txt
        elif [ "$ppin" = "Gavin" ]; then
            p=eval/Gavin.txt
        elif [ "$ppin" = "KroganCore" ]; then
            p=eval/KroganCore.txt
        elif [ "$ppin" = "KroganExt" ]; then
            p=eval/KroganExt.txt
        elif [ "$ppin" = "BIM" ]; then
            p=eval/BIM.txt
        fi
        # P5COMP
        ./pipeline2.sh -p $p -r $r -o $o \
            -n data/Negatome/negatome_2_mix_mapped.txt -f perpair -a "P5COMP-${gldstd}-${ppin}"
        # PC2P, CUBCO+, and ClusterOne
        for method in "${methods[@]}"; do
            ./pipeline3.sh -p $p -r $r -o $o -f perprotein -a "${method}-${gldstd}-${ppin}"
        done
    done
done