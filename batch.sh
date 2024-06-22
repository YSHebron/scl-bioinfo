#!/bin/bash
echo "Method,GldStd,PPIN,Predicts,Refs,Precision,Recall,F1-score,F2-score,AUC-PR,MMR,Sensitivity,PPP,Accuracy,F-Match,Separation" > results.csv

for i in {1..1}
do
    # metrics.sh -resultsfile -recreate
    ./metrics.sh results.csv false
done