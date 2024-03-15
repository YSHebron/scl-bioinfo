# Program that evaluates predicted complexes by getting the precision and recall
# python code/PC2P/eval_PC2P.py code/PC2P/results/KroganExt/G_PredictedComplexes_iter0.txt
import argparse
import os
import sys
from helper import printc
from typing import Tuple

# To emulate Yong and Wong for predicts, we also add number of correct matches
class Cluster:
    def __init__(self, proteins = set(), score = 0, id = None):
        self.proteins = set(proteins)
        self.score = score
        self.id = id
    
    def __str__(self):
        return "(%d_%.6f): %s" % (len(self.proteins), self.score, " ".join(self.proteins))

# Calculate Jaccard similarity between P and C
# P: Predicted cluster
# C: Reference complex
def Jaccard(P, C):
    return len(P.intersection(C)) / len(P.union(C))

if __name__ == '__main__':
    # gldstd: reference complexes in the gold standard
    gldstd = []
    with open("code\PC2P\Yeast\CYC2008_complexes.txt") as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the complex id
            proteins = line.split()
            gldstd.append(Cluster(proteins, 1, lineno))
     
    # predicts: predicted clusters
    # NOTE: What if we reformat each line in the predicts file as "score p1 p2 ..."?
    # Clusters are defined here as objects, with predicts as a set of clusters
    # Alternatively, predicts: { cid: { proteins: set(p1, p2, ...), score: float } }
    predicts = []
    with open("code\PC2P\Results\Collins_CYC_Predicted_iter0.txt") as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[1].strip("):"))
            predicts.append(Cluster(proteins, score, lineno))
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    TP = 0
    for i in range(0, len(predicts)):
        for j in range(0, len(gldstd)):
            P, C = predicts[i], gldstd[j]
            # Predicted cluster and reference complex is a match if Jaccard(P,C) > 0.5
            if Jaccard(P.proteins, C.proteins) > 0.5:
                TP += 1
                break; # found a match
    
    precision = TP/len(predicts)    # Positive Predictive Value / Accuracy / Quality
    recall = TP/len(gldstd)         # True Positive Rate / Quantity
    
    print("Precision: %.6f" % precision)
    print("Recall: %.6f" % recall)