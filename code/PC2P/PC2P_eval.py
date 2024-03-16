# Program that evaluates predicted complexes by getting the precision and recall
# python code/PC2P/eval_PC2P.py code/PC2P/results/KroganExt/G_PredictedComplexes_iter0.txt
import argparse
import os
import sys
from helper import printc
from typing import Tuple
import PredictedClusters_Analysis as pc

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
    # Set parameters
    match_thresh = 0.5
    
    # refs: reference complexes in the gold standard
    refs = []
    with open("code\PC2P\Yeast\CYC_complexes.txt") as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the complex id
            proteins = line.split()
            ref = Cluster(proteins, score=None, id=lineno)
            ref.matched = False
            refs.append(ref)
     
    # clusters: predicted clusters
    # NOTE: What if we reformat each line in the predicts file as "score p1 p2 ..."?
    # Clusters are defined here as objects, with predicts as a set of clusters
    # Alternatively, predicts: { cid: { proteins: set(p1, p2, ...), score: float } }
    clusters = []
    with open("code\PC2P\ResultsNew\Collins_CYC_Predicted_iter0.txt") as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters.append(cluster)
    clusters.sort(key = lambda x: x.score, reverse=True)
    
    ### Compute how many correct ref matches is made by a cluster
    for cluster in clusters:
        for ref in refs:
            # Predicted cluster and reference complex is a match if Jaccard(P,C) > 0.5
            if Jaccard(cluster.proteins, ref.proteins) > match_thresh:
                cluster.matches += 1
    
    ### Calculate Precision and Recall for Score Threshold score_thresh
    # precision = TP/(TP+FP) = corrects/predicts, recall = TP/(TP+FN) = corrects/len(refs)
    predicts = 0    # Complexes predicted so far
    corrects = 0
    score_thresh = -1
    rec_thresh = 0.01   # acceptable recall
    auc = 0
    auc_prevrecall = 0

    printc("Threshold\tPreds\tPrecision\tRecall")
    for cluster in clusters:
        if cluster.score != score_thresh:
            recall = corrects/len(refs)   # Without train-test split, all our refs are technically test complexes
            if (score_thresh != -1 and recall > rec_thresh):
                precision = corrects/predicts
                print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, corrects, precision, recall))
                score_thresh = cluster.score
                
                # Prepare for next calculations
                while (rec_thresh < recall): rec_thresh += 0.01
                auc += (recall - auc_prevrecall) * (corrects/predicts)
                auc_prevrecall = recall
            score_thresh = cluster.score
            
        if (cluster.matches > 0): corrects += 1
        
        predicts += 1
            
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    printc("Overall precision and recall")
    TP = 0
    for P in clusters:
        for C in refs:
            if Jaccard(P.proteins, C.proteins) > 0.5:
                TP += 1
                break; # found a match
    
    precision = TP/len(clusters)    # Positive Predictive Value / Accuracy / Quality
    recall = TP/len(refs)           # True Positive Rate / Quantity
    
    print("Precision: %.6f" % precision)
    print("Recall: %.6f" % recall)
    
    printc("=========OMRANIAN==========")
    print("Precision: %.6f" % pc.precision_Jaccard(refs, clusters))
    print("Recall: %.6f" % pc.recall_Jaccard(refs, clusters))
    print("F-Measure: %.6f" % pc.F_measure_Jaccard(refs, clusters))
        
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity