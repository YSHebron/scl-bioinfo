# Program that evaluates predicted complexes by getting the precision and recall
# python code/PC2P/eval_PC2P.py code/PC2P/results/KroganExt/G_PredictedComplexes_iter0.txt
import argparse
import os
import sys
from helper import printc
from typing import Tuple

# Calculate Jaccard similarity between P and C
# P: Preducted cluster
# C: Reference complex

def jaccard_similarity(P, C):
    intersection_size = len(P.intersection(C))
    union_size = len(P.union(C))
    
    if union_size == 0:
        return 0.0
    
    return intersection_size / union_size

if __name__ == '__main__':
    ### complexes will contain gold standard
    complexes = dict()
    complexfile = "code/PC2P/Yeast/CYC2008_scl.txt"
    
    with open(complexfile) as f:
        for line in f:
            temp = line.split(' ', 1)
            cid, pids = temp[0], temp[1].strip().split()
            complexes[cid] = pids
            
    print(complexes)
    
    # pairplot for precision-recall
    # plan out results format (using set of pairplots colorcoded according to either gold standard/ppi dataset/approach/clustering algorithm used used)
