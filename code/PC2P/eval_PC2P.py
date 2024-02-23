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
    complexes = dict()
    complexfile = "Yeast/complexes_yeast.txt"
    
    ### Read complexes
    with open(complexfile) as f:
        for line in f:
            temp = (line.split())[:2]
            pid, cid = temp[0], temp[1]
            if cid not in complexes.keys():
                complexes[cid] = [pid]
            else:
                complexes[cid].append(pid)
    
    print(complexes)
    
    complexfile_new = "Yeast/CYC2008_scl.txt"
    with open(complexfile_new, "w") as f:
        for cid in complexes:
            f.write("{0} {1}\n".format(cid, ' '.join(complexes[cid])))
    
