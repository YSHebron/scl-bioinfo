# Program that evaluates predicted complexes by getting the precision and recall
# Do eval.py --help for help on running this program.
# Sample run: python code/eval.py data/Results/Dummy/Trial/ClusterOne_postprocessed.txt data/Results/Dummy/Trial/CUBCO+_postprocessed.txt data/Results/Dummy/Trial/PC2P_postprocessed.txt data/Results/Dummy/Trial/P5COMP_clusters.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/Trial
# python code/eval.py data/Results/Dummy/Raw/Collins_CYC/ClusterOne_postprocessed.txt data/Results/Dummy/Raw/Collins_CYC/CUBCO+_postprocessed.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/Raw
# python code/eval.py data/Results/Dummy/P5COMP/ClusterOne_postprocessed.txt data/Results/Dummy/P5COMP/CUBCO+_postprocessed.txt data/Results/Dummy/P5COMP/PC2P_postprocessed.txt data/Results/Dummy/P5COMP/P5COMP_clusters.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/P5COMP

import argparse
import PredictedClusters_Analysis as pc
import pandas as pd
from utils import printc
from classes import Cluster
from pathlib import Path

parser = argparse.ArgumentParser(description='Evaluate results of clustering algorithm. Must follow line format `(size_score): p1 p2 ...`')
parser.add_argument('predictsfile', type=Path, help='path to predicted clusters file')
parser.add_argument('complexfile', type=Path, help='path to gold standard complex file')
parser.add_argument('outfile_M', type=Path, help='path to output file for results, in .csv format')
parser.add_argument('outfile_A', type=Path, help='path to output file for AUC-PR, in .csv format')
parser.add_argument('--attribs', type=str, help='evaluation attribute of format "algo-goldstd-ppin", ex: "P5COMP-CYC-Collins"')
args = parser.parse_args()

# Calculate Jaccard similarity between P and C
# P: Predicted cluster
# C: Reference complex
def Jaccard(P, C):
    return len(P.intersection(C)) / len(P.union(C))

# Calculate C matches P
def CmatchesP(C, P, match_thresh):
    lenC = len(C.proteins)
    lenP = len(P.proteins)
    # require exact match for small comps
    if (lenC <= 3 and lenP <= 3):
        if (lenC != lenP):
            return 0
        num_intersect = len(C.proteins.intersection(P.proteins))
        if num_intersect == lenC:
            return 1
        else:
            return 0

    # if one is small comp, return 0
    if lenC <= 3 or lenP <= 3:
        return 0
    
    ### Jaccard
    if (lenC > lenP and lenP/lenC < match_thresh) or (lenP > lenC and lenC/lenP < match_thresh):
        return 0
    
    return Jaccard(C.proteins,P.proteins)

def calc_prec_rec_comp_pred(clusters, refs, matched_complexes_ref, outputfile: Path, matchscore_thr=0.5, quiet=True) -> float:
    results = []

    correct_clusters = {}
    correct_smallclusters = {}

    clusters_copy = clusters.copy()

    for cluster in clusters:
        cluster_matches = [cluster.score, 0, cluster, {}]
        for ref in refs:
            matchscore = CmatchesP(cluster, ref, matchscore_thr)
            if matchscore >= matchscore_thr:
                # print(f"{str(ref.proteins)} and {str(cluster.proteins)} are correct")
                correct_clusters[cluster] = 1
                if len(cluster.proteins) <= 3:
                    correct_smallclusters[cluster] = 1

                cluster_matches[1] += 1
                cluster_matches[3][ref] = 1
                matched_complexes_ref[ref] = 1
            
        # Add to the results if the cluster has a match
        if cluster_matches[1] > 0:
            results.append(cluster_matches)

    for cluster in correct_clusters:
        clusters_copy.remove(cluster)

    # append to the results array the incorrect clusers
    for cluster in clusters_copy:
        tmparray = [cluster.score, 0, cluster, {}]
        results.append(tmparray)

    ### Calculate Precision and Recall for Score Threshold score_thresh output to outputdir
    # precision = TP/(TP+FP) = corrects/predicts, recall = TP/(TP+FN) = matched/len(refs)
    predicts = 0    # Complexes predicted so far
    corrects = 0
    matched = 0
    score_thresh = -1
    rec_thresh = 0.01   # acceptable recall
    auc = 0
    auc_prevrecall = 0
    matched_complexes = {}

    auc_points = pd.DataFrame(columns=["Threshold", "Preds", "Precision", "Recall"])
    if not quiet:
        print("Threshold\tPreds\tPrecision\tRecall")
    for cluster in sorted(results, key=lambda x: x[0], reverse=True):
        if cluster[0] != score_thresh:
            recall = matched/len(refs)   # Without train-test split, all our refs are technically test complexes
            if (score_thresh != -1 and recall > rec_thresh):
                precision = corrects/predicts
                if not quiet:
                    print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, predicts, precision, recall))
                auc_points.loc[len(auc_points)] = [score_thresh, predicts, precision, recall]
                
                # Prepare for next calculations
                while (rec_thresh < recall):
                    rec_thresh += 0.01
                auc += (recall - auc_prevrecall) * (corrects/predicts)
                auc_prevrecall = recall
            score_thresh = cluster[0]
            
        if (cluster[1] > 0):
            corrects += 1
            # Mark the predicted complexes as matched in a dictionary
            for comp in cluster[3]:
                matched_complexes[comp] = 1
    
            # Update the count of matched complexes
            matched = len(matched_complexes)
        
        predicts += 1

    if predicts > 0:
        recall = matched / len(refs)
        if not quiet:
            print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, predicts, corrects/predicts, recall))
        auc_points.loc[len(auc_points)] = [score_thresh, predicts, precision, recall]
    else:
        recall = matched / len(refs)
        precision = 0
        if not quiet:
            print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, predicts, corrects/predicts, recall))
        auc_points.loc[len(auc_points)] = [score_thresh, predicts, precision, recall]

    auc += (recall - auc_prevrecall) * (corrects/predicts)
    auc_points.loc[len(auc_points)] = [score_thresh, predicts, precision, recall]
    
    return auc


if __name__ == '__main__':
    predictsfile: Path = args.predictsfile
    complexfile: Path = args.complexfile
    outfile_M: Path = args.outfile_M
    outfile_A: Path = args.outfile_A
    printc("Starting evaluation of %s..." % args.attribs, "red")
    printc("predictsfile:\t%s" % predictsfile)
    printc("complexfile:\t%s" % complexfile)
    printc("outfile_M:\t%s" % outfile_M)    # all metrics, aside from AUC-PR plot points
    printc("outfile_A:\t%s" % outfile_A)    # AUC-PR plot points

    # Set parameters
    match_thresh = 0.75

    # Note: pc (imported) measures accept format of list of clusters (sets)
    # We generate these as follows
    predicts_set, refs_set = [], []
    with open(predictsfile) as f1, open(complexfile) as f2:
        for line in f1:
            cluster = line.strip(' \n').split(':')[1].split()
            predicts_set.append(set(cluster))
        for line in f2:
            complex = line.strip(' \n').split()
            refs_set.append(set(complex))
    
    # Goal: Create a list of metrics that we can append to outfile
    # Format: "Method, GldStd, PPIN, Precision, Recall, F1-score, F2-score,
    #          AUC-PR, Separation, F-Match, Accuracy, Sensitivity, PPP, MMR"
    # if outfile_M still doesn't exist, recreate it
    df = pd.read_csv(outfile_M, index_col=False)
    
    # Setup for AUC-PR (requires scored clusters)
    ## refs: reference complexes in the gold standard
    refs = []
    with complexfile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the complex id
            proteins = line.split()
            ref = Cluster(proteins, score=None, id=lineno)
            ref.matched = False # default assumption
            refs.append(ref)

    ## clusters: predicted clusters by the algorithm
    clusters = []
    with predictsfile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters.append(cluster)
    clusters.sort(key = lambda x: x.score, reverse=True)

    # print(calc_prec_rec_comp_pred(clusters, refs, {}, outfile_A, quiet=True))
    
    # print("Precision:\t%.6f" % pc.precision_Jaccard(refs, clusters))
    # print("Recall:\t\t%.6f" % pc.recall_Jaccard(refs, clusters))
    # print("F-score:\t%.6f" % pc.F_measure_Jaccard(refs, clusters))
    
    # ID VARIABLES
    method, gldstd, ppin = args.attribs.split('-')
    num_predicts = len(clusters)
    num_refs = len(refs)
    
    # METRICS
    precision = pc.precision(refs_set, predicts_set)
    recall = pc.recall(refs_set, predicts_set)
    f1 = pc.F_measure(refs_set, predicts_set, F=1)
    f2 = pc.F_measure(refs_set, predicts_set, F=2)
    auc = calc_prec_rec_comp_pred(clusters, refs, {}, outfile_A, quiet=True) # also prints auc_points
    mmr = pc.maximum_matching_ratio(refs_set, predicts_set)
    sensitivity = pc.clusteringwise_sensitivity(refs_set, predicts_set)
    ppp = pc.positive_predictive_value(refs_set, predicts_set)
    accuracy = pc.accuracy(refs_set, predicts_set)
    f_match = pc.fraction_matched(refs_set, predicts_set)
    separation = pc.clusteringwise_separation(refs_set, predicts_set)
    
    # Collect Results and Append to results.csv Masterfile
    results = [method, gldstd, ppin, num_predicts, num_refs,
                       precision, recall, f1, f2, auc,
                       mmr, sensitivity, ppp, accuracy, f_match, separation]
    df.loc[len(df)] = results
    with outfile_M.open('a') as f:
        f.write(",".join([str(res) for res in results]))
        f.write("\n")
    
    print(df)
    