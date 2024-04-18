# Program that evaluates predicted complexes by getting the precision and recall
# Do PC2P_eval.py --help for help on running this program.
# Sample run: python code/PC2P/PC2P_eval.py data/Results/Dummy/Dummy_CYC_testonly_predicted.txt data/Yeast/CYC_complexes.txt data/Analysis/Dummy

import argparse
import PredictedClusters_Analysis as pc
from helper import printc
from pathlib import Path

parser = argparse.ArgumentParser(description='Evaluate results of PC2P (or other clustering algorithm). Must follow line format (size_score): p1 p2 ...')
parser.add_argument('predictsfile', type=str, help='relpath to predicted clusters file')
parser.add_argument('complexfile', type=str, help='relpath to gold standard complex file')
parser.add_argument('outputdir', type=str, help='relpath to output dir for evaluation results')
args = parser.parse_args()

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

def calc_prec_rec_comp_pred(matchscore_thr, clusters, refs, matched_complexes_ref, outputfile: Path, quiet=True):
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

    # print("Num predicted clusters = ", len(clusters))
    # print("\tNum correct clusters = ", len(correct_clusters))
    # print("\tNum correct small clusters = ", len(correct_smallclusters))
    # print("\tNum complexes matched = ", len(matched_complexes_ref))
    # print("\tNum clusters not correct = ", len(clusters_copy))

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

    f = outputfile.open("w")
    filtering = outputfile.stem.split("_")[2]
    filename = outputfile.name
    f.write(f"{filtering}\n{filename}\n")
    if not quiet:
        print("Threshold\tPreds\tPrecision\tRecall")
    for cluster in sorted(results, key=lambda x: x[0], reverse=True):
        if cluster[0] != score_thresh:
            recall = matched/len(refs)   # Without train-test split, all our refs are technically test complexes
            if (score_thresh != -1 and recall > rec_thresh):
                precision = corrects/predicts
                if not quiet:
                    print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, predicts, precision, recall))
                f.write("%.6f\t%d\t%.6f\t%.6f\n" % (score_thresh, predicts, precision, recall))
                
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
        f.write("%.6f\t%d\t%.6f\t%.6f\n" % (score_thresh, predicts, corrects/predicts, recall))
    else:
        recall = matched / len(refs)
        precision = 0
        if not quiet:
            print("%.6f\t%d\t%.6f\t%.6f" % (score_thresh, predicts, corrects/predicts, recall))
        f.write("%.6f\t%d\t%.6f\t%.6f\n" % (score_thresh, predicts, corrects/predicts, recall))
    f.close()

    auc += (recall - auc_prevrecall) * (corrects/predicts)
    print("AUC:\t\t%.6f" % auc)


if __name__ == '__main__':
    predictsfile, complexfile, outputdir = Path(args.predictsfile), Path(args.complexfile), Path(args.outputdir)
    printc("Predicts File:\t%s" % predictsfile)
    printc("Complex File:\t%s" % complexfile)
    printc("Output File:\t%s" % Path(outputdir, predictsfile.stem.removesuffix("_predicted") + "_eval.txt"))
    # Set parameters
    match_thresh = 0.5

    # refs: reference complexes in the gold standard
    refs = []
    with complexfile.open() as f:
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
    with predictsfile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters.append(cluster)
    clusters.sort(key = lambda x: x.score, reverse=True)

    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    outputfile = Path(outputdir, predictsfile.stem.removesuffix("_predicted") + "_eval.txt")
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity
    calc_prec_rec_comp_pred(match_thresh, clusters, refs, {}, outputfile, quiet=True)
    
    print("Precision:\t%.6f" % pc.precision_Jaccard(refs, clusters))
    print("Recall:\t\t%.6f" % pc.recall_Jaccard(refs, clusters))
    print("F-score:\t%.6f" % pc.F_measure_Jaccard(refs, clusters))
        
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity