# Program that evaluates predicted complexes by getting the precision and recall
# Do eval.py --help for help on running this program.
# Sample run: python code/eval.py data/Results/Dummy/Trial/ClusterOne_postprocessed.txt data/Results/Dummy/Trial/CUBCO+_postprocessed.txt data/Results/Dummy/Trial/PC2P_postprocessed.txt data/Results/Dummy/Trial/P5COMP_clusters.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/Trial
# python code/eval.py data/Results/Dummy/Raw/Collins_CYC/ClusterOne_postprocessed.txt data/Results/Dummy/Raw/Collins_CYC/CUBCO+_postprocessed.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/Raw
# python code/eval.py data/Results/Dummy/P5COMP/ClusterOne_postprocessed.txt data/Results/Dummy/P5COMP/CUBCO+_postprocessed.txt data/Results/Dummy/P5COMP/PC2P_postprocessed.txt data/Results/Dummy/P5COMP/P5COMP_clusters.txt data/Yeast/CYC_complexes.txt data/Results/Dummy/P5COMP

import argparse
import PredictedClusters_Analysis as pc
from utils import printc
from classes import Cluster
from pathlib import Path

parser = argparse.ArgumentParser(description='Evaluate results of PC2P (or other clustering algorithm). Must follow line format (size_score): p1 p2 ...')
parser.add_argument('c1file', type=str, help='relpath to predicted clusters file by using CluserOne')
parser.add_argument('cubcofile', type=str, help='relpath to predicted clusters file by using CUBCO+')
parser.add_argument('pc2pfile', type=str, help='relpath to predicted clusters file by using PC2P')
parser.add_argument('p5compfile', type=str, help='relpath to predicted clusters file by using P5COMP')
parser.add_argument('complexfile', type=str, help='relpath to gold standard complex file')
parser.add_argument('outputdir', type=str, help='relpath to output dir for evaluation results')
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
    c1file, cubcofile, pc2pfile, p5compfile, complexfile, outputdir = Path(args.c1file), Path(args.cubcofile), Path(args.pc2pfile), Path(args.p5compfile), Path(args.complexfile), Path(args.outputdir)
    #c1file, cubcofile, complexfile, outputdir = Path(args.c1file), Path(args.cubcofile), Path(args.complexfile), Path(args.outputdir)
    printc("CluserOne clusters File:\t%s" % c1file)
    printc("CUBCO+ cluters File:\t%s" % cubcofile)
    printc("PC2P clusters File:\t%s" % pc2pfile)
    printc("P5COMP clusters File:\t%s" % p5compfile)
    printc("Complex File:\t%s" % complexfile)
    #printc("Output File:\t%s" % Path(outputdir, "P5COMP_eval.txt"))

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
    # Clusters are defined here as objects, with predicts as a set of clusters
    # Alternatively, predicts: { cid: { proteins: set(p1, p2, ...), score: float } }
    printc("Evaluating PC2P...\n")
    clusters_PC2P = []
    with pc2pfile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters_PC2P.append(cluster)
    clusters_PC2P.sort(key = lambda x: x.score, reverse=True)

    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    eval_PC2P = Path(outputdir, "PC2P_eval.txt")
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity
    calc_prec_rec_comp_pred(match_thresh, clusters_PC2P, refs, {}, eval_PC2P, quiet=True)
    
    print("Precision:\t%.6f" % pc.precision(refs, clusters_PC2P))
    print("Recall:\t\t%.6f" % pc.recall(refs, clusters_PC2P))
    print("F-score:\t%.6f" % pc.F_measure(refs, clusters_PC2P))

    printc("Evaluating CUBCO+...\n")
    clusters_CUBCO = []
    with cubcofile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters_CUBCO.append(cluster)
    clusters_CUBCO.sort(key = lambda x: x.score, reverse=True)

    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    eval_CUBCO = Path(outputdir, "CUBCO+_eval.txt")
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity
    calc_prec_rec_comp_pred(match_thresh, clusters_CUBCO, refs, {}, eval_CUBCO , quiet=True)
    
    print("Precision:\t%.6f" % pc.precision(refs, clusters_CUBCO))
    print("Recall:\t\t%.6f" % pc.recall(refs, clusters_CUBCO))
    print("F-score:\t%.6f" % pc.F_measure(refs, clusters_CUBCO))

    printc("Evaluating ClusterOne...\n")
    clusters_ClusterOne = []
    with c1file.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters_ClusterOne.append(cluster)
    clusters_ClusterOne.sort(key = lambda x: x.score, reverse=True)

    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    eval_ClusterOne = Path(outputdir, "ClusterOne_eval.txt")
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity
    calc_prec_rec_comp_pred(match_thresh, clusters_ClusterOne, refs, {}, eval_ClusterOne, quiet=True)
    
    print("Precision:\t%.6f" % pc.precision(refs, clusters_ClusterOne))
    print("Recall:\t\t%.6f" % pc.recall(refs, clusters_ClusterOne))
    print("F-score:\t%.6f" % pc.F_measure(refs, clusters_ClusterOne))

    printc("Evaluating P5COMP...\n")
    clusters_P5COMP = []
    with p5compfile.open() as f:
        for lineno, line in enumerate(f, 1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            cluster = Cluster(proteins, score, lineno)
            cluster.matches = 0
            clusters_P5COMP.append(cluster)
    clusters_P5COMP.sort(key = lambda x: x.score, reverse=True)

    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    eval_P5COMP = Path(outputdir, "P5COMP_eval.txt")
    
    ### Calculate Precision and Recall for Score Threshold s
    # precision = TP/(TP+FP) = TP/len(predicts), recall = TP/(TP+FN) = TP/len(gldstd)
    # Positive Predictive Value / Accuracy / Quality
    # True Positive Rate / Quantity
    calc_prec_rec_comp_pred(match_thresh, clusters_P5COMP, refs, {}, eval_P5COMP, quiet=True)
    
    print("Precision:\t%.6f" % pc.precision(refs, clusters_P5COMP))
    print("Recall:\t\t%.6f" % pc.recall(refs, clusters_P5COMP))
    print("F-score:\t%.6f" % pc.F_measure(refs, clusters_P5COMP))