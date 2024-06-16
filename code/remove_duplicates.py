# Program that removes duplicates in a predicted clusters file
# Do remove_duplicates.py --help for help on running this program.

import argparse
import PredictedClusters_Analysis as pc
import pandas as pd
from utils import printc
from classes import Cluster
from pathlib import Path

parser = argparse.ArgumentParser(description='Removes duplicate clusters from the predicted clusters file')
parser.add_argument('predictsfile', type=Path, help='path to predicted clusters file')
args = parser.parse_args()

# Function to read the clusters in a predicted clusters file
def read_Clusters(file):
    clusters = []
    with file.open() as f:
        for lineno, line in enumerate(f,1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            clus = Cluster(proteins, score, lineno)
            clus.matches = 0
            clusters.append(clus)

    return clusters

def match_clusters_lb(clus1, clus2):
    # Ensure clus1 is the larger cluster for efficiency
    if len(clus1) < len(clus2):
        clus1, clus2 = clus2, clus1
    
    # Compute the number of intersecting elements
    num_intersect = len(set(clus2).intersection(set(clus1)))
    
    # Compute the union of elements
    num_union = len(clus1) + len(clus2) - num_intersect
    
    return num_intersect / num_union

def remove_duplicate_clusters(clusters, match_thresh):
    clusters_to_keep = []
    clusters_to_delete = set()

    for idx1 in range(len(clusters)):
        if clusters[idx1].id in clusters_to_delete:
            continue
        for idx2 in range(idx1 + 1, len(clusters)):
            thresh = match_thresh
            clus1 = clusters[idx1]
            clus2 = clusters[idx2]
            if clus2.id in clusters_to_delete:
                continue
            if len(clus1.proteins) <= 3 and len(clus2.proteins) <= 3:
                thresh=1 # If small clusters, they should be equal
            matchscore = match_clusters_lb(clus1.proteins, clus2.proteins)
            if matchscore >= thresh:
                clusters_to_delete.add(clus2.id)
        
        clusters_to_keep.append(clusters[idx1])

    printc(f"In remove_duplicate_clusters, num to delete = {len(clusters_to_delete)}")

    return [clus for clus in clusters_to_keep if clus.id not in clusters_to_delete]

if __name__ == '__main__':
    predictsfile = Path(args.predictsfile)

    # Set parameters
    match_thresh = 0.75
    printc("Match threshold:\t%.2f" % match_thresh)

    clusters = read_Clusters(predictsfile)

    print("Number of clusters before removing duplicates:\t", len(clusters))
    clusters = remove_duplicate_clusters(clusters, match_thresh)
    print("Number of clusters after removing duplicates:\t", len(clusters))

    # Write clusters to the file again
    with predictsfile.open("w") as f:
        for cluster in clusters:
            # protein === node
            # Score the complex by their weighted density
            # Each line: (len(complex)_score): p1 p2 p3 ...
            score = cluster.score
            length = len(cluster.proteins)
            f.write(f"({length}_{score}): ")
            for protein in cluster.proteins:
                f.write("%s " % protein)
            f.write("\n")