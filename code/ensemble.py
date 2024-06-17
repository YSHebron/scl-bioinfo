# Program that combines all of the predicted clusters
# Do ensemble.py --help for help on running this program.
# Sample run: python code/ensemble.py data/Results/Dummy/Trial/ClusterOne_postprocessed.txt data/Results/Dummy/Trial/CUBCO+_postprocessed.txt data/Results/Dummy/Trial/PC2P_postprocessed.txt data/Results/Dummy/Trial
# python code/ensemble.py data/Results/Dummy/P5COMP/ClusterOne_postprocessed.txt data/Results/Dummy/P5COMP/CUBCO+_postprocessed.txt data/Results/Dummy/P5COMP/PC2P_postprocessed.txt data/Results/Dummy/P5COMP/P5COMP_clusters.txt


import argparse
from utils import printc
from classes import Cluster
from pathlib import Path

parser = argparse.ArgumentParser(description='Combining the clusters. Must follow line format (size_score): p1 p2 ...')
parser.add_argument('c1file', type=str, help='relpath to predicted clusters file by using CluserOne')
parser.add_argument('cubcofile', type=str, help='relpath to predicted clusters file by using CUBCO+')
parser.add_argument('pc2pfile', type=str, help='relpath to predicted clusters file by using PC2P')
parser.add_argument('outfile', type=str, help='relpath to output dir for evaluation results')
args = parser.parse_args()


# Function to read the clusters in a predicted clusters file
def read_Clusters(file, initial_lineno):
    clusters = []
    with file.open() as f:
        for lineno, line in enumerate(f,1):
            # Let lineno be the cluster id
            raw = line.split()
            proteins, score = raw[1:], float(raw[0].split("_")[-1].strip("):"))
            clus = Cluster(proteins, score, lineno+initial_lineno)
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
            if matchscore >= 1:
                clusters_to_delete.add(clus2.id)
        
        clusters_to_keep.append(clusters[idx1])

    printc(f"In remove_duplicate_clusters, num to delete = {len(clusters_to_delete)}")

    return [clus for clus in clusters_to_keep if clus.id not in clusters_to_delete]

def remove_duplicate_clusters_id(clusters):
    unique_clusters = {}
    for cluster in clusters:
        if cluster.id not in unique_clusters:
            unique_clusters[cluster.id] = cluster
        # else:
        #     print(f"Duplicate ID found and removed: {cluster.id}")
    
    printc(f"In remove_duplicate_clusters_id, num to delete = {len(clusters) - len(unique_clusters)}")
    return list(unique_clusters.values())

def ensemble_clustering(cluster1, cluster2, cluster3, match_thresh):
    clusters_to_keep = []

    # Check if a cluster in cluster1 is in cluster2 and cluster3 based on a match threshold (match_thresh)
    # If match_clusters_lb >= match_thresh, it is a match and store the cluster with the higher score in clusters_to_keep

    for index1 in range(len(cluster1)):
        votes = [] # This stores the similar clusters.
        clus1 = cluster1[index1]
        votes.append(clus1)

        # Check if the cluster is in cluster2
        for index2 in range(len(cluster2)):
            clus2 = cluster2[index2]
            matchscore = match_clusters_lb(clus1.proteins, clus2.proteins)
            if matchscore >= match_thresh:
                # It is in cluster2. +1 to the vote
                votes.append(clus2)
                # Done checking cluster 2
                break
        
        # Check if the cluster is in cluster3
        for index3 in range(len(cluster3)):
            clus3 = cluster3[index3]
            matchscore = match_clusters_lb(clus1.proteins, clus3.proteins)
            if matchscore >= match_thresh:
                # It is in cluster3. +1 to the vote
                votes.append(clus3)
                # Done checking cluster 3
                break

        # Check how many votes. If >=2 then the cluster with the highest score will be kept
        if len(votes) >= 2:
            # Sort votes by score in descending order and keep the highest score cluster
            votes.sort(key=lambda x: x.score, reverse=True)
            clusters_to_keep.append(votes[0])

    print("In ensemble_clustering, num to keep:\t", len(clusters_to_keep))

    return clusters_to_keep


def remove_low_reliability_clusters(clusters):
    clusters_to_keep = []
    for cluster in clusters:
        cluster: Cluster
        if cluster.score > 1.0:
            clusters_to_keep.append(cluster)
    return clusters_to_keep


if __name__ == '__main__':
    c1file, cubcofile, pc2pfile, outfile = Path(args.c1file), Path(args.cubcofile), Path(args.pc2pfile), Path(args.outfile)
    printc("CluserOne clusters File:\t%s" % c1file)
    printc("CUBCO+ cluters File:\t%s" % cubcofile)
    printc("PC2P clusters File:\t%s" % pc2pfile)
    printc("Output File:\t%s" % outfile)
    
    # Set parameters
    match_thresh = 0.75
    printc("Match threshold:\t%.2f" % match_thresh)

    c1_clusters = read_Clusters(c1file, 0)
    # print("ClusterOne clusters before removing duplicates:\t", len(c1_clusters))
    # c1_clusters = remove_duplicate_clusters(c1_clusters, match_thresh)
    # print("ClusterOne clusters after removing duplicates:\t", len(c1_clusters))

    cubco_clusters = read_Clusters(cubcofile, len(c1_clusters))
    print("CUBCO+ clusters before removing duplicates:\t", len(cubco_clusters))
    cubco_clusters = remove_duplicate_clusters(cubco_clusters, match_thresh)
    print("CUBCO+ clusters after removing duplicates:\t", len(cubco_clusters))

    pc2p_clusters = read_Clusters(pc2pfile, len(c1_clusters)+len(cubco_clusters))
    # print("PC2P clusters before removing duplicates:\t", len(pc2p_clusters))
    # pc2p_clusters = remove_duplicate_clusters(pc2p_clusters, match_thresh)
    # print("PC2P clusters after removing duplicates:\t", len(pc2p_clusters))
    
    combined_clusters = c1_clusters + cubco_clusters + pc2p_clusters
    printc("Total number of clusters:\t%d" %len(combined_clusters))
    # NOTE: Found out the best configuration is a simple concatenation of all the results
    # final_clusters = ensemble_clustering(c1_clusters, cubco_clusters, pc2p_clusters, match_thresh)
    # printc("Number of possible clusters:\t%d" %len(final_clusters))
    # final_clusters += ensemble_clustering(cubco_clusters, pc2p_clusters, c1_clusters, match_thresh)
    # printc("Number of possible clusters:\t%d" %len(final_clusters))
    # final_clusters += ensemble_clustering(pc2p_clusters, c1_clusters, cubco_clusters, match_thresh)
    # printc("Number of possible clusters:\t%d" %len(final_clusters))
    # final_clusters = remove_duplicate_clusters_id(final_clusters)
    # final_clusters = remove_low_reliability_clusters(combined_clusters)
    # print("P5COMP clusters:\t", len(final_clusters))

    # Output results
    with outfile.open('w+') as f:
        for cluster in combined_clusters:
            f.write(str(cluster) + "\n")

    


