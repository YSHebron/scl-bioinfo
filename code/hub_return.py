import networkx as nx
import argparse
import utils
from utils import printc
from classes import Cluster
from pathlib import Path

parser = argparse.ArgumentParser(description="Hub Decomposition 2: Returning Hub Proteins")
parser.add_argument('predictsfile', type=Path, help='path to predicted clusters')
parser.add_argument('adjustfile', type=Path, help='path to iAdjustCD scored PPIN')
parser.add_argument('hubfile', type=Path, help='path to hub proteins')
parser.add_argument('ppinfile', type=Path, 
                    help='path to original possibly preprocessed PPIN, for regaining scores')
parser.add_argument('outfile', type=Path, help='writepath for postprocessed PPIN')
args = parser.parse_args()


def connectivity(u: str, C: Cluster, w: dict) -> float:
    ret = 0
    for v in C.proteins:
        key = (u, v) if u < v else (v, u)
        if key not in w:
            continue
        ret += w[key]    # where w(u, v) is calculated from iterative adjust CD
    return ret/len(C.proteins)

# function to score the cluster
def score_cluster(C: Cluster, w: dict):
    total = 0
    for u in C.proteins:
        for v in C.proteins:
            if u == v:
                continue
            key = (u, v) if u < v else (v, u)
            # BUG: hotfix for PC2P
            if key not in w:
                continue
            total += w[key]
    score = total*2 / ((len(C.proteins) * (len(C.proteins)-1)))
    return score

if __name__ == '__main__':
    # NOTE: Although P == clusters, C == complexes in our eqns, here C == clusters to follow lit.
    printc("Performing return of hub proteins...", "green")
    predicts = []
    adjust_ppin = utils.read_ppin_to_dict(args.adjustfile, weighted=True)
    hubs = []
    orig_ppin = utils.read_ppin_to_dict(args.ppinfile, weighted=True)
    with open(args.predictsfile) as f1, open(args.hubfile) as f2:
        cid = 1
        for line in f1:
            cluster = line.split()
            predicts.append(Cluster(cluster, id=cid))
            cid += 1
        hubs = f2.readline().split()
        
    hub_add_thres = 0.3     # literature-based
    
    num_hub_proteins_added = 0
    for u in hubs:
        for C in predicts:
            C: Cluster
            if connectivity(u, C, adjust_ppin) > hub_add_thres:
                C.proteins.add(u)
                num_hub_proteins_added += 1
                # print(f"Added hub {u} to Cluster {C.id}")
    printc(f"num_hub_proteins_added: {num_hub_proteins_added}")
                
    for C in predicts.copy():
        if len(C.proteins) < 2:
            predicts.remove(C)
            continue
        C.score = score_cluster(C, orig_ppin)
        if C.score == 0:
            predicts.remove(C)
                
    # Writing postprocessed clusters (with hubs)
    args.outfile.parent.mkdir(exist_ok=True, parents=True)
    with open(args.outfile, 'w') as f:
        for C in predicts:
            f.write(str(C))
            f.write("\n")