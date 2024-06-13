# Do cluster_one_scoring.py --help for help on running this program.
# Sample run: python code/ClusterOne/cluster_one_scoring.py data/Yeast/FilteredPPINs/Collins_CYC_perpair_weighted.txt data/Results/Dummy/clusters_clusterone.py

import networkx as nx
import pandas as pd
import argparse
from pathlib import Path
from utils import printc

parser = argparse.ArgumentParser(description='Perform PC2P on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('inputfile', type=Path, help='relpath to PPI dataset to be processed')
parser.add_argument('clusterfile', type=Path, help='relpath to clusters file to be scored')
parser.add_argument('outfile', type=Path, help='writepath for ClusterOne scored clusters' )
args = parser.parse_args()

# function to score the complex
def get_score(complex, edgesref):
    totalweight = 0
    for id1 in range(len(complex)):
        for id2 in range(1,len(complex)):
           id_a = complex[id1]
           id_b = complex[id2]
           key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
           if key in edgesref:
               totalweight += edgesref[key]
    score = totalweight*2 / ((len(complex) * (len(complex)-1)))
    return score


if __name__ == '__main__':
    inputfile, clusterfile, outfile = args.inputfile, args.clusterfile, args.outfile
    
    
    ### Read the PPIN (given as a PPI dataset or edgelist)
    # Assumes the inputfile has scores, which may or may not be used
    G = nx.Graph()
    if inputfile.suffix == ".csv":
        # for .csv with header p1, p2, and weight
        df = pd.read_csv(inputfile)
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.suffix == ".tsv":
        # for .tsv with header p1, p2, and weight
        df = pd.read_csv(inputfile, sep="\t")
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.suffix == ".txt":
        # for .txt with no header edge lists
        G = nx.read_weighted_edgelist(inputfile, create_using = nx.Graph, nodetype = str)
    
    printc("CWD:\t%s " % Path.cwd())
    printc("Input:\t%s" % inputfile)
    printc("Cluster file:\t%s" % clusterfile)

    ### Prepare for Scoring Detected Clusters (done simultaneously with writing to result)
    scorededges = {}
    for id_a, id_b, score in G.edges(data=True):
        # keys are tuples of protein pairs, and values are their scores
        key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
        scorededges[key] = score["weight"]
    
    with clusterfile.open("r") as f:
        lines = f.readlines()

    # Process each line to remove whitespace and split by tab characters
    cmplx = [line.strip().split('\t') for line in lines if line.strip()]
    cmplx.sort(key=len, reverse=True)

    # print(cmplx)
    outfile.parent.mkdir(exist_ok=True, parents=True)
    with outfile.open("w") as f:
        # complex === line
        for complex in cmplx:
            # protein === node
            # Score the complex by their weighted density
            # Each line: (len(complex)_score): p1 p2 p3 ...
            score = get_score(complex, scorededges)
            f.write(f"({len(complex)}_{score}): ")
            for protein in complex:
                f.write("%s " % protein)
            f.write("\n")

    