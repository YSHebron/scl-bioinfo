import networkx as nx
import pandas as pd
import argparse
import utils
from pathlib import Path

parser = argparse.ArgumentParser(description="Hub Decomposition 1: Hub Removal")
parser.add_argument('ppinfile', type=Path, help='path to PPIN')
parser.add_argument('outfile', type=Path, help='writepath for decomposed PPIN')
parser.add_argument('hubfile', type=Path, help='writepath for hub proteins')
args = parser.parse_args()

if __name__ == '__main__':
    G = utils.read_ppin_to_graph(args.ppinfile)
    print(sorted(list(G.degree), key=lambda a:a[1]))