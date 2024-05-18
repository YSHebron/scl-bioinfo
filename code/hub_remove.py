import networkx as nx
import time
import argparse
import utils
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('ppinfile', type=Path)
parser.add_argument('--nhub', type=int)
args = parser.parse_args()

if __name__ == '__main__':
    ppinfile, nhub = args.ppinfile, args.nhub
    G = utils.read_ppin_to_graph(ppinfile)
    
    