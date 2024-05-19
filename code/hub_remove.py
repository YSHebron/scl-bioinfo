import networkx as nx
import time
import argparse
import utils
import matplotlib.pyplot as plt
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('ppinfile', type=Path)
args = parser.parse_args()

if __name__ == '__main__':
    ppinfile = args.ppinfile
    G = utils.read_ppin_to_graph(ppinfile)
    
    