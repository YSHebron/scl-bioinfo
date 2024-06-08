import networkx as nx
import argparse
import utils
from utils import printc
from pathlib import Path

parser = argparse.ArgumentParser(description="Iterative AdjustCD")
parser.add_argument('ppinfile', type=Path, help='path to pre-decomposistion PPIN')
parser.add_argument('outfile', type=Path, help='writepath for supplementary iAdjustCD PPIN')
args = parser.parse_args()

def size_of_all_neighborhoods_over_V(G: nx.Graph) -> float:
    num = sum([len(G[x]) for x in G.nodes])
    return num / G.number_of_nodes()

def _lambda(u: str, Nx: set, k: int, G: nx.Graph) -> float:
    global penalty_minuend
    if k == 1:
        return max(0, penalty_minuend - len(Nx))
    ret = 0
    for x in G.nodes:
        for y in G[x]:
            ret += iAdjustCD(x, y, k-1, G)
    penalty = ret/G.number_of_nodes() - sum([iAdjustCD(x, u, k-1, G) for x in Nx])
    return max(0, penalty)

def iAdjustCD(u: str, v: str, k: int, G: nx.Graph) -> float:
    Nu, Nv = set(G[u]), set(G[v])
    # Base case
    if k == 1:
        num = 2*len(Nu.intersection(Nv))
        den = len(Nu) + _lambda(u, Nu, 1, G) + len(Nv) + _lambda(v, Nv, 1, G)
        return num/den
    # For each node in the intersection of Nu and Nv, compute the ff
    num = sum([iAdjustCD(x, u, k-1, G) + iAdjustCD(x, v, k-1, G) for x in Nu.intersection(Nv)])
    den = sum([iAdjustCD(x, u, k-1, G) for x in Nu]) + _lambda(u, Nu, k-1, G) \
            + sum([iAdjustCD(x, v, k-1, G) for x in Nv]) + _lambda(v, Nv, k-1, G)
    if den == 0:    # BUG: unsure if denominator is allowed to be 0, consult iAdjustCD and CMC paper
        return 0
    return num/den

if __name__ == '__main__':
    printc("Performing iterative AdjustCD scoring...")
    
    G = utils.read_ppin_to_graph(args.ppinfile)
    global penalty_minuend
    penalty_minuend = size_of_all_neighborhoods_over_V(G)
    
    adjust_ppin = {}
    for u, v in G.edges:
        key = (u, v) if u < v else (v, u)
        adjust_ppin[key] = iAdjustCD(u, v, 2, G)
        
    utils.write_ppin_dict_to_txt(adjust_ppin, args.outfile, weighted=True)
    
