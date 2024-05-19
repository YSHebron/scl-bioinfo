import networkx as nx
import argparse
import utils
from utils import printc
from pathlib import Path

parser = argparse.ArgumentParser(description="Hub Decomposition 1: Hub Removal")
parser.add_argument('ppinfile', type=Path, help='path to PPIN')
parser.add_argument('outfile', type=Path, help='writepath for decomposed PPIN')
parser.add_argument('hubfile', type=Path, help='writepath for hub proteins')
args = parser.parse_args()

def rel_connectivity(G: nx.Graph):
    largest_cc = max(nx.connected_components(G), key=len)
    H = G.subgraph(largest_cc)
    f = H.number_of_nodes() / G.number_of_nodes()
    return f

if __name__ == '__main__':
    G = utils.read_ppin_to_graph(args.ppinfile)
    
    printc("Performing Hub Removal...")
    
    # Topologically identify nhub
    ranked_list = []
    for item in G.degree:
        ranked_list.append(item)
    ranked_list = sorted(ranked_list, key=lambda a:a[1], reverse=True)
    ranked_nodes = [node[0] for node in ranked_list]
    
    print(f"Relative Connectivity of PPIN: {rel_connectivity(G)}")
    
    sharp_increase = 0
    deg_at_sharp = 0
    f_prev = 1.0
    
    Gn = nx.Graph()
    for n in range(len(ranked_list)):
        nbunch = ranked_nodes[0:n+1]
        Gn = G.subgraph(nbunch)
        f = rel_connectivity(Gn)
        df = f - f_prev
        # print(df)
        if df > 0 and df > sharp_increase:
            sharp_increase = df
            deg_at_sharp = G.degree(ranked_nodes[n])
        f_prev = f
    
    print(f"Largest df:\t{sharp_increase}")
    print(f"nhub cutoff:\t{deg_at_sharp}")
    nhub = deg_at_sharp
    
    # Filter out hub proteins from ranked list
    hubs = [item[0] for item in filter(lambda node:node[1] >= nhub, ranked_list)]
    
    # Write hub proteins to hubfile
    args.hubfile.parent.mkdir(exist_ok=True, parents=True)
    with open(args.hubfile, 'w') as f:
        for node in hubs:
            f.write(f"{node} ")
        f.write('\n')

    # Write decomposed PPIN to decompfile
    G.remove_nodes_from(hubs)
    args.outfile.parent.mkdir(exist_ok=True, parents=True)    
    nx.write_weighted_edgelist(G, args.outfile)
    print(f"Decomposed PPIN size(E): {G.number_of_edges()}")
        