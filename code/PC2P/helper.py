import networkx as nx
import sys
from termcolor import colored

def printc(str, *args, **kwargs):
    print(colored(str, "red"), end = kwargs.get('end', None))
    
def positive_int(x):
    i = int(x)
    if i < 1:
        raise ValueError('Nonpositive values are not allowed')
    return i
    
# Tranforms PPIN from Integrated to something readable by our PC2P code
def write_new_ppin(ppin, ppin_new):
    with open(ppin, "r") as f1, open(ppin_new, "w+") as f2:
        for line in f1:
            temp = line.split()
            if temp[2] == "PPIREL":
                f2.write("{} {} {}\n".format(temp[0], temp[1], temp[3]))
    
# Tranforms Gold Standard from Integrated to something readable by our PC2P code
def write_new_complexfile(complexfile, complexfile_new):
    complexes = dict()
    
    ### Read complexes
    with open(complexfile) as f:
        for line in f:
            temp = (line.split())[:2]
            pid, cid = temp[0], temp[1]
            if cid not in complexes.keys():
                complexes[cid] = [pid]
            else:
                complexes[cid].append(pid)
    
    print(complexes)
    
    with open(complexfile_new, "w") as f:
        for cid in complexes:
            f.write("{0} {1}\n".format(cid, ' '.join(complexes[cid])))
            
def graph_stats(G: nx.Graph):
    printc("Properties of G")
    print("Size of V(G):", G.number_of_nodes())
    print("Size of E(G):", G.number_of_edges())
    print("Number of Components:", nx.number_connected_components(G))
    print("Average clustering coefficient:", round(nx.average_clustering(G), 2))
    print("Max Diameter:", max(map(nx.diameter, [G.subgraph(c).copy() for c in nx.connected_components(G)])))
    
def graph_memory(G: nx.Graph):
    edge_mem = sum([sys.getsizeof(e) for e in G.edges])
    node_mem = sum([sys.getsizeof(n) for n in G.nodes])

    print("Edge memory:", edge_mem)
    print("Node memory:", node_mem)
    print("Total memory:", edge_mem + node_mem)