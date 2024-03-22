import networkx as nx
import sys
from termcolor import colored

def printc(str, color="cyan", *args, **kwargs):
    "`print` wrapper that adds color to the console output."
    print(colored(str, color), end = kwargs.get('end', None))
    
def positive_int(x: str):
    """Returns whether the given string is a positive integer."""
    i = int(x)
    if i < 1:
        raise ValueError('Nonpositive values are not allowed')
    return i

def is_numeric(x: str):
    """Returns whether the given string can be interpreted as a number."""
    try:
        float(x)
        return True
    except:
        return False
    
def write_new_ppin(ppinpath: str, writepath: str, feat: str):
    """
    Writes a simple PPIN with lines `(p1 p2 score)` from 
    a composite PPIN with lines `(p1 p2 feat feat_score)`.
    Filters the protein pairs by `feat`.
    """
    with open(ppinpath) as f1, open(writepath, "w+") as f2:
        for line in f1:
            temp = line.split()
            if temp[2] == feat:
                f2.write("{} {} {}\n".format(temp[0], temp[1], temp[3]))
    
def write_new_complexfile(complexpath: str, writepath: str, with_cid=False):
    """
    Writes a complex file with lines `(p1 p2 ...)` from
    a reference file with lines `(pid cid cname)`.
    """
    complexes = dict()    
    with open(complexpath) as f1, open(writepath, "w+") as f2:
        for line in f1:
            temp = (line.split())[:2]
            pid, cid = temp[0], temp[1]
            if cid not in complexes.keys():
                complexes[cid] = [pid]
            else:
                complexes[cid].append(pid)
                    
        for cid in complexes:
            if with_cid: f2.write("{} {}\n".format(" ".join(complexes[cid])))
            else: f2.write("{}\n".format(" ".join(complexes[cid])))

def graph_memory(G: nx.Graph, units="MB"):
    """
    Returns graph memory usage in the specified byte units:
        units: [B|KB|MB] (Defaults to MB)
    """
    edge_mem = sum([sys.getsizeof(e) for e in G.edges])
    node_mem = sum([sys.getsizeof(n) for n in G.nodes])
    size = edge_mem + node_mem
    printer = lambda size, units: print("Total memory: %.6f %s" % (size, units))
    if units == "KB":
        printer(size/10**3, "KB")
    elif units == "MB":
        printer(size/10**6, "MB")
    else:
        printer(size, "bytes")

def graph_stats(G: nx.Graph, with_diameter=True):
    """Prints some graph statistics"""
    printc("Properties of Graph G")
    print("Size of V(G):", G.number_of_nodes())
    print("Size of E(G):", G.number_of_edges())
    print("Number of Components:", nx.number_connected_components(G))
    print("Average clustering coefficient:", round(nx.average_clustering(G), 2))
    if with_diameter:
        print("Max Diameter:", max(map(nx.diameter, [G.subgraph(c).copy() for c in nx.connected_components(G)])))

def graph_memory(G: nx.Graph, units="MB"):
    """
    Prints graph memory usage in the specified byte units:
        units: [B|KB|MB] (Defaults to MB)
    """
    edge_mem = sum([sys.getsizeof(e) for e in G.edges])
    node_mem = sum([sys.getsizeof(n) for n in G.nodes])
    size = edge_mem + node_mem
    printer = lambda size, units: print("Total memory: %.6f %s" % (size, units))
    if units == "KB":
        printer(size/10**3, "KB")
    elif units == "MB":
        printer(size/10**6, "MB")
    else:
        printer(size, "bytes")

def canonical_protein_name(name: str):
    """Returns the canonical name of a protein by performing a few simple
    transformations on the name."""
    return name.strip().upper()

def read_network(fname: str):
    """Returns a list of all the protein pairs `{p1, p2}` in a weighted or unweighted edge file."""
    known_proteins = list()
    for line in open(fname):
        parts = [canonical_protein_name(part) for part in line.strip().split() if not is_numeric(part)]
        known_proteins.append(set(parts))
    return known_proteins