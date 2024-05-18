import networkx as nx
import pandas as pd
from pathlib import Path
import sys
import os
from termcolor import colored


def printc(str, color="cyan", *args, **kwargs):
    "`print` wrapper that adds color to the console output."
    print(colored(str, color), end=kwargs.get("end", None))


def positive_int(x: str):
    """Returns whether the given string is a positive integer."""
    i = int(x)
    if i < 1:
        raise ValueError("Nonpositive values are not allowed")
    return i


def is_numeric(x: str):
    """Returns whether the given string can be interpreted as a number."""
    if isinstance(x, float):
        return True
    
def read_ppin_to_graph(ppinfile: Path) -> nx.Graph:
    """
    Reads an input PPIN and converts it to a NetworkX Graph.
    Recommended: .txt without header where each row is (u v w).
    Allowed: .csv/.tsv with header (u, v, w) where each row is (u, v, w).

    Args:
        ppinfile (Path): path to PPIN
    """

    if ppinfile.suffix != ".txt":
        df = pd.read_csv(ppinfile)
        return nx.from_pandas_edgelist(df, source = "u", target = "v", create_using = nx.Graph, edge_attr = "w")
    else:
        return nx.read_weighted_edgelist(ppinfile, create_using = nx.Graph, nodetype = str)


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
            if with_cid:
                f2.write("{}\n".format(" ".join(complexes[cid])))
            else:
                f2.write("{}\n".format(" ".join(complexes[cid])))


def graph_memory(G: nx.Graph, units="MB"):
    """
    Returns graph memory usage in the specified byte units:
        units: [B|KB|MB] (Defaults to MB)
    """
    edge_mem = sum([sys.getsizeof(e) for e in G.edges])
    node_mem = sum([sys.getsizeof(n) for n in G.nodes])
    size = edge_mem + node_mem

    def printer(size, units):
        print("Total memory: %.6f %s" % (size, units))

    if units == "KB":
        printer(size / 10**3, "KB")
    elif units == "MB":
        printer(size / 10**6, "MB")
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
        print(
            "Max Diameter:",
            max(
                map(
                    nx.diameter,
                    [G.subgraph(c).copy() for c in nx.connected_components(G)],
                )
            ),
        )


def canonical_protein_name(name: str):
    """Returns the canonical name of a protein by performing a few simple
    transformations on the name."""
    return name.strip().upper()


def read_network(fname: str):
    """Returns a list of all the protein pairs `{p1, p2}` in a weighted or unweighted edge file."""
    known_proteins = list()
    for line in open(fname):
        parts = [
            canonical_protein_name(part)
            for part in line.strip().split()
            if not is_numeric(part)
        ]
        known_proteins.append(set(parts))
    return known_proteins


def filterPerProtein(ppinname, ppinfile, gldstd, outputdir, species="Yeast"):
    filtering = "perprotein"
    rawPPIN = "code/PC2P/{}/{}/{}".format(species, ppinname, ppinfile)
    reffile = "code/PC2P/{}/{}_complexes.txt".format(species, gldstd)
    outputfile = outputdir + "/{}_{}_{}_weighted.txt".format(
        ppinname, gldstd, filtering
    )
    present = set()
    scorededges = {}
    with open(reffile) as f1, open(rawPPIN) as f2:
        for line in f1:
            for pid in line.split():
                present.add(pid)
        for line in f2:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    with open(outputfile, "w+") as f:
        for key in scorededges:
            if key[0] in present or key[1] in present:
                f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))
                
                
def filterDirect(ppinname, ppinfile, gldstd, outputdir, species="Yeast"):
    filtering = "direct"
    rawPPIN = "code/PC2P/{}/{}/{}".format(species, ppinname, ppinfile)
    directfile = "code/PC2P/{}/{}/{}_{}_Graph.txt".format(
        species, ppinname, ppinname, gldstd
    )
    outputfile = outputdir + "/{}_{}_{}_weighted.txt".format(
        ppinname, gldstd, filtering
    )
    scorededges = {}
    with open(rawPPIN) as f1:
        for line in f1:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)
    with open(directfile) as interm, open(outputfile, "w+") as f:
        for line in interm:
            p1, p2 = line.split()
            linepair = (p1, p2) if p1 < p2 else (p2, p1)
            for key in scorededges:
                if linepair == key:
                    f.write("%s %s %f\n" % (linepair[0], linepair[1], scorededges[key]))