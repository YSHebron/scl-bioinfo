import networkx as nx
import pandas as pd
from pathlib import Path
from termcolor import colored

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

    
def read_ppin_to_dict(ppinfile: Path, weighted=False) -> dict:
    ppin = {}
    with open(ppinfile) as f:
        if weighted:
            for line in f:
                u, v, s = line.split()
                key = (u, v) if u < v else (v, u) # lexical ordering
                ppin[key] = float(s)
        else:
            for line in f:
                u, v = line.split()
                key = (u, v) if u < v else (v, u) # lexical ordering
                ppin[key] = None
            
    return ppin


def write_ppin_dict_to_txt(ppin: dict, outfile: Path, weighted=False):
    outfile.parent.mkdir(exist_ok=True, parents=True)
    with open(outfile, 'w') as f:
        if weighted:
            for ppi in ppin:
                u, v = ppi
                s = ppin[ppi]
                f.write(f"{u} {v} {s}\n")
        else:
            for ppi in ppin:
                u, v = ppi
                f.write(f"{u} {v}\n")


def printc(str, color="cyan", *args, **kwargs):
    "`print` wrapper that adds color to the console output."
    print(colored(str, color), end=kwargs.get("end", None))