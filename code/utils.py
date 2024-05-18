import networkx as nx
import pandas as pd
from pathlib import Path

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
