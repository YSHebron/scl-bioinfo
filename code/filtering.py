import argparse
import utils
import networkx as nx
from utils import printc
from pathlib import Path

parser = argparse.ArgumentParser(description='Gold standard filtering and optional Negatome filtering.')
parser.add_argument('ppinfile', type=Path, help='path to PPIN')
parser.add_argument('reffile', type=Path, help='path to gold standard')
parser.add_argument('outfile', type=Path, help='writepath for filtered PPIN')
parser.add_argument('--negfile', type=Path, help='path to negatome', nargs='?', required=False, const=None)
parser.add_argument('--confidence', type=float, metavar='[0.0-1.0]', default=0.0, help='interaction score threshold')
parser.add_argument('--filtering', type=str, choices='[perpair,perprotein,direct]', default='perpair', help='initial filtering')
parser.add_argument('--directfile', type=Path, help='required when using direct filtering', nargs='?', required=False)
args = parser.parse_args()


def filter_perpair(ppin: dict, reffile: Path) -> dict:
    # NOTE: May improve detection of small complexes, but worsen detection of sparse and overlapping complexes. [Verify]
    pairs_in_ref = set()    # NOTE: Unrealistically assumes that all proteins in a complex interact
    with open(reffile) as f:
        for line in f:
            ref_complex = line.split()
            for u in ref_complex:
                for v in ref_complex:
                    if u == v:
                        continue
                    pair = (u, v) if u < v else (v, u)
                    pairs_in_ref.add(pair)
            
    # Filtering
    for ppi in ppin.copy():
        if ppi not in pairs_in_ref or ppin[ppi] < args.confidence:
            del ppin[ppi]
    
    return ppin
            

def filter_perprotein(ppin: dict, reffile: Path) -> dict:
    # NOTE: May improve detection of all complexes, but not as much as perpair.
    present = set()
    with open(reffile) as f:
        for line in f:
            complex = line.split()
            for pid in complex:
                present.add(pid)
    
    # Filtering
    for ppi in ppin.copy():
        if not(ppi[0] in present and ppi[1] in present) or ppin[ppi] < args.confidence:
            del ppin[ppi]
        
    return ppin


def filter_direct(ppin: nx.Graph, directfile: nx.Graph) -> nx.Graph:
    """
    Concept: Take the unweighted edges from directfile and score them using ppin
    Effectively filters ppin to only the edges present in directfile
    """
    for e in directfile.edges:
        directfile[e[0]][e[1]]['weight'] = ppin[e[0]][e[1]]['weight']
    return directfile
    

def filter_with_negatome(ppin: dict, negppin: dict) -> dict:
    filtered_ppin = {}
    for ppi in ppin.keys() - negppin.keys():
        filtered_ppin[ppi] = ppin[ppi]
    
    return filtered_ppin

def reffile_to_graph(reffile: Path) -> Path:
    # y1 y2 y3
    # y1 y2
    # -->
    # y1 y2
    # y1 y3
    edgelist = set()
    with reffile.open() as f:
        for line in f:
            complex = line.strip(' \n').split()
            for i, ui in enumerate(complex):
                for j, uj in enumerate(complex):
                    if i == j: continue
                    edgelist.add(frozenset([ui, uj]))
    outfile = Path("data/Interm/reffile_edge_list.txt")
    with outfile.open('w') as f:
        for edge in edgelist:
            e = list(edge)
            f.write(f"{e[0]} {e[1]}")
            f.write('\n')
    return(outfile)

if __name__ == '__main__':
    # TODO: Add utils that could read_ppin_to_graph.
    # Might be easier to work with this representation.
    ppin = utils.read_ppin_to_dict(args.ppinfile, weighted=True)
    printc("Filtering applied\tsize(E)")
    print(f"None\t\t\t{len(ppin)}")
    
    # Filtering Directives
    if args.negfile is not None:
        ppin = filter_with_negatome(ppin, utils.read_ppin_to_dict(args.negfile))
        print(f"Negatome 2.0\t\t{len(ppin)}")
        
    if args.filtering == "perpair":
        ppin = filter_perpair(ppin, args.reffile)
        print(f"Per Pair\t\t{len(ppin)}")
    elif args.filtering == "perprotein":
        ppin = filter_perprotein(ppin, args.reffile)
        print(f"Per Protein\t\t{len(ppin)}")
    elif args.filtering == "direct" and args.directfile is None:
        raise FileNotFoundError("Direct file is required when using direct filtering.")
    elif args.filtering == "direct" and args.directfile:
        ppin = nx.read_weighted_edgelist(args.ppinfile, create_using=nx.Graph)
        direct = nx.read_edgelist(args.directfile, create_using=nx.Graph)
        # reffile = nx.read_edgelist(reffile_to_graph(args.reffile), create_using=nx.Graph)
        direct = filter_direct(ppin, direct)
        print(f"Direct\t\t\t{len(direct.edges)}")
        nx.write_weighted_edgelist(direct, args.outfile)

    if args.filtering != "direct":
        utils.write_ppin_dict_to_txt(ppin, args.outfile, weighted=True)
    