# Do CUBCO.py --help for help on running this program.
# Sample run: python code/CUBCO+/CUBCO.py data/Yeast/FilteredPPINs/Collins_CYC_perpair_weighted.txt data/Results
# Sample run: python code/CUBCO+/CUBCO.py data/Dummy_CYC_testonly.txt data/Results 

import networkx as nx
import math
import pandas as pd
import argparse
import time
from pathlib import Path
from utils import printc, positive_int, graph_stats, graph_memory
from operator import itemgetter

parser = argparse.ArgumentParser(description='Perform CUBCO+ on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('inputfile', type=Path, help='relpath to PPI dataset to be processed')
parser.add_argument('outputdir', type=Path, help='relpath to output dir for predicted complexes')
args = parser.parse_args()

def ave_number_of_p3(e,G):
    ws = [p for p in nx.all_simple_paths(G,e[0],e[1],cutoff=3) if len(p)==4]
    s = 0
    if ws:
        for w in ws:
            s += 1/(math.sqrt(nx.degree(G,w[1])*nx.degree(G,w[2])))
    return s

def yield_cmp_len1(*args, **kwargs):
    global round
    v = set(G_bar.nodes()).difference(G_bar[reachable[0]])
    G_bar.remove_nodes_from(v)
    yield v

def yield_cmp_len1_nonreach(*args, **kwargs):
    global round
    v = set(G_bar.nodes()).difference(G_bar[reachable[0]])
    G_bar.remove_nodes_from(v)
    yield v

def score(H,cmp_nodes):
    degree = H.degree(cmp_nodes)
    sum_dg1 = sum([d[1] for d in degree ])

    h = nx.induced_subgraph(H,cmp_nodes)
    degree = h.degree(cmp_nodes)
    sum_dg2 = sum([d[1] for d in degree ])
    
    return (sum_dg2/2)/(sum_dg1-sum_dg2)

def yield_cmp_len_greater_1(*args, **kwargs):
    global round
    f1 = False
    
    node_set1 = set()
    node_set1.update(nbr for n in reachable for nbr in G_bar[n] if nbr in non_reachable)
    if node_set1 == set(non_reachable):
        f1 = True
    else:
        cmp_nodes1 = G_bar.nodes() ^ node_set1
        score1 = score(G,cmp_nodes1)
    
    node_set2 = set()
    node_set2.update(nbr for n in non_reachable for nbr in G_bar[n] if nbr in reachable)
    if node_set2 == set(reachable):
        if f1:
            degree = G_bar.degree()
            sorted_nodes = sorted(degree, key=itemgetter(1))
            mindeg = sorted_nodes[0][0]
            v = set(G_bar.nodes()).difference(G_bar[mindeg])
            #v = set(G[mindeg])
            #v.add(mindeg)
            G_bar.remove_nodes_from(v)
            yield v
        else:
            G_bar.remove_nodes_from(cmp_nodes1)
            yield cmp_nodes1
        
    elif f1:
        cmp_nodes2 = G_bar.nodes() ^ node_set2
        G_bar.remove_nodes_from(cmp_nodes2)
        yield cmp_nodes2
    else:
        cmp_nodes2 = G_bar.nodes() ^ node_set2
        score2 = score(G,cmp_nodes2)
        if score1 > score2:
            G_bar.remove_nodes_from(cmp_nodes1)
            yield cmp_nodes1
        else:
            G_bar.remove_nodes_from(cmp_nodes2)
            yield cmp_nodes2

def min_cut(*args, **kwargs):
    global round
    while len(G_bar.nodes()) > 1:
        if len(G_bar.nodes()) <= 3:
            yield set(G_bar.nodes())
            break
        elif not nx.is_connected(G_bar):
            yield set(G_bar.nodes())
            G_bar.clear()
        else:
            global reachable, non_reachable
            cut_value, partition = nx.stoer_wagner(G_bar)
            reachable, non_reachable = partition
            if len(reachable) == 1:
                yield yield_cmp_len1().__next__()
            elif len(non_reachable) == 1:
                yield yield_cmp_len1_nonreach().__next__()
            else:
                yield yield_cmp_len_greater_1().__next__()
    else:
        yield set(G_bar.nodes())

if __name__ == '__main__':
    inputfile, outputdir = args.inputfile, args.outputdir
    global reachable, non_reachable, G_bar, G, node_deg1
    start_time = time.time()
    rounds = 1
    
    ### Read the PPIN (given as a PPI dataset or edgelist)
    # Assumes the inputfile has scores, which may or may not be used
    G = nx.Graph()
    if inputfile.suffix == ".csv":
        # for .csv with header p1, p2, and weight
        df = pd.read_csv(inputfile)
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.suffix == ".tsv":
        # for .tsv with header p1, p2, and weight
        df = pd.read_csv(inputfile, sep="\t")
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.suffix == ".txt":
        # for .txt with no header edge lists
        G = nx.read_weighted_edgelist(inputfile, create_using = nx.Graph, nodetype = str)
    
    outputfile = Path(outputdir, "clusters_cubcoplus.txt")
    gcomplement = Path(outputdir, inputfile.stem.removesuffix("_weighted") + "_complement.txt")
    printc("CWD:\t%s " % Path.cwd())
    printc("Input:\t%s" % inputfile)
    printc("G Complement: \t%s" % gcomplement)
    printc("Output:\t%s" % outputfile)
    graph_stats(G)
    graph_memory(G)
    print(nx.to_numpy_array(G))

    ### PreProcessing
    new_edges = list()
    printc('Getting the complement of G ......\n')
    for g in nx.connected_components(G):
        g_cmp = nx.complement(nx.induced_subgraph(G,g))
        # graph_stats(g_cmp)
        for edge in g_cmp.edges():
            rescore = ave_number_of_p3(edge, G)
            new_edges.append([edge[0], edge[1], rescore])

    ### Writing new edge to gcomplement file
    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    with gcomplement.open("w") as f:
        for e in new_edges:
            f.writelines("{0} {1} {2}\n".format(e[0],e[1],e[2]))

    ###CUBCO+
    G_bar = nx.Graph()
    G_cmp = nx.read_weighted_edgelist(gcomplement, create_using = nx.Graph(), nodetype = str)
    printc('Performing CUBCO+ ......')
    for g in nx.connected_components(G):

        if len(g)>3:
            G_bar = nx.induced_subgraph(G_cmp,g).copy()
            result = min_cut()
        else:
            result = [g]
        
        i = 1
        edges = 0
        with outputfile.open("a") as f:
            if G_bar.number_of_edges()>0:
                for r in result:
                    if r:
                        sbg = nx.induced_subgraph(G,r)
                        sbg_cmp = list(nx.connected_components(sbg))
                        if len(sbg_cmp) > 1:
                            for cmp in sbg_cmp:
                                f.writelines("%s " % c for c in cmp)
                                f.write("\n")
                            edges += len(sbg.edges())
                        else: 
                            f.writelines("%s " % v for v in r)
                            f.write("\n")
                            edges += len(sbg.edges())
            else:
                sbg = nx.induced_subgraph(G,g)
                sbg_cmp = list(nx.connected_components(sbg))
                # printc(sbg_cmp)
                if len(sbg_cmp) > 1:
                    for cmp in sbg_cmp:
                        f.writelines("%s " % c for c in cmp)
                        f.write("\n")
                    edges += len(sbg.edges())
                else: 
                    f.writelines("%s " % v for v in g)
                    f.write("\n")
                    edges += len(sbg.edges())
                print('round after one result in a file: ', rounds, end="\r")
                rounds += 1
        if not nx.is_empty(G_bar):
            G_bar.clear()
        #G_bar.clear()
    printc('Number of removed edges is: %s' % (len(G.edges()) - edges))
    
    
    with outputfile.open("r") as f:
        cmplx = f.readlines()
    cmplx.sort(key=len, reverse=True)
    cmplx = [x.strip() for x in cmplx] 
    cmplx = {tuple(x.split(' ')) for x in cmplx}
    
    cmplx = list(cmplx)
    cmplx.sort(key=len, reverse=True)
    
    with outputfile.open("w") as f:
        for item in cmplx:
            for node in item:
                f.write("%s " % node)
            f.write("\n")


    printc("--- %s seconds ---" % (time.time() - start_time))
 