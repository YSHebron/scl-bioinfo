# Do PC2P.py --help for help on running this program.
# Sample run: python PC2P.py dummy_CYC.txt results -p

import os
import sys
import ray
import networkx as nx
import pandas as pd
import multiprocessing as mp
import time
from helper import printc, positive_int, graph_stats, graph_memory
from typing import Tuple    # For explicit typesetting and hints

import argparse
parser = argparse.ArgumentParser(description='Perform PC2P on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('inputfile', type=str, help='relpath to PPI dataset')
parser.add_argument('outputdir', type=str, help='relpath to output dir for predicted complexes (CNPs)')
parser.add_argument('-i', metavar='ITERS', type=int, default=1, help='num of clustering iterations to produce (default 1)')
parser.add_argument('-p', action='store_true', help='assert for parallel mode (default sequential)')
parser.add_argument('-f', action='store_true', help='assert to force mp')
parops = parser.add_argument_group('parallel options', description='if parallel mode is enabled with -p, the following options may be set')
parops.add_argument('--pool_thresh', nargs='?', type=positive_int, default=100, const=100, help='num of graph components to selectively trigger parallelization (for mp only, default 100)')
parops.add_argument('--num_procs', nargs='?', type=positive_int, default=16, const=8, help='num of processes created by each call to pool (for mp only, default 8)')
args = parser.parse_args()

# function to score the complex
def get_score(complex, edgesref):
    totalweight = 0
    for id1 in range(len(complex)):
        for id2 in range(1,len(complex)):
           id_a = complex[id1]
           id_b = complex[id2]
           key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
           if key in edgesref:
               totalweight += edgesref[key]
    score = totalweight*2 / ((len(complex) * (len(complex)-1)))
    return score

# Generate minimumm node cut (edge_cut) and apply it to G, scores components (complexes) of G
# and evaluates G
def perform_cnp(args: Tuple[nx.Graph, int, str, str, bool, int, int]):
    start_time = time.time()
    G, i, inputfile, outputdir, is_parallel, pool_thresh, num_procs = args
    osname = sys.platform
    printc("Started process (iteration) %d..." % i)
    edge_cut = []   # store edge_cut that emulates effect of min node cut (why not just use node_cut directly)
    if not is_parallel:
        import sequential
        printc("Now running PC2P_Sequential.py! Current cwd :: " + os.getcwd())
        edge_cut = sequential.Find_CNP(G)
    else:
        if osname == "linux":
            # num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment_ray.yml"
            runtime_env = {"conda": conda_env}
            ray.init(runtime_env=runtime_env)
            import parallel_ray
            printc("Now running parallel_ray.py! Current cwd :: " + os.getcwd())
            edge_cut = parallel_ray.Find_CNP(G)
        else:
            import parallel_multiprocess
            printc("Now running parallel_multiprocess.py! Current cwd :: " + os.getcwd())
            edge_cut = parallel_multiprocess.Find_CNP(G, pool_thresh, num_procs)
    
    ### Apply minimum node cut (use pyplot)
    G_cnp = G.copy()
    G_cnp.remove_edges_from(edge_cut) # G_CNP now contains G's predicted CNP

    # We save each predicted cluster (complex) in the predicted CNP of G in one line
    # Each remaining connected component (predicted BSsG coherent partitions) is treated as its own cluster
    # We write one line per predicted cluster
    # Each cluster is represented as a set in G_cnp_components
    # TODO: Sort output by score, not by size
    G_cnp_components = list(nx.connected_components(G_cnp))
    G_cnp_components.sort(key=len, reverse=True)
    G_cnp_components = [sorted(list(cplx)) for cplx in G_cnp_components]
    printc("First 10 complexes, sorted by decreasing number of proteins:")
    print(*G_cnp_components[:10], sep="\n")

    ### Prepare for Scoring Detected Clusters (done simultaneously with writing to result)
    scorededges = {}
    for id_a, id_b, score in G.edges(data=True):
        key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
        scorededges[key] = score["weight"]
    
    ### Writing to results
    outputdir = outputdir + "/"
    if not os.path.exists(outputdir):
        os.mkdirs(outputdir)
    if "/" in inputfile:
        temp = inputfile.split("/")[-1].split("_")
    else:
        temp = inputfile.split("\\")[-1].split("_")
    if len(temp) > 3:
        filename = f"{temp[0]}_{temp[1]}_{temp[2]}_Predicted"
    else:
        filename = f"{temp[0]}_{temp[1]}_Predicted"
    with open(outputdir + '{}_iter{}.txt'.format(filename, i), 'w') as f:
        # complex === line
        for complex in G_cnp_components:
            # protein === node
            if len(complex) < 2:
                continue   # filter out single-protein complexes
            # Score the complex by their weighted density
            # Each line: (len(complex)_score): p1 p2 p3 ...
            score = get_score(complex, scorededges)
            f.write("(" + str(len(complex)) + "_" + str(score) + "): ")
            for protein in complex:
                f.write("%s " % protein)
            f.write("\n")
    
    printc("Iteration %d took %d seconds to finish." % (i, time.time() - start_time))
    
    ### Return G_cnp_components for analysis phase
    return G_cnp_components


if __name__ == '__main__':
    inputfile, outputdir, iters, is_parallel, pool_thresh, num_procs = str(args.inputfile), str(args.outputdir), int(args.i), bool(args.p), int(args.pool_thresh), int(args.num_procs)
    
    start_time = time.time()
    
    ### Read the PPIN (given as a PPI dataset or edgelist)
    # Assumes the inputfile has scores, which may or may not be used
    G = nx.Graph()
    if inputfile.endswith(".csv"):
        # for .csv with header p1, p2, and weight
        df = pd.read_csv(inputfile)
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.endswith(".tsv"):
        # for .tsv with header p1, p2, and weight
        df = pd.read_csv(inputfile, sep="\t")
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph, edge_attr = "weight")
    elif inputfile.endswith(".txt"):
        # for .txt with no header edge lists
        G = nx.read_weighted_edgelist(inputfile, create_using = nx.Graph, nodetype = str)
    
    printc("Processing inputfile %s" % inputfile)
    graph_stats(G)
    graph_memory(G)
    print(nx.to_numpy_array(G))

    ### Clustering (Coherent Network Partitioning)
    # To run sequential code, we call Find_CNP from sequential.py
    # To run parallelized code, we call either parallel_multiprocess.py or parallel_ray.py 
    
    # NOTE: Found an MCL implementation in Python, might be useful
    # TODO: Justify iterations by using a stochastic approach
    #   (like by applying SWC before the parameter-free approach)
    if not is_parallel:
        with mp.Pool(processes=mp.cpu_count() if mp.cpu_count() <= iters else iters) as pool:
            args_list = [(G, i, inputfile, outputdir, is_parallel, pool_thresh, num_procs) for i in range(0, iters)]
            pool.map_async(perform_cnp, args_list)  # removed spin lock to improve performance, terminate instead through taskmgr
            pool.close()
            pool.join()
    else:
        for i in range(0, iters):
            perform_cnp([G, i, inputfile, outputdir, is_parallel, pool_thresh, num_procs])
            
    printc("Algorithm took %d seconds to finish." % (time.time() - start_time))