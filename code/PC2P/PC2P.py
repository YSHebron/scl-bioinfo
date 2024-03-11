# Do PC2P.py --help for help on running this program.
# Sample run: python PC2P.py dummy_CYC.txt results -p

import os
import sys
import networkx as nx
import pandas as pd
import multiprocessing as mp
import time
from helper import printc
import eval_PC2P

def positive_int(x):
    i = int(x)
    if i < 1:
        raise ValueError('Nonpositive values are not allowed')
    return i

import argparse
parser = argparse.ArgumentParser(description='Perform PC2P on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('inputfile', type=str, help='relpath to PPI dataset')
parser.add_argument('outputdir', type=str, help='relpath to output dir for predicted complexes (CNPs)')
parser.add_argument('-i', metavar='ITERS', type=int, default=1, help='num of clustering iterations to produce (default 1)')
parser.add_argument('-p', action='store_true', help='assert for parallel mode (default sequential)')
parops = parser.add_argument_group('parallel options', description='if parallel mode is enabled with -p, the following options may be set')
parops.add_argument('--pool_thresh', nargs='?', type=positive_int, default=100, const=100, help='num of graph components to selectively trigger parallelization (for mp only, default 100)')
parops.add_argument('--num_procs', nargs='?', type=positive_int, default=8, const=8, help='num of processes created by each call to pool (for mp only, default 8)')
args = parser.parse_args()

# function to score the complex
def get_score(complex, edgesref):
    totalweight = 0
    for id1 in range(len(complex)):
        for id2 in range(1,len(complex)):
           id_a = complex[id1]
           id_b = complex[id2]
           key = f"{id_a}|{id_b}" if id_a < id_b else f"{id_b}|{id_a}"
           if key in edgesref:
               totalweight += edgesref[key]
    score = totalweight*2 / ((len(complex) * (len(complex) -1)))
    return score

def perform_cnp(args):
    G, i, inputfile, is_parallel, pool_thresh, num_procs, outputdir = args
    osname = sys.platform
    print("Iteration", i)
    if not is_parallel:
        import sequential
        printc("Hello from PC2P_Sequential.py! Current cwd :: " + os.getcwd())
        edge_cut = sequential.Find_CNP(G)
    else:
        if (osname == "linux"):
            import psutil
            import ray
            num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment.yml"
            runtime_env = {"conda": conda_env, "working_dir": "code/PC2P"}  # NOTE: This is why this should be run at the root
            ray.init(num_cpus=num_cpus, runtime_env=runtime_env)
            import parallel_ray
            edge_cut = parallel_ray.Find_CNPs_V2(G)
        else:
            import parallel_multiprocess
            printc("Now running parallel_multiprocess.py! :: " + os.getcwd())
            edge_cut = parallel_multiprocess.Find_CNP(G, pool_thresh, num_procs)

    ### Minimum Edge Cut
    G_copy = G.copy()
    G_copy.remove_edges_from(edge_cut) # G_copy now contains G's predicted CNP
    
    ### Writing to results
    outputdir = outputdir + "/"
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)

    # We save each predicted cluster (complex) in the predicted CNP of G in one line
    # Each remaining connected component (predicted BSsG coherent partitions) is treated as its own cluster
    # We write one line per predicted cluster
    # Each cluster is represented as a set in G_cnp_components
    G_cnp_components = list(nx.connected_components(G_copy))
    G_cnp_components.sort(key=len, reverse=True)
    G_cnp_components = [sorted(list(cplx)) for cplx in G_cnp_components]
    printc("First 10 complexes, sorted by decreasing number of proteins:")
    print(*G_cnp_components[:10], sep="\n")

    scorededges = {}
    neighbours = {}
    with open(inputfile) as edges_file:
        for line in edges_file:
            edge = line.split()
            id_a = edge[0]
            id_b = edge[1]
            score = edge[2]
            key = f"{id_a}|{id_b}" if id_a < id_b else f"{id_b}|{id_a}"
            scorededges[key] = float(score)
            
    with open(outputdir + 'G_PredictedComplexes_iter{}.txt'.format(i), 'w') as f:
        # complex === line
        for complex in G_cnp_components:
            # protein === node
            if len(complex) < 2: continue   # filter out single-protein complexes

            # score the complex by their weighted density
            score = get_score(complex, scorededges)
            f.write("(" + str(len(complex)) + "_" + str(score) + "): ")
            for protein in complex:
                f.write("%s " % protein)
            f.write("\n")
    
    return i

if __name__ == '__main__':
    inputfile, outputdir, iters, is_parallel, pool_thresh, num_procs = args.inputfile, args.outputdir, args.i, args.p, args.pool_thresh, args.num_procs
    
    start_time = time.time()
    
    ### Read the PPIN (given as a PPI dataset or edgelist)
    # Assumes the inputfile has scores, which may or may not be used
    G = nx.Graph()
    if inputfile.endswith(".csv"):
        # for .csv with header edge lists
        df = pd.read_csv(inputfile)
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph(), edge_attr = "score")
    elif inputfile.endswith(".tsv"):
        # for .tsv with header edge lists
        df = pd.read_csv(inputfile, sep="\t")
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph(), edge_attr = "score")
    elif inputfile.endswith(".txt"):
        # for .txt with no header edge lists
        G = nx.read_weighted_edgelist(inputfile, create_using = nx.Graph(), nodetype = str)
    printc(G)
    
    print("Size of V(G):", G.number_of_nodes())
    print("Size of E(G):", G.number_of_edges())

    ### Clustering
    # To run sequential code, we need to call Find_CNP from sequential.py
    # To run parallelized code in Windows and Unix, we need to call parallel_multiprocess.py 
    # To run parallelized code in Linux and Mac, we need to call Find_CNP from parallel_ray.py
    
    # Parallel iterations for sequential mode: Create a Pool for parallel execution
    # Sequential iterations for parallel mode
    # TODO: Justify iterations by using a stochastic approach
    #   (like by applying SWC before the parameter-free approach)
    if not is_parallel:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            # Prepare arguments for each iteration
            args_list = [(G, i, inputfile, is_parallel, pool_thresh, num_procs, outputdir) for i in range(0, iters)]
            # Use the pool to map the function to the arguments
            result = pool.map_async(perform_cnp, args_list)            
            while not result.ready():
                time.sleep(1)
            result = result.get()
            pool.close()
            pool.join()
    else:
        for i in range(0, iters):
            args_list = (G, i, inputfile, is_parallel, pool_thresh, num_procs, outputdir)
            perform_cnp(args_list)
             
    ### Evaluation (call Analysis)
    
    printc("Algorithm took %s seconds to finish." % (time.time() - start_time))