# This program must be run from scl-bioinfo root
# py code/PC2P/main.py data/intermediate/data_yeast_rand.csv code/PC2P/results 1 2         (accepts .csv PPIN)
# python code/PC2P/main.py code/PC2P/Human/PIPS/PIPS_Corum_Graph.txt code/PC2P/results 1 2     (accepts .txt PPIN)
# NOTE: Programming this must be generalizable so when using a different clustering algorithm it would be a quick transition

import os
import sys
import networkx as nx
import pandas as pd
import multiprocessing as mp
import time
from helper import printc
from typing import Tuple
import eval_PC2P

# Parameters (unrelated to the actual clustering!)
# argv[1] inputfile: relpath to PPI dataset csv file (with header p1, p2, score) or txt file (no header)
# argv[2] outputdir: relpath to output dir for predicted complexes (CNPs)
# argv[3] iters: number of clustering iterations to produce (to match run_xval.bat) 
# argv[4] mode: 1 for sequential, 2 for parallel
# argv[5] pool_thresh: only if mode 2, default 100, tells if a round will need to parallelize based on number of components
# argv[6] num_procs: only if mode 2, default 8, tells how many processes will be created for each call to pool

# Get user-given and OS-based arguments (no ValueError validation)
def get_opts() -> Tuple[str, str, str, int, int, None|int, None|int]:
    
    osname = sys.platform
    printc("Detected platform: {}".format(osname))
    
    # TODO: Use argparse for better scaling, combine with converting run_PC2P.bat to powershell script
    try: inputfile = sys.argv[1]
    except IndexError: print("No inputfile path provided. Exiting."); sys.exit()
    try: os.path.isfile(inputfile)
    except: print("Invalid inputfile. Exiting."); sys.exit()
    
    try: outputdir = sys.argv[2]
    except IndexError: print("No outputdir path provided. Exiting."); sys.exit()
    
    try: iters = int(sys.argv[3])   # Note that PC2P is deterministic, hence each iteration produce the same clustering
    except IndexError: print("No number of iterations provided. Exiting."); sys.exit()
    except ValueError: print("Invalid mode. Exiting."); sys.exit()
    
    try: mode = int(sys.argv[4])
    except IndexError: print("No mode provided. Exiting."); sys.exit()
    except ValueError: print("Invalid mode. Exiting."); sys.exit()
    
    if mode==1:
        printc("Run mode: sequential")
        return osname, inputfile, outputdir, iters, mode, None, None
    elif mode==2:
        printc("Run mode: parallel")
        # Return for parallel_multiprocess
        if osname=="win32":
            try: pool_thresh = int(sys.argv[5])
            except IndexError: pool_thresh = 100
            except ValueError: print("Invalid pool_thresh. Exiting."); sys.exit()
            try: num_procs = int(sys.argv[6])
            except IndexError: num_procs = 8
            except ValueError: print("Invalid num_procs. Exiting."); sys.exit()
            printc("pool_thresh: " + str(pool_thresh))
            printc("num_procs: " + str(num_procs))
            return osname, inputfile, outputdir, iters, mode, pool_thresh, num_procs
        # Return for parallel_ray
        return osname, inputfile, outputdir, iters, mode, None, None
    else: print("Mode argument is not 1 or 2. Exiting."); sys.exit()

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
    G, i, osname, mode, pool_thresh, num_procs, outputdir = args
    iter_time = time.time()
    print("Iteration", i)
    if (mode == 1):
        import sequential
        edge_cut = sequential.Find_CNP(G)
    else:
        if (osname == "linux"):
            import psutil
            import ray
            num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment.yml"
            runtime_env = {"conda": conda_env, "working_dir": "code/PC2P"}
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

    # # For unweighted inputfile (default)
    # nx.write_edgelist(G_copy, outputdir + "G_PredictedCNP_edgelist.txt", data=False)
    # nx.write_edgelist(G_copy, outputdir + "G_PredictedCNP_edgelist.gz", data=False)
    # # For weighted inputfile (not really used)
    # nx.write_weighted_edgelist(G_copy, outputdir + "G_PredictedCNP_weightededgelist.txt")
    # nx.write_weighted_edgelist(G_copy, outputdir + "G_PredictedCNP_weightededgelist.gz")
    # # Note: Use gzip -dk file.gz to extract .gz in linux

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
    
    printc("Iteration {} took {} seconds to finish.".format(i, time.time() - iter_time))

if __name__ == '__main__':
    osname, inputfile, outputdir, iters, mode, pool_thresh, num_procs = get_opts()
    
    start_time = time.time()
    G = nx.Graph()
    if inputfile.endswith(".csv"):
        # for .csv with header inputs
        df = pd.read_csv(inputfile)
        print(df)
        G = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph(), edge_attr = "score")
    elif inputfile.endswith(".txt"):
        # for .txt no header inputs
        G = nx.read_weighted_edgelist(inputfile, create_using = nx.Graph(), nodetype = str)
    print(G)
    
    print("Size of V(G):", G.number_of_nodes())
    print("Size of E(G):", G.number_of_edges())
    
    ### Read the score edges file
    # Assumes the inputfile has scores

    ### Clustering
    # To run sequential code, we need to call Find_CNP from sequential.py
    # To run parallelized code in Windows and Unix, we need to call parallel_multiprocess.py 
    # To run parallelized code in Linux and Mac, we need to call Find_CNP from parallel_ray.py
    
    # Parallel iterations for sequential mode: Create a Pool for parallel execution
    # Sequential iterations for parallel mode
    # TODO: Justify iterations by using a stochastic approach
    #   (like by applying SWC before the parameter-free approach)
    if mode==1:
        with mp.Pool(processes=mp.cpu_count()) as pool:
            # Prepare arguments for each iteration
            args_list = [(G, i, osname, mode, pool_thresh, num_procs, outputdir) for i in range(0, iters)]
            # Use the pool to map the function to the arguments
            pool.map(perform_cnp, args_list)
    elif mode==2:
        for i in range(0, iters):
            args_list = (G, i, osname, mode, pool_thresh, num_procs, outputdir)
            perform_cnp(args_list)
            
    ### Evaluation
    
    
    printc("Algorithm took %s seconds to finish." % (time.time() - start_time))