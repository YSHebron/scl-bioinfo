# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:27:04 2020

@author: somranian
"""

# This program must be run from .../scl-bioinfo/
import os
import sys
import networkx as nx
import numpy as np
import pandas as pd
import time
from helper import printc
from multiprocessing import Value

def platformer():
    osname = sys.platform
    printc("Detected platform: {}".format(osname))
    # printc("Enter 1 for sequential mode, Enter 2 for parallel mode: ", end = "")
    mode = int(sys.argv[1])
    if (osname == "win32" and mode == 2):
        try:
            pool_thresh = sys.argv[2]
        except:
            pool_thresh = 100
        printc("pool_thresh (default 100): " + str(pool_thresh))
        try:
            num_procs = sys.argv[3]
        except:
            num_procs = 8
        printc("num_procs (default 8): " + str(num_procs))
        return osname, mode, int(pool_thresh), int(num_procs)
    return osname, mode, 100, 8

if __name__ == '__main__':
    osname, mode, pool_thresh, num_procs = platformer()

    start_time = time.time()

    # Now using relative paths, linux path also happens to work for win32
    """As an example here the network PIPS_Corum_Graph is called!"""
    sample_path = "data/intermediate/data_yeast_rand.csv"
    df = pd.read_csv(sample_path)
    print(df)
    network = nx.from_pandas_edgelist(df, source = "p1", target = "p2", create_using = nx.Graph(), edge_attr = None)
    # network = nx.read_weighted_edgelist(sample_path, create_using = nx.Graph(), nodetype = str)
    G = network.copy()
    # nx.draw_networkx(G)
    print(G)
    
    print("V:", G.number_of_nodes())
    print("E:", G.number_of_edges())

    """ To run sequential code, we need to call Find_CNP from sequential.py
        To run parallelized code in Windows and Unix, we need to call parallel_multiprocess.py 
        To run parallelized code in Linux and Mac, we need to call Find_CNP from parallel_ray.py """
    if (mode == 1):
        import sequential
        edge_cut = sequential.Find_CNP(G)
    else:
        if (osname == "linux"):
            import psutil
            import ray
            num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment.yml"
            runtime_env = {"conda": conda_env, "working_dir": "PC2P"}
            ray.init(num_cpus=num_cpus, runtime_env=runtime_env)
            import parallel_ray
            edge_cut = parallel_ray.Find_CNPs_V2(G)
        else:
            import parallel_multiprocess
            printc("Now running parallel_multiprocess.py! :: " + os.getcwd())
            edge_cut = parallel_multiprocess.Find_CNP(G, pool_thresh, num_procs)

    """ To save the result clusters in Graph format"""
    G_copy = G.copy()
    G_copy.remove_edges_from(edge_cut)
    output_dir = "code/PC2P/results/"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    #--- If the edges are unweighted ----------
    nx.write_edgelist(G_copy, output_dir + "STRING_CNPPredicted_V4.edgelist.gz", data=False)
    #--- If the edges are weighted ----------
    nx.write_weighted_edgelist(G_copy, output_dir + 'STRING_CNPPredicted_V4.weighted.edgelist.gz')
    # Note: Use gzip -dk file.gz to extract .gz

    """ We save each predicted cluster in one line """
    # I.e. one line per predicted cluster
    G_cnp_components = list(nx.connected_components(G_copy))
    G_cnp_components.sort(key=len, reverse=True)

    with open(output_dir + 'G_PredictedClusters.txt', 'w') as f:
        for item in G_cnp_components:
            for node in item:
                f.write("%s " % node)
            f.write("\n")

    printc("Algorithm took %s seconds to finish." % (time.time() - start_time))