# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 20:27:04 2020

@author: somranian
"""

# This program must be roon from .../scl-bioinfo/
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
    printc("Enter 1 for sequential mode, Enter 2 for parallel mode: ", end = "")
    mode = int(input())
    if (osname == "win32" and mode == 2):
        printc("Please enter pool_thresh (default 100): ", end = "")
        pool_thresh = input()
        if pool_thresh == "": pool_thresh = 100
        printc("Please enter num_procs (default 8): ", end = "")
        num_procs = input()
        if num_procs == "": num_procs = 8
        return osname, mode, int(pool_thresh), int(num_procs)
    return osname, mode, 100, 8

if __name__ == '__main__':
    osname, mode, pool_thresh, num_procs = platformer()

    start_time = time.time()

    # Now using relative paths, linux path also happens to work for win32
    """As an example here the network PIPS_Corum_Graph is called!"""
    sample_path = "PC2P/Human/PIPS/PIPS_Corum_Graph.txt"
    network = nx.read_weighted_edgelist(sample_path, create_using = nx.Graph(), nodetype = str)
    G = network.copy()

    """ To run code sequentially, we need to call Find_CNP from PC2P_Sequential.
        To run code parallel in Windows and Unix, we nee to call Find_CNP from PC2P_ParallelMultiprocess
        To run code parallel in Linux and Mac, we nee to call Find_CNP from PC2P_ParallelRay """
    if (mode == 1):
        import PC2P_Sequential
        edge_cut = PC2P_Sequential.Find_CNP(G)
    else:
        if (osname == "linux"):
            import psutil
            import ray
            num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment.yml"
            runtime_env = {"conda": conda_env, "working_dir": "PC2P"}
            ray.init(num_cpus=num_cpus, runtime_env=runtime_env)
            import PC2P_ParallelRay
            edge_cut = PC2P_ParallelRay.Find_CNPs_V2(G)
        else:
            import PC2P_ParallelMultiprocess
            printc("Now running PC2P_ParallelMultiprocess! :: " + os.getcwd())
            edge_cut = PC2P_ParallelMultiprocess.Find_CNP(G, pool_thresh, num_procs)

    """ To save the result clusters in Graph format"""
    G_copy = G.copy()
    G_copy.remove_edges_from(edge_cut)
    output_dir = "PC2P/scl/"
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