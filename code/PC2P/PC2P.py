# Do PC2P.py --help for help on running this program.
# Sample run: python code/PC2P/PC2P.py data/Yeast/FilteredPPINs/Collins_CYC_direct_weighted.txt data/Results/Dummy -p

import ray
import networkx as nx
import time
import argparse
from pathlib import Path
from utils import read_ppin_to_graph, printc, positive_int, graph_stats, graph_memory

parser = argparse.ArgumentParser(description='Perform PC2P on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('ppinfile', type=Path, help='path to PPI dataset to be processed')
parser.add_argument('outfile', type=Path, help='writepath for PC2P predicted clusters')
parser.add_argument('-p', type=str, choices=["ray", "mp"], nargs='?', const='mp',
                    help='assert for parallel mode and select parallelism package (const selected is "mp")')
parops = parser.add_argument_group('parallel options', description='if parallel mode is enabled with -p mp, the following options may be set')
parops.add_argument('--pool_thresh', nargs='?', type=positive_int, default=100, const=100, help='num of graph components to selectively trigger parallelization (for mp only, default 100)')
parops.add_argument('--num_procs', nargs='?', type=positive_int, default=16, const=16, help='num of processes created by each call to pool (for mp only, default 16)')
args = parser.parse_args()

# Generate minimumm node cut (edge_cut) and apply it to G
def perform_cnp(G):
    parallel = args.p
    edge_cut = []   # store edge_cut that emulates effect of min node cut (why not just use node_cut directly)
    if parallel is None:
        import sequential
        printc("Now running PC2P_Sequential.py!")
        edge_cut = sequential.Find_CNP(G)
    else:
        printc("Running in parallel mode: %s" % parallel, "red")
        if parallel == "ray":
            # num_cpus = psutil.cpu_count(logical=False)
            conda_env = "environment_ray.yml"
            runtime_env = {"conda": conda_env, "working_dir": "code/PC2P"}
            ray.init(runtime_env=runtime_env)
            print("Ray resources: ", ray.available_resources())
            import parallel_ray
            printc("Now running parallel_ray.py!")
            edge_cut = parallel_ray.Find_CNP(G)
        elif parallel == "mp":
            import parallel_multiprocess
            printc("Now running parallel_multiprocess.py!")
            edge_cut = parallel_multiprocess.Find_CNP(G, args.pool_thresh, args.num_procs)
    
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
    print(*G_cnp_components[:10], sep='\n')
    
    ### Writing to results
    outfile = args.outfile
    outfile.parent.mkdir(exist_ok=True, parents=True)
    with open(outfile, 'w') as f:
        # complex === line
        for complex in G_cnp_components:
            # protein === node
            # Score the complex by their weighted density
            # Each line: p1 p2 p3 ...
            for protein in complex:
                f.write(f"{protein} ")
            f.write('\n')
        
    ### Return G_cnp_components for analysis phase
    return G_cnp_components


if __name__ == '__main__':    
    start_time = time.time()
    
    G = read_ppin_to_graph(args.ppinfile)
    
    graph_stats(G)
    graph_memory(G)
    print(nx.to_numpy_array(G))

    ### Clustering (Coherent Network Partitioning)
    # To run sequential code, we call Find_CNP from sequential.py
    # To run parallelized code, we call either parallel_multiprocess.py or parallel_ray.py 
    perform_cnp(G)
            
    printc("Algorithm took %d seconds to finish." % (time.time() - start_time))