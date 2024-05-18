# Do PC2P.py --help for help on running this program.
# Sample run: python code/PC2P/PC2P.py data/Yeast/FilteredPPINs/Collins_CYC_direct_weighted.txt data/Results/Dummy -p
# TODO: Hub proteins

import ray
import networkx as nx
import pandas as pd
import time
import argparse
from pathlib import Path
from helper import printc, positive_int, graph_stats, graph_memory

parser = argparse.ArgumentParser(description='Perform PC2P on PPI dataset and evaluate results. Can work with csv (with header p1, p2, score) or txt (no header) inputs. For scored clusters, use weighted edge list.')
parser.add_argument('inputfile', type=Path, help='relpath to PPI dataset to be processed')
parser.add_argument('outputdir', type=Path, help='relpath to output dir for predicted complexes')
parser.add_argument('-p', type=str, choices=["ray", "mp"], nargs='?', const='mp',
                    help='assert for parallel mode and select parallelism package (const selected is "mp")')
parops = parser.add_argument_group('parallel options', description='if parallel mode is enabled with -p mp, the following options may be set')
parops.add_argument('--pool_thresh', nargs='?', type=positive_int, default=100, const=100, help='num of graph components to selectively trigger parallelization (for mp only, default 100)')
parops.add_argument('--num_procs', nargs='?', type=positive_int, default=16, const=16, help='num of processes created by each call to pool (for mp only, default 16)')
args = parser.parse_args()

# Generate minimumm node cut (edge_cut) and apply it to G, scores components (complexes) of G
# and evaluates G
def perform_cnp(G, inputfile, outputdir, parallel, pool_thresh, num_procs):
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
    
    ### Writing to results
    if not outputdir.is_dir():
        outputdir.mkdir(parents=True)
    outputfile = Path(outputdir, inputfile.stem.removesuffix("_weighted") + "_predicted.txt")
    with outputfile.open("w") as f:
        # complex === line
        for complex in G_cnp_components:
            # protein === node
            # Score the complex by their weighted density
            # Each line: p1 p2 p3 ...
            for protein in complex:
                f.write("%s " % protein)
            f.write("\n")
        
    ### Return G_cnp_components for analysis phase
    return G_cnp_components


if __name__ == '__main__':
    inputfile, outputdir, parallel, pool_thresh, num_procs = args.inputfile, args.outputdir, args.p, args.pool_thresh, args.num_procs
    
    start_time = time.time()
    
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
    
    printc("CWD:\t%s " % Path.cwd())
    printc("Input:\t%s" % inputfile)
    printc("Output:\t%s" % Path(outputdir, inputfile.stem.removesuffix("_weighted") + "_predicted.txt"))
    graph_stats(G)
    graph_memory(G)
    print(nx.to_numpy_array(G))

    ### Clustering (Coherent Network Partitioning)
    # To run sequential code, we call Find_CNP from sequential.py
    # To run parallelized code, we call either parallel_multiprocess.py or parallel_ray.py 
    
    perform_cnp(G, inputfile, outputdir, parallel, pool_thresh, num_procs)
            
    printc("Algorithm took %d seconds to finish." % (time.time() - start_time))
    
    # TODO: Add other clustering directives