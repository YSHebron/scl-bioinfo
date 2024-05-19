import argparse
import utils
from utils import printc
from pathlib import Path

parser = argparse.ArgumentParser(description='PerProteinPair and Negatome filtering.')
parser.add_argument('ppinfile', type=Path, help='path to PPIN')
parser.add_argument('reffile', type=Path, help='path to gold standard')
parser.add_argument('outputdir', type=Path, help='directory to output preprocessed PPIN')
parser.add_argument('--negfile', type=Path, help='path to negatome', nargs='?', required=False, const=None)
parser.add_argument('--confidence', type=float, metavar='[0.0-1.0]', default=0.0, help='interaction score threshold')
args = parser.parse_args()


def filter_per_pair(ppin: dict, reffile: Path) -> dict:
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
            

def filter_per_protein(ppin: dict, reffile: Path) -> dict:
    # NOTE: May improve detection of all complexes, but not as much as per_pair.
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


def filter_with_negatome(ppin: dict, negppin: dict) -> dict:
    filtered_ppin = {}
    for ppi in ppin.keys() - negppin.keys():
        filtered_ppin[ppi] = ppin[ppi]
    
    return filtered_ppin

if __name__ == '__main__':
    ppin = utils.read_ppin_to_dict(args.ppinfile, weighted=True)
    printc("Filtering applied\tsize(E)")
    print(f"None\t\t\t{len(ppin)}")
    
    # Filtering Directives
    if args.negfile is not None:
        ppin = filter_with_negatome(ppin, utils.read_ppin_to_dict(args.negfile))
        print(f"Negatome 2.0\t\t{len(ppin)}")
    ppin = filter_per_pair(ppin, args.reffile)
    print(f"Per Protein Pair\t{len(ppin)}")
    # ppin = filter_per_protein(ppin, args.reffile)
    # print(f"Per Protein\t{len(ppin)}")
    
    outfile = args.outputdir / "ppin_filtered.txt"
    utils.write_ppin_dict_to_txt(ppin, outfile, weighted=True)
    
