import os
import argparse
from pathlib import Path
"""Filter then rescore:
This program only works for the Yeast datasets.
Demonstrated that the perpair filtering is the best filtering.
"""

# TODO: Negatome


parser = argparse.ArgumentParser(description='Preprocess PPIN before clustering. Currently only does per protein pair filtering.')
parser.add_argument('ppin_file', type=Path, help='relpath to raw PPIN')
parser.add_argument('ppin_name', type=str, help='name of the raw PPIN, e.g. BioGRID, IntAct, Collins, Gavin, etc.')
parser.add_argument('reference_file', type=Path, help='relpath to gold standard to assist in filtering')
parser.add_argument('reference_name', type=str, help='name of the gold standard, e.g. CYC, SGD, Corum, etc.')
parser.add_argument('output_dir', type=Path, help='directory to output preprocessed PPIN and other preprocessing files')
args = parser.parse_args()


def filterPerPair(ppin: Path, ppin_name, gldstd: Path, gldstd_name, outputdir: Path):
    filtered_file = Path(outputdir, f"{ppin_name}_{goldstd_name}_perpair_weighted.txt")
    present = set()
    scorededges = {}
    with open(gldstd) as f1, open(ppin) as f2:
        for line in f1:
            for p1 in line.split():
                for p2 in line.split():
                    if p1 == p2:
                        continue
                    present.add(frozenset([p1, p2]))
        for line in f2:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    with open(filtered_file, "w+") as f:
        for key in scorededges:
            if frozenset(key) in present:
                f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))


if __name__ == "__main__":
    # Entry point for preprocessing. Also the Python entry point for the entire pipeline. 
    ppin, ppin_name, goldstd, goldstd_name, output_dir = args
    filterPerPair(Path(ppin), ppin_name, Path(goldstd), goldstd_name, Path(output_dir))
            
    # Add other preprocessing directives here
