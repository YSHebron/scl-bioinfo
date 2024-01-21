# Program that evaluates predicted complexes by getting the precision and recall
# python code/PC2P/eval_PC2P.py code/PC2P/results/KroganExt/G_PredictedComplexes_iter0.txt
import argparse
import os
import sys
from helper import printc
from typing import Tuple

# Calculate Jaccard similarity between P and C
# P: Preducted cluster
# C: Reference complex

def jaccard_similarity(P, C):
    intersection_size = len(P.intersection(C))
    union_size = len(P.union(C))
    
    if union_size == 0:
        return 0.0  # Handle the case where both sets are empty
    
    return intersection_size / union_size

if __name__ == '__main__':
    # Reference file for the complexes of yeast
    complexfile = "code/PC2P/Yeast/CYC2008_complexes.txt"
    
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Evaluate PC2P Arguments')

    # Add command line arguments
    parser.add_argument('arg1', type=str, help='Input file')
    #parser.add_argument('arg2', type=str, help='Output directory')

    # Parse the command line arguments
    args = parser.parse_args()

    # Parameters
    # argv[0] = inputfile: file that contains the predicted complexes to be evaluated
    # argv[1] = outputdir: relpath for the output 
    try: inputfile = args.arg1
    except IndexError: print("No inputfile path provided. Exiting."); sys.exit()
    try: os.path.isfile(inputfile)
    except: print("Invalid inputfile. Exiting."); sys.exit()
    
    # try: outputdir = args.arg2
    # except IndexError: print("No outputdir path provided. Exiting."); sys.exit()

    ### Read complexes
    with open(complexfile) as complex_file:
        complexes= [line.strip().split(" ") for line in complex_file.readlines()]

    # print("num complexes read =", len(complexes))
    # print(complexes)

    
    # Read the predicted complexes
    # clusters is the list
    
    with open(inputfile, 'r') as file:
        predicted = [line.strip().split(" ") for line in file.readlines()]
    

    # Determine if any of the predicted clusters matches with any of the reference complexes by using Jaccard similarity
    match_Complexes = 0
    for comp in predicted:
        a = set(comp)
        for ref in complexes:
            b = set(ref)
            jaccard = jaccard_similarity(a,b)
            # Predicted cluster and reference complex is a match if jaccard > 0.5
            if jaccard > 0.5:
                match_Complexes += 1
                break; # found a match

    print("Number of matches: ", match_Complexes, end="\n")

    # Calculate Precision and Recall
    precision = match_Complexes/len(predicted)
    recall = match_Complexes/len(complexes)

    print("Precision: ", precision, "\n")
    print("Recall: ", recall)

    