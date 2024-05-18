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

### Prepare for Scoring Detected Clusters (done simultaneously with writing to result)
scorededges = {}
for id_a, id_b, score in G.edges(data=True):
    # keys are tuples of protein pairs, and values are their scores
    key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
    scorededges[key] = score["weight"]
    
# ### Writing to results
# if not outputdir.is_dir():
#     outputdir.mkdir(parents=True)
# outputfile = Path(outputdir, inputfile.stem.removesuffix("_weighted") + "_predicted.txt")
# with outputfile.open("w") as f:
#     # complex === line
#     for complex in G_cnp_components:
#         # protein === node
#         # Score the complex by their weighted density
#         # Each line: (len(complex)_score): p1 p2 p3 ...
#         score = get_score(complex, scorededges)
#         f.write(f"({len(complex)}_{score}): ")
#         for protein in complex:
#             f.write("%s " % protein)
#         f.write("\n")