rawPPI = "code/PC2P/Yeast/KroganExt/krogan2006_extended.txt"

goldStandard = "code/PC2P/Yeast/SGD_complexes.txt"

present = set()
scorededges = {}
with open(rawPPI) as f1, open(goldStandard) as f2:
    present = set()
    for line in f2:
        if line != "\n":
            for pid in line.split():
                present.add(pid)
    print(present)
    
    for line in f1:
        if line != "\n":
            print(line.split())
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)
                
    print(scorededges)
    
import os
outputdir = "code/PC2P/Yeast/KroganExt"
if not os.path.isdir(outputdir):
    os.mkdir(outputdir)
with open(outputdir + '/KroganExt_SGD_weighted.txt', 'w') as f:
    for key in scorededges:
        if key[0] in present or key[1] in present:
            f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))   
