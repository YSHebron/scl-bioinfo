PPINNAME = "Collins"
PPINFILE = "collins2007.txt"
GLDSTDNAME = "CYC2008"
OUTPUTDIR = "code/PC2P/NewDatasets"

rawPPI = "code/PC2P/Yeast/%s/%s" % (PPINNAME, PPINFILE)
goldStandard = "code/PC2P/Yeast/%s_complexes.txt" % GLDSTDNAME

present = set()
scorededges = {}
with open(goldStandard) as gldstd, open(rawPPI) as raw:
    for line in gldstd:
        complex = line.split()
        for i in range(0, len(complex)):
            for j in range(0, len(complex)):
                if i == j: continue
                p1, p2 = complex[i], complex[j]
                present.add((p1, p2) if p1 < p2 else (p2, p1))
    print(present)
    
    for line in raw:
        p1, p2, weight = line.split()
        key = (p1, p2) if p1 < p2 else (p2, p1)
        scorededges[key] = float(weight)
                
    print(scorededges)
    
import os
if not os.path.isdir(OUTPUTDIR):
    os.mkdir(OUTPUTDIR)
with open(OUTPUTDIR + '/%s_%s_Weighted.txt' % (PPINNAME, GLDSTDNAME), 'w') as f:
    for key in scorededges:
        if key in present:
            f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))   
