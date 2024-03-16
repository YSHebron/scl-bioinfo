PPINNAME = "Collins"
PPINFILE = "collins2007.txt"
GLDSTDNAME = "CYC2008"
OUTPUTDIR = "code/PC2P/NewDatasets"

def filterPerProtein():
    pass

def filterPerProteinPair():
    pass

def directFilter():
    pass

# rawPPI = "code/PC2P/Yeast/%s/%s" % (PPINNAME, PPINFILE)
# goldStandard = "code/PC2P/Yeast/%s_complexes.txt" % GLDSTDNAME

# present = set()
# scorededges = {}
# with open(goldStandard) as gldstd, open(rawPPI) as raw:
#     for line in gldstd:
#         complex = line.split()
#         for i in range(0, len(complex)):
#             for j in range(0, len(complex)):
#                 if i == j: continue
#                 p1, p2 = complex[i], complex[j]
#                 present.add((p1, p2) if p1 < p2 else (p2, p1))
#     print(present)
    
#     for line in raw:
#         p1, p2, weight = line.split()
#         key = (p1, p2) if p1 < p2 else (p2, p1)
#         scorededges[key] = float(weight)
                
#     print(scorededges)
    
# import os
# if not os.path.isdir(OUTPUTDIR):
#     os.mkdir(OUTPUTDIR)
# with open(OUTPUTDIR + '/%s_%s_Weighted.txt' % (PPINNAME, GLDSTDNAME), 'w') as f:
#     for key in scorededges:
#         if key in present:
#             f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))   


def directfilter(rawPPI, intermediate):
    scorededges = {}
    with open(rawPPI) as raw:
        for line in raw:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)
    
    lineno = 0
    with open(intermediate) as interm, open("code/PC2P/Collins_CYC_Weighted_Direct.txt", "w+") as f:
        for line in interm:
            p1, p2 = line.split()
            for pair in scorededges:
                if pair == (p1, p2) or pair == (p2, p1):
                    f.write("%s %s %f\n" % (p1, p2, scorededges[pair]))
                    lineno += 1
    
    print("Number of pairs", lineno)
                    
directfilter("code/PC2P/Yeast/Collins/collins2007.txt", "code/PC2P/Yeast/Collins/Collins_CYC_Graph.txt")