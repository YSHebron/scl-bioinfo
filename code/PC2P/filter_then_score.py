import os
"""Filter then score:
Note that this program only works for the Yeast datasets.
This has effectively demonstrated that the perpair filtering is the best filtering.
"""
# TODO: Negatome

def filterPerProtein(ppinname, ppinfile, gldstd, outputdir, species="Yeast"):
    filtering = "perprotein"
    rawPPIN = "code/PC2P/{}/{}/{}".format(species, ppinname, ppinfile)
    reffile = "code/PC2P/{}/{}_complexes.txt".format(species, gldstd)
    outputfile = outputdir + "/{}_{}_{}_weighted.txt".format(
        ppinname, gldstd, filtering
    )
    present = set()
    scorededges = {}
    with open(reffile) as f1, open(rawPPIN) as f2:
        for line in f1:
            for pid in line.split():
                present.add(pid)
        for line in f2:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)

    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    with open(outputfile, "w+") as f:
        for key in scorededges:
            if key[0] in present or key[1] in present:
                f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))


def filterPerPair(ppinname, ppinfile, gldstd, outputdir, species="Yeast"):
    filtering = "perpair"
    rawPPIN = "code/PC2P/{}/{}/{}".format(species, ppinname, ppinfile)
    reffile = "code/PC2P/{}/{}_complexes.txt".format(species, gldstd)
    outputfile = outputdir + "/{}_{}_{}_weighted.txt".format(
        ppinname, gldstd, filtering
    )
    present = set()
    scorededges = {}
    with open(reffile) as f1, open(rawPPIN) as f2:
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
    with open(outputfile, "w+") as f:
        for key in scorededges:
            if frozenset(key) in present:
                f.write("%s %s %f\n" % (key[0], key[1], scorededges[key]))


def filterDirect(ppinname, ppinfile, gldstd, outputdir, species="Yeast"):
    filtering = "direct"
    rawPPIN = "code/PC2P/{}/{}/{}".format(species, ppinname, ppinfile)
    directfile = "code/PC2P/{}/{}/{}_{}_Graph.txt".format(
        species, ppinname, ppinname, gldstd
    )
    outputfile = outputdir + "/{}_{}_{}_weighted.txt".format(
        ppinname, gldstd, filtering
    )
    scorededges = {}
    with open(rawPPIN) as f1:
        for line in f1:
            p1, p2, weight = line.split()
            key = (p1, p2) if p1 < p2 else (p2, p1)
            scorededges[key] = float(weight)
    with open(directfile) as interm, open(outputfile, "w+") as f:
        for line in interm:
            p1, p2 = line.split()
            linepair = (p1, p2) if p1 < p2 else (p2, p1)
            for key in scorededges:
                if linepair == key:
                    f.write("%s %s %f\n" % (linepair[0], linepair[1], scorededges[key]))


if __name__ == "__main__":
    # Raw, unfiltered Omranian PPINs with scored edges
    # NOTE: Because Omranian's Human PPINs are completely unweighted, we defer their inclusion.
    species = "Yeast"
    ppins = {
        "Collins": "collins2007.txt",
        "Gavin": "gavin2006_socioaffinities_rescaled.txt",
        "KroganCore": "krogan2006_core.txt",
        "KroganExt": "krogan2006_extended.txt",
    }
    gldstds = ["CYC", "SGD"]
    outputdir = "data/{}/FilteredPPINs".format(species)
    for ppinname in ppins.keys():
        for gldstd in gldstds:
            filterPerProtein(ppinname, ppins[ppinname], gldstd, outputdir)
            filterPerPair(ppinname, ppins[ppinname], gldstd, outputdir)
            filterDirect(ppinname, ppins[ppinname], gldstd, outputdir)
