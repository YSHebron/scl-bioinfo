ppin = "code/PC2P/Yeast/data_yeast.txt"
ppin_new = "code/PC2P/Yeast/Integrated_CYC_ppirel.txt"      # this is the necessary filename format
complex = "code/PC2P/Yeast/complexes_yeast.txt"
complex_new = "code/PC2P/Yeast/CYC_complexes_integ.txt"

def write_new_ppin(ppin, ppin_new):
    with open(ppin, "r") as f1, open(ppin_new, "w+") as f2:
        for line in f1:
            temp = line.split()
            if temp[2] == "PPIREL":
                f2.write("{} {} {}\n".format(temp[0], temp[1], temp[3]))

def write_new_complexfile(complexfile, complexfile_new):
    complexes = dict()
    
    ### Read and generate complexes
    with open(complexfile) as f:
        for line in f:
            temp = (line.split())[:2]
            pid, cid = temp[0], temp[1]
            if cid not in complexes.keys():
                complexes[cid] = [pid]
            else:
                complexes[cid].append(pid)

    with open(complexfile_new, "w+") as f:
        for cid in complexes:
            f.write("{}\n".format(' '.join(complexes[cid])))

write_new_ppin(ppin, ppin_new)
write_new_complexfile(complex, complex_new)