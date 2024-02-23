from termcolor import colored

def printc(str, *args, **kwargs):
    print(colored(str, "red"), end = kwargs.get('end', None))
    
def write_new_complexfile(complexfile, complexfile_new):
    complexes = dict()
    
    ### Read complexes
    with open(complexfile) as f:
        for line in f:
            temp = (line.split())[:2]
            pid, cid = temp[0], temp[1]
            if cid not in complexes.keys():
                complexes[cid] = [pid]
            else:
                complexes[cid].append(pid)
    
    print(complexes)
    
    with open(complexfile_new, "w") as f:
        for cid in complexes:
            f.write("{0} {1}\n".format(cid, ' '.join(complexes[cid])))