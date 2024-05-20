# To emulate Yong and Wong for predicts, we also add number of correct matches
class Cluster:
    def __init__(self, proteins = set(), score = 0, id = None):
        self.proteins = set(proteins)
        self.score = score
        self.id = id
    
    def __str__(self):
        return "(%d_%.6f): %s" % (len(self.proteins), self.score, " ".join(self.proteins))