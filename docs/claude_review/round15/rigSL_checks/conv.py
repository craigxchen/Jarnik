# conv.py -- convert a base-pipeline witness (labels with lifts) to the slverify format (adds g1, top, A, F)
import sys, pickle
from sl import *
def convert(S, labels, deep, phi):
    out = {}
    for q, lab in labels.items():
        A = len(lab["lifts"]) + (deep.get(q) or 1) - 1
        out[q] = dict(g1=tuple(Torus(q).g), top=tuple(lab["lifts"][-1]), A=A)
    return dict(P=S.P, K=S.K, F=sorted(S.F), phi=phi, labels=out)
if __name__ == '__main__':
    W = pickle.load(open(sys.argv[1], "rb")); name = sys.argv[2]
    S = setup(name, W["K"])
    Wc = convert(S, W["labels"], W["deep"], W["phi"]); Wc["Q"] = W["Q"]
    pickle.dump(Wc, open(sys.argv[3], "wb"))
