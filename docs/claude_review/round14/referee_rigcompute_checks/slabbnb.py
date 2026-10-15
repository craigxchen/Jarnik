# slabbnb.py -- complete interval branch-and-bound for the slab system of a joint witness at norm <= K:
#   exists phi in [0,pi)^r with |sin(<v,phi> - s pi/4)| >= a_(v,s) := nu_s m_(v,s)/sqrt(Norm Z_v) for all v in V_K, s?
# Box pruned only if some (v,s) is violated on the WHOLE box (with a 1e-9 safety margin in favour of feasibility);
# box accepted if every constraint holds on the whole box (centre then re-verified at 50 digits).
# Exhausting the tree = certificate of infeasibility (floating point with margin; labelled as such).
import sys, math, pickle
import numpy as np
import mpmath as mp
from reflib import *
from rverify import R
from slab import setup, verify_mp

def bnb(A, a, r, cap=3_000_000):
    if (a >= 1).any(): return "INFEASIBLE(a>=1)", None, 0
    w = np.arcsin(a)  # forbidden half-widths
    z0 = np.arange(4) * np.pi / 4
    absA = np.abs(A)
    stack = [(np.full(r, np.pi / 2), np.full(r, np.pi / 2))]
    nodes = 0
    while stack:
        c, h = stack.pop(); nodes += 1
        if nodes > cap: return "ABORT", None, nodes
        xc = A @ c; xh = absA @ h
        d = (xc[:, None] - z0[None, :]) % np.pi
        d = np.minimum(d, np.pi - d)  # distance to nearest forbidden-zone centre
        if (d + xh[:, None] < w - 1e-9).any():
            continue  # some constraint violated on the whole box
        if ((d - xh[:, None] >= w) | (w == 0)).all() and (xh < np.pi / 2).all():
            return "FOUND", c, nodes
        j = int(np.argmax(h)); h2 = h.copy(); h2[j] /= 2
        c1 = c.copy(); c1[j] -= h2[j]; c2 = c.copy(); c2[j] += h2[j]
        stack.append((c1, h2)); stack.append((c2, h2))
    return "INFEASIBLE", None, nodes

if __name__ == "__main__":
    fn = sys.argv[1]
    W = R(open(fn, "rb")).load(); P, K = W["P"], W["K"]; r = len(P)
    labs = {q: [int(x) for x in W["labels"][q]["lifts"][-1]] for q in W["labels"]}
    V, A, a, dem = setup(P, K, labs)
    st, phi, nodes = bnb(A, a, r)
    extra = ""
    if phi is not None: extra = f" mp-verified slack {mp.nstr(verify_mp([float(x) for x in phi], V, P, dem), 6)} phi {list(phi)}"
    print(fn.split("/")[-1], "slab system at norm <=", K, ":", st, "nodes", nodes, extra, flush=True)
