# slab usage of an 8-point Hadamard block placed on 8 non-base points of a regime box:
# B = {k in {-1,0,1}^8 : sum k = 0, |<k,chi_a>| <= 4 (a != 0)}  (M8 variant with all 8 coordinates boxed)
import itertools, numpy as np
from collections import Counter
X = list(range(8))
chi = np.array([[(-1)**bin(a & x).count("1") for x in X] for a in range(8)])
allk = np.array(list(itertools.product((-1,0,1), repeat=8)))
B = allk[(allk.sum(1) == 0) & np.all(np.abs(allk @ chi[1:].T) <= 4, axis=1)]
print("|B|", len(B))
signs = np.array(list(itertools.product((-1,1), repeat=8)))
use = np.abs(B @ signs.T).max(0)            # max |<k,beta>| for each sign vector beta
nonconst = [i for i, s in enumerate(signs) if len(set(s)) > 1]
c = Counter(int(use[i]) for i in nonconst)
print("usage distribution over the 254 nonconstant sign vectors:", sorted(c.items()))
print("mean usage", np.mean([use[i] for i in nonconst]), "max", max(use[i] for i in nonconst), "min", min(use[i] for i in nonconst))
