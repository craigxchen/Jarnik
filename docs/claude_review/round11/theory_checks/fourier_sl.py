"""Boolean SL on 'Fourier-uniform' sets (evidence; symmetric model with all 2^k - 1 Hadamard slabs).

Coordinates = points of F_2^k (M = 2^k).  u_q(y) = (-1)^{q.y}, q != 0 (all nonzero characters, r = M - 1).
X = {Y subset F_2^k : |sum_{y in Y} (-1)^{q.y}| <= s for all q != 0}; the empty set is in X.
SL_bool(m): exists phi on X, phi(empty) != 0, sum_{Y in X, Y contains J} phi(Y) = 0 for all |J| < m.
X and the constraint family are invariant under AGL(k,2); averaging reduces to invariant phi:
unknowns = AGL-orbits of X, equations = AGL-orbits of J (|J| < m).  Exact rational arithmetic (sympy).
Also: product certificates at the empty set = families of >= m disjoint nonempty blocks all of whose unions lie in X.
"""
import itertools, sys
from fractions import Fraction
import numpy as np
import sympy

def agl_generators(k):
    M = 2 ** k
    pts = list(itertools.product((0, 1), repeat=k))
    idx = {p: i for i, p in enumerate(pts)}
    gens = []
    # translations by unit vectors
    for j in range(k):
        e = tuple(1 if t == j else 0 for t in range(k))
        gens.append([idx[tuple((a + b) % 2 for a, b in zip(p, e))] for p in pts])
    # GL generators: transvection x_0 += x_1, and cyclic shift of coordinates
    gens.append([idx[tuple(((p[0] + p[1]) % 2,) + p[1:])] for p in pts])
    gens.append([idx[tuple(p[1:] + p[:1])] for p in pts])
    return pts, gens

def apply(perm, mask):
    out = 0
    i = 0
    while mask:
        if mask & 1:
            out |= 1 << perm[i]
        mask >>= 1; i += 1
    return out

def orbits(masks, gens):
    masks = set(masks)
    seen = {}
    reps = []
    for mk in sorted(masks):
        if mk in seen:
            continue
        orb = [mk]; seen[mk] = len(reps)
        stack = [mk]
        while stack:
            x = stack.pop()
            for g in gens:
                y = apply(g, x)
                if y not in seen:
                    seen[y] = len(reps); orb.append(y); stack.append(y)
        reps.append(orb)
    return reps, seen

def main(k, s, mlist):
    M = 2 ** k
    pts, gens = agl_generators(k)
    chars = np.array([[(-1) ** (sum(a * b for a, b in zip(q, y)) % 2) for y in pts] for q in pts[1:]])
    X = []
    for mask in range(2 ** M):
        v = np.array([(mask >> i) & 1 for i in range(M)])
        if np.all(np.abs(chars @ v) <= s):
            X.append(mask)
    Xorb, Xseen = orbits(X, gens)
    print(f"k={k} M={M} s={s}: |X|={len(X)}  X-orbits={len(Xorb)}"); sys.stdout.flush()
    empty_orb = Xseen[0]
    for m in mlist:
        Js = [sum(1 << i for i in J) for j in range(m) for J in itertools.combinations(range(M), j)]
        Jorb, _ = orbits(Js, gens)
        Jreps = [o[0] for o in Jorb]
        N = sympy.zeros(len(Jreps), len(Xorb))
        for a, J in enumerate(Jreps):
            for b, orb in enumerate(Xorb):
                N[a, b] = sum(1 for Y in orb if (Y & J) == J)
        ker = N.nullspace()
        sl = any(vec[empty_orb] != 0 for vec in ker)
        print(f"  m={m}: J-orbits={len(Jreps)} dim ker (invariant V)={len(ker)}  SL(char0) at empty set: {sl}")
        sys.stdout.flush()

if __name__ == '__main__':
    k = int(sys.argv[1]); s = int(sys.argv[2])
    mlist = [int(a) for a in sys.argv[3:]]
    main(k, s, mlist)
