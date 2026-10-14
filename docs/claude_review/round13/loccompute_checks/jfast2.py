"""J(P) new generators, selection mod p, exact lift (ratrec over a 62-bit prime) + exact verification.
Every returned form is an exact element of J(P)_i (verified: its full kernel vector vanishes on P over Z).
Selected forms are independent mod PR modulo y*J(P)_{i-1}-reductions (computed mod PR) -> also over Q."""
import math, numpy as np
from fractions import Fraction
from flint import nmod_mat
from loclib import bmonos, PR
from echelon import Echelon
from jfast import P62, ratrec

def _rank(rows, ncols, p):
    if not rows: return 0
    return nmod_mat(len(rows), ncols, [int(x) % p for r in rows for x in r], p).rank()

def jproj_gens_fast2(P, b, maxdeg, mind=1):
    P = np.array([tuple(int(z) for z in p) for p in P], dtype=np.int64)
    out = {}; rho = None; fails = 0
    prev_span_p = []     # mod-PR spanning set of J(P)_{i-1} (top parts)
    for i in range(1, maxdeg + 1):
        top = bmonos(b, i); low = [e for t in range(i) for e in bmonos(b, t)]
        mons = low + top
        E = np.ones((len(P), len(mons)), dtype=np.int64)
        for j, e in enumerate(mons):
            for z in range(b):
                if e[z]: E[:, j] *= P[:, z] ** e[z]
        Em = nmod_mat(E.shape[0], E.shape[1], [int(x) % P62 for x in E.ravel()], P62)
        X, nul = Em.nullspace()
        topidx = {e: a for a, e in enumerate(top)}
        prevm = bmonos(b, i - 1)
        shifted = []
        for v in prev_span_p:
            for q in range(b):
                w = [0] * len(top)
                for a, e in enumerate(prevm):
                    if v[a]:
                        e2 = list(e); e2[q] += 1; w[topidx[tuple(e2)]] = (w[topidx[tuple(e2)]] + v[a]) % PR
                shifted.append(w)
        # basis of shifted mod PR
        ech = Echelon(len(top)); cur = []
        for w in shifted:
            if ech.add(w): cur.append(w)
            if ech.rank == len(top): break
        r = ech.rank
        new = []
        cand = [[int(X[k, c]) for k in range(len(mons))] for c in range(nul)]
        for vec in cand:
            if r == len(top): break
            topv = vec[len(low):]
            if not any(topv): continue
            ft = [ratrec(x, P62) for x in topv]
            if any(f is None for f in ft): fails += 1; continue
            dt = 1
            for f in ft: dt = dt * f.denominator // math.gcd(dt, f.denominator)
            if not any(ech.reduce([int(f * dt) % PR for f in ft])): continue
            fr = [ratrec(x, P62) for x in vec]
            if any(f is None for f in fr): fails += 1; continue
            den = 1
            for f in fr: den = den * f.denominator // math.gcd(den, f.denominator)
            iv = [int(f * den) for f in fr]
            g = 0
            for z in iv: g = math.gcd(g, z)
            iv = [z // g for z in iv]
            tv = iv[len(low):]
            tvp = [z % PR for z in tv]
            if not any(ech.reduce(tvp)): continue
            # exact verification E @ iv == 0
            mx = max(abs(z) for z in iv)
            if mx * int(np.abs(E).max()) * len(mons) < 2 ** 62:
                ok = not (E @ np.array(iv, dtype=np.int64)).any()
            else:
                Eo = E.astype(object); ok = not any(Eo.dot(np.array(iv, dtype=object)))
            if not ok: fails += 1; continue
            ech.add(tvp); cur.append(tvp); new.append(tv); r += 1
        out[i] = [{e: v[a] for a, e in enumerate(top) if v[a]} for v in new]
        prev_span_p = cur
        if r == len(top) and rho is None: rho = i - 1
    return out, rho, fails
