"""NUMERICAL: the general 8-dimensional family of balanced (class beta_6) rational curves in Mbar_{0,6},
anchors+movers model in a coordinate y on P^1.  Seeded from a Klein-ansatz solution (klein6.py).

Special points: P[x][B] (anchor x in {0,1,i}, B nonempty subset of movers {a,b,c}): exactly the movers in B
equal x;  Q[B] (|B|>=2): exactly the movers in B coincide at a non-anchor value.
Mover m = kappa_m prod_{B ni m}(y - P[0][B]) / prod_{B ni m}(y - P[i][B]),  and m-1 vanishes on P[1][B].
"""
import numpy as np, itertools
from klein6 import movers_t, collisions

MOV = ['a', 'b', 'c']
ANC = ['0', '1', 'i']
SUBS = [frozenset(s) for r in (1, 2, 3) for s in itertools.combinations(MOV, r)]
QSUB = [S for S in SUBS if len(S) >= 2]
KEYS = [(x, S) for x in ANC for S in SUBS] + [('q', S) for S in QSUB]   # 25
KIDX = {k: n for n, k in enumerate(KEYS)}
NU = 25 + 3

def from_klein(u, zc):
    M = movers_t(u, zc)
    names, ev = collisions(M)
    y = np.zeros(25, complex)
    for t0, big in ev:
        assert len(big) == 1
        g = set(big[0])
        anc = [x for x in ANC if x in g]
        movs = frozenset(m for m in MOV if m in g)
        if anc:
            assert len(anc) == 1
            y[KIDX[(anc[0], movs)]] = t0
        else:
            y[KIDX[('q', movs)]] = t0
    kap = []
    for m in MOV:
        N, D = M[m]
        N = np.trim_zeros(N, 'f'); D = np.trim_zeros(D, 'f')
        kap.append(N[0]/D[0])
    return np.concatenate([y, kap]), M

def mover_polys(v):
    y = v[:25]; kap = v[25:28]
    out = {}
    for n, m in enumerate(MOV):
        q = {x: np.poly([y[KIDX[(x, S)]] for S in SUBS if m in S]) for x in ANC}
        out[m] = (kap[n], q)
    return out

def mval(mp, m, w):
    k, q = mp[m]
    return k*np.polyval(q['0'], w)/np.polyval(q['i'], w)

def residual6(v):
    mp = mover_polys(v)
    R = []
    for m in MOV:
        k, q = mp[m]
        dep = k*q['0'] - q['i'] - (k-1)*q['1']
        R += list(dep[1:])
    y = v[:25]
    for S in QSUB:
        w = y[KIDX[('q', S)]]
        ms = sorted(S)
        for m2 in ms[1:]:
            R.append(mval(mp, ms[0], w) - mval(mp, m2, w))
    return np.array(R)

def cluster_values(v):
    """value of the cluster at each of the 25 special points (inf encoded as None)."""
    mp = mover_polys(v); y = v[:25]
    vals = []
    for (x, S) in KEYS:
        if x == '0': vals.append(0.0)
        elif x == '1': vals.append(1.0)
        elif x == 'i': vals.append(None)
        else: vals.append(mval(mp, sorted(S)[0], y[KIDX[(x, S)]]))
    return vals

if __name__ == '__main__':
    u = np.load('klein6_sol_3.npy'); zc = np.load('klein6_zc.npy')
    v, M = from_klein(u, zc)
    print('residual of Klein seed in general formulation:', np.linalg.norm(residual6(v)))
    # Jacobian rank -> dimension of the family at this point (minus PGL_2 = 3)
    h = 1e-7; n = len(v); F0 = residual6(v)
    J = np.zeros((len(F0), n), complex)
    for k in range(n):
        dv = np.zeros(n, complex); dv[k] = h
        J[:, k] = (residual6(v+dv) - residual6(v-dv))/(2*h)
    s = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(s > 1e-6*s[0]))
    print('equations', len(F0), 'unknowns', n, 'Jacobian rank', rank, ' local dim =', n - rank, ' minus PGL2 =', n - rank - 3)
    np.save('gen6_seed.npy', v)
