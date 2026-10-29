"""T4: syl8 (L=2), e=6: is the invariant obstruction F (t2) already in the FIRST-ORDER closures J1(A,l)_6 for
l = e0, e0-e1?  (mod p).  Reports dims of invariant parts of J(A)_6, J1(A,l)_6, J(A_l)_6."""
import sys, time, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
exec(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/t3_firstorder.py').read().split("for name, l in")[0].replace("model, L, emax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])", "model, L, emax = 'syl8', 2, 6"))
D = pickle.load(open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks/t2_syl8_L2_e6.pkl', 'rb'))
orbits = D['orbits']; e = 6
N = T(d, e); idx = {a: i for i, a in enumerate(smonos(d, e))}
Fvecs = []
for orb in orbits:
    v = [0] * N
    for a in orb: v[idx[a]] = 1
    Fvecs.append(v)
def J1rows(l, e):
    rows = [list(r) for r in JA[e]]
    for f in range(1, e + 2):
        for vec in JA[f]:
            g = deriv(vec, f, l)
            if not g: continue
            for mono in smonos(d, e - f + 1):
                h = mult(g, mono)
                if h:
                    r = [0] * N
                    for a, c in h.items(): r[idx[a]] = c
                    rows.append(r)
    for x in range(d):
        if l[x]:
            for mono in smonos(d, e - 2):
                b = list(mono); b[x] += 2
                if max(b) <= 2:
                    r = [0] * N; r[idx[tuple(b)]] = 1; rows.append(r)
    return rows
def inv_part(rows):
    """dimension of {a in F_p^8 : sum a_i F_i in span(rows)} and a basis"""
    R = flint.nmod_mat(len(rows), N, [x % p for r in rows for x in r], p)
    Rb, rk = R.rref()
    rowsb = [[int(Rb[i, j]) for j in range(N)] for i in range(rk)]
    # columns: rowsb^T (N x rk) | -Fvecs^T (N x 8): nullspace -> projections onto last 8 coords
    M = flint.nmod_mat(N, rk + 8, [ (rowsb[j][i] if j < rk else (-Fvecs[j - rk][i]) % p) for i in range(N) for j in range(rk + 8)], p)
    X, nul = M.nullspace()
    vecs = [[int(X[rk + t, c]) for t in range(8)] for c in range(nul)]
    if not vecs: return 0, []
    V = flint.nmod_mat(len(vecs), 8, [x for v in vecs for x in v], p)
    Vb, r = V.rref()
    return r, [[int(Vb[i, j]) for j in range(8)] for i in range(r)]
dJA, bJA = inv_part([list(r) for r in JA[e]])
print("dim J(A)_6^G =", dJA)
for name, l in [('e0', (1, 0, 0, 0, 0, 0, 0)), ('e0-e1', (1, -1, 0, 0, 0, 0, 0))]:
    Al = shiftset(A, l, Aset)
    d1, b1 = inv_part(J1rows(l, e))
    # J(A_l)_6 rows
    M = Me_modp(Al, e, p, d)
    X, nul = M.transpose().nullspace()
    JAl = [[int(X[i, j]) for i in range(N)] for j in range(nul)]
    d2, b2 = inv_part(JAl)
    print("l=%s: dim J1(A,l)_6^G = %d, dim J(A_l)_6^G = %d" % (name, d1, d2), flush=True)
    # is the obstruction (t2 representative) in J1?  rep coefficients: orbit index of reps
    print("   J1^G basis:", b1)
    print("   J(A_l)^G basis:", b2)
print("J(A)^G basis:", bJA)
