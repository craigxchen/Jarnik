"""T3: first-order derivative closure.  For f in I(A): tau_l f - f in I(A_l), leading form d_l LF(f).  So
J1(A,l) := J(A) + (d_l J(A)) + (k_x^2 : l_x != 0)  is contained in J(A_l).  Test equality degree by degree (mod p):
dim J1(A,l)_e  vs  dim J(A_l)_e = T_d(e) - H_{A_l}(e).  Usage: t3_firstorder.py model L emax"""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L, emax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
p = P1
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
A = build_box(n, [(S, L) for S in pats]); Aset = set(A)
def Jbasis(S, e):
    """rows = basis of J(S)_e (std forms, coefficient vectors over smonos(d,e)) mod p"""
    M = Me_modp(S, e, p, d)
    N = T(d, e)
    if M is None:
        return [[1 if i == j else 0 for j in range(N)] for i in range(N)]
    X, nul = M.transpose().nullspace()      # vectors v with v . (columns of M) = 0
    return [[int(X[i, j]) for i in range(N)] for j in range(nul)]
JA = {e: Jbasis(A, e) for e in range(0, emax + 2)}
def deriv(vec, e, l):
    """d_l of a std form of degree e (vector over smonos(d,e)) -> std form of degree e-1 (dict)"""
    out = {}
    for c, a in zip(vec, smonos(d, e)):
        if c == 0: continue
        for x in range(d):
            if l[x] and a[x]:
                b = list(a); b[x] -= 1; b = tuple(b)
                out[b] = (out.get(b, 0) + c * a[x] * l[x]) % p
    return out
def mult(form, mono):
    out = {}
    for a, c in form.items():
        b = tuple(i + j for i, j in zip(a, mono))
        if max(b) <= 2: out[b] = (out.get(b, 0) + c) % p
    return out
def rank_rows(rows, N):
    if not rows: return 0
    return flint.nmod_mat(len(rows), N, [x % p for r in rows for x in r], p).rank()
for name, l in [('e0', (1, 0, 0, 0, 0, 0, 0)), ('e0-e1', (1, -1, 0, 0, 0, 0, 0))]:
    Al = shiftset(A, l, Aset)
    for e in range(0, emax + 1):
        N = T(d, e); idx = {a: i for i, a in enumerate(smonos(d, e))}
        rows = [list(r) for r in JA[e]]
        for f in range(1, e + 2):          # d_l J(A)_f lands in degree f-1, times R_{e-f+1}
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
        r1 = rank_rows(rows, N)
        HAl = H_modp(Al, e, p, d)
        print("%s L=%d l=%s e=%d: dim J1 = %d, dim J(A_l) = %d, dim J(A) = %d  %s" % (model, L, name, e, r1, N - HAl, len(JA[e]),
              "EQUAL" if r1 == N - HAl else "J1 SMALLER by %d" % (N - HAl - r1)), flush=True)
