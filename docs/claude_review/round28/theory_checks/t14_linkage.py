"""T14: linkage form of (SHIFT_e) (mod p).  C = {-1,0,1}^d, B = C \ A, and for l in L_pair: C_l = C cap (C - l) (a grid),
B_l = C_l \ A_l, mu_{e_x} = k_x, mu_{e_x - e_y} = k_x k_y.  Gorenstein duality in R/J(C) (socle degree 2d) predicts
   (SHIFT_e)  <=>  J(B)_{2d-e}  =  sum_l mu_l J(B_l)_{2d-e-deg mu_l}   (mod J(C)),
and dim J(B)_{2d-e}/J(C) = H_A(e).  We compute both sides (C-standard monomials; J(S)_D = M_D(S)^perp) and also the
linkage identities H_{A_l}(e) = T_{C_l}(e) - H_{B_l}(s_l - e).  Usage: model L e"""
import sys, time
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/theory_checks')
from tlib import *
model, L, e = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
p = P1
n, pats = syl(3) if model == 'syl8' else sdc(7)
d = n - 1
A = build_box(n, [(S, L) for S in pats]); Aset = set(A)
C = list(itertools.product((-1, 0, 1), repeat=d))
B = [k for k in C if k not in Aset]
def Jperp_rows(S, D):
    """basis rows (over smonos(d,D)) of J(S)_D mod J(C) = annihilator of M_D(S)"""
    N = T(d, D)
    M = Me_modp(S, D, p, d) if S else None
    if M is None:
        return [[1 if i == j else 0 for j in range(N)] for i in range(N)]
    X, nul = M.transpose().nullspace()
    return [[int(X[i, j]) for i in range(N)] for j in range(nul)]
D = 2 * d - e
JB = Jperp_rows(B, D)
print(model, "L=%d e=%d: dim J(B)_%d mod J(C) = %d ; H_A(e) = %d" % (L, e, D, len(JB), H_modp(A, e, p, d)), flush=True)
idxD = {a: i for i, a in enumerate(smonos(d, D))}
rows = []
t0 = time.time()
for l in Lpair_half(d):
    supp = [x for x in range(d) if l[x]]
    Cl = shiftset(C, l)
    Al = shiftset(A, l, Aset); Als = set(Al)
    Bl = [k for k in Cl if k not in Als]
    Dl = D - len(supp)
    # J(B_l)_{Dl} mod J(C): annihilator of M_{Dl}(B_l) within C-standard monomials; but we need it mod J(C_l):
    # monomials with exponent 2 at x in supp are in J(C_l); they are killed by mu_l (k_x^3 in J(C)), so harmless.
    Jl = Jperp_rows(Bl, Dl)
    mons = smonos(d, Dl)
    for r in Jl:
        v = [0] * T(d, D)
        for c, a in zip(r, mons):
            if c == 0: continue
            b = list(a)
            ok = True
            for x in supp:
                b[x] += 1
                if b[x] > 2: ok = False
            if ok: v[idxD[tuple(b)]] = (v[idxD[tuple(b)]] + c) % p
        if any(v): rows.append(v)
    # linkage identity for A_l in the grid C_l (socle degree s_l = 2d - |supp|)
    if l == Lpair_half(d)[0] or l == Lpair_half(d)[d]:
        sl = 2 * d - len(supp)
        ok = all(H_modp(Al, t, p, d) == H_modp(Cl, t, p, d) - (H_modp(Bl, sl - t, p, d) if Bl else 0) for t in range(0, sl + 1))
        print("  linkage identity in C_l for l=%s: %s" % (l, ok), flush=True)
R = flint.nmod_mat(len(rows), T(d, D), [x for r in rows for x in r], p)
rk = R.rank()
both = flint.nmod_mat(len(rows) + len(JB), T(d, D), [x for r in rows + JB for x in r], p).rank()
print("  dim sum_l mu_l J(B_l) = %d ; dim (that + J(B)) = %d ; dim J(B) = %d  -> containment: %s ; linkage-SHIFT_e: %s (%.0fs)" % (
      rk, both, len(JB), both == len(JB), rk == len(JB), time.time() - t0))
