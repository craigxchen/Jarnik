# Local solubility of depth sets D: exist integers a, y_d (d in D) with n = a^2+y_0^2 and
# n - (a-d)^2 = y_d^2 for all d in D (d=0 included).  Necessary: modulo M, there are residues
# a, n such that every n-(a-d)^2 is a square mod M, and n is a sum of two squares mod M.
# (We test a few prime-power moduli; passing is necessary, not sufficient.)
import itertools, sys
def squares(M): return {x*x % M for x in range(M)}
MODS = [16, 32, 9, 27, 5, 25, 7, 49, 11, 13, 17]
SQ = {M: squares(M) for M in MODS}
def locally_ok(D):
    for M in MODS:
        sq = SQ[M]
        ok = False
        for a in range(M):
            for n in range(M):
                if all((n - (a - d)**2) % M in sq for d in D):
                    ok = True; break
            if ok: break
        if not ok: return M
    return 0
res = []
for m in (3, 4, 5, 6):
    for dmax in range(m-1, 13):
        for mid in itertools.combinations(range(1, dmax), m-2):
            D = (0,) + mid + (dmax,)
            ob = locally_ok(D)
            if ob == 0:
                res.append((m, dmax, D))
for m in (3, 4, 5, 6):
    L = [r for r in res if r[0] == m]
    print("m=%d locally soluble depth sets with smallest dmax:" % m, [r[2] for r in L if r[1] == min(x[1] for x in L)] if L else None, " count(dmax<=12)=%d" % len(L))
import pickle; pickle.dump(res, open('local_ok.pkl', 'wb'))
