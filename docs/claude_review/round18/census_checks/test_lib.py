"""Cross-check census_lib contacts against the min-formula of ptolemy.md Thm B(iv)."""
import random
from census_lib import *
def minformula(X, p, T, m):
    # labels 0..m; g_ij = v_p(det(P_i,P_j)), P_0=(1,0), P_i=(X_i,1): det(P_0,P_i)=1 -> g_0i=0
    def g(a, b):
        if a == b: raise
        if a == 0 or b == 0: return 0
        return vp(X[a-1]-X[b-1], p)
    best = None
    Tc = [x for x in range(0, m+1) if x not in T]
    for a, b in combinations(sorted(T), 2):
        for c, dd in combinations(Tc, 2):
            v = min(g(a,b)+g(c,dd)-g(a,c)-g(b,dd), g(a,b)+g(c,dd)-g(a,dd)-g(b,c))
            best = v if best is None else min(best, v)
    return max(0, best)
random.seed(1)
nbad = 0; ntest = 0
for trial in range(300):
    m = random.choice([4,5,6,7])
    base = random.choice([2,3,5,6,12,30])
    X = random.sample(range(-base**4, base**4), m)
    assert check_boolean_identity(X)
    l, delta = contacts_integral(X)
    ps = set()
    for i in range(m):
        for j in range(i+1, m): ps.update(factor_primes(X[i]-X[j]))
    for p in ps:
        for T in visible_blocks(m):
            ntest += 1
            a = l.get(p, {}).get(T, 0)
            if a != minformula(X, p, T, m): nbad += 1
print("min-formula comparisons:", ntest, "mismatches:", nbad)
