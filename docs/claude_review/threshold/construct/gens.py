"""point-set generators (independent)."""
import itertools

def Qslice(M, p, n):
    """Q^M(p,n) = {k in Z^M: sum k = p-n, |k|_1 <= p+n}"""
    s, t = p - n, p + n
    out = []
    def rec(pref, used):
        if len(pref) == M - 1:
            last = s - sum(pref)
            if used + abs(last) <= t:
                out.append(tuple(pref) + (last,))
            return
        for a in range(-(t - used), t - used + 1):
            rec(pref + [a], used + abs(a))
    rec([], 0)
    return out

def orbit(v):
    return sorted(set(itertools.permutations(v)))

def SU(M, s, U, box):
    """{k in Z^M : sum k = s, k(P) <= U(P) for all proper nonempty P}, U: function of the frozenset P.
    enumerated inside the box [-box, box]^M."""
    out = []
    subsets = [frozenset(c) for r in range(1, M) for c in itertools.combinations(range(M), r)]
    Uv = {P: U(P) for P in subsets}
    for pref in itertools.product(range(-box, box + 1), repeat=M - 1):
        last = s - sum(pref)
        if abs(last) > box:
            continue
        k = pref + (last,)
        ok = True
        for P in subsets:
            if sum(k[j] for j in P) > Uv[P]:
                ok = False
                break
        if ok:
            out.append(k)
    return out

def SUsym(M, s, U):
    """{k in Z^M : sum k = s, k(P) <= U[|P|-1] for all proper nonempty P} (U symmetric profile).
    Uses: k(P) <= U_j for all |P|=j  <=>  sum of the j largest entries <= U_j."""
    lo, hi = s - U[M - 2], U[0]
    out = []
    def rec(pref):
        if len(pref) == M - 1:
            last = s - sum(pref)
            if last < lo or last > hi:
                return
            k = pref + [last]
            ks = sorted(k, reverse=True)
            acc = 0
            for j in range(M - 1):
                acc += ks[j]
                if acc > U[j]:
                    return
            out.append(tuple(k))
            return
        for a in range(lo, hi + 1):
            pref2 = pref + [a]
            ks = sorted(pref2, reverse=True)
            acc = 0
            ok = True
            for j in range(len(ks)):
                acc += ks[j]
                if acc > U[j]:
                    ok = False
                    break
            if ok:
                rec(pref2)
    rec([])
    return out
