# Exhaustive check of Thm 2.5(a) combinatorics: for every class structure on S (zero class size z,
# nonzero class sizes = integer partition of s-z), every 1<=L<=s-1, the max total degree of disjoint gadgets
# {S-single (nonzero class, cost1, deg2), pair (diff classes, cost0, deg2), hexpair (two nonzero diff classes, cost1, deg3)}
# with total cost <= L is >= tau = min(2 sigma, 2L+2t, s+L).  Also report whether equality holds (tightness).
from functools import lru_cache
def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0: yield (); return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n-p, p): yield (p,)+rest
def best(labels, L):
    labels = tuple(sorted(labels))
    @lru_cache(maxsize=None)
    def rec(state, budget):
        # state: tuple of counts per class index (0 = zero class)
        state = list(state)
        # pick first nonempty class i, coordinate x in it: either unused, or used in gadget
        i = next((k for k,c in enumerate(state) if c>0), None)
        if i is None: return 0
        res = -1
        st = state[:]; st[i]-=1
        res = max(res, rec(tuple(st), budget))           # unused
        if i != 0 and budget >= 1: res = max(res, 2 + rec(tuple(st), budget-1))  # S-single
        for k in range(len(state)):
            if k != i and st[k] > 0:
                st2 = st[:]; st2[k]-=1
                res = max(res, 2 + rec(tuple(st2), budget))  # pair
                if i != 0 and k != 0 and budget >= 1:
                    res = max(res, 3 + rec(tuple(st2), budget-1))  # hexpair
        return res
    return rec
bad = []; cnt = 0; tight = 0
for s in range(2, 10):
    for z in range(0, s+1):
        for part in partitions(s-z):
            state = (z,) + part
            sigma = s - z; M = max((z,)+part); t = s - M
            for L in range(1, s):
                tau = min(2*sigma, 2*L+2*t, s+L)
                rec = best(None, L) if False else None
                # build rec fresh
                from functools import lru_cache
                def mk():
                    @lru_cache(maxsize=None)
                    def rec(stt, budget):
                        stt = list(stt)
                        i = next((k for k,c in enumerate(stt) if c>0), None)
                        if i is None: return 0
                        st = stt[:]; st[i]-=1
                        r = rec(tuple(st), budget)
                        if i != 0 and budget >= 1: r = max(r, 2 + rec(tuple(st), budget-1))
                        for k in range(len(stt)):
                            if k != i and st[k] > 0:
                                st2 = st[:]; st2[k]-=1
                                r = max(r, 2 + rec(tuple(st2), budget))
                                if i != 0 and k != 0 and budget >= 1:
                                    r = max(r, 3 + rec(tuple(st2), budget-1))
                        return r
                    return rec
                val = mk()(state, L)
                cnt += 1
                if val < tau: bad.append((s, state, L, val, tau))
                if val == tau: tight += 1
print('instances', cnt, 'violations', len(bad), 'tight', tight)
print(bad[:10])
