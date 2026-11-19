# Test VI.3: full hyperplane designs with translation B (conj(Y) = Y + B), Y = T + theta*B, B in V.
# For every side choice (DFS with prefix gcds) compute the full-gcd reduced degree D', criticality
# (with monic normalisation), and check:
#  I1: D' >= (1/2) sum_{variable pairs} h_v            (claimed for every instance; no criticality used)
#  I2: D' >= (Q/8 - 1/2) deg B                          (content one)
#  For critical instances additionally: margin identity sum delta <= s'D'/4 on the common-unit subfamily,
#  difference identity F_i - F_j = P_J(B) P_J(Y+a), and delta_ij >= sum_{odd p not | B} ord_p(P_J(B)) deg p.
import sys, random, itertools, math, json
from f4lib import *
from clusters import span

def coords(basis):
    # map v -> coordinate bitmask
    out = {0: 0}
    for i, b in enumerate(basis):
        for v, c in list(out.items()):
            out[v ^ b] = c | (1 << i)
    return out

def parity(x):
    return bin(x).count('1') & 1

def odd_primes_upto(n):
    return [p for p in range(2, 1 << (n + 1)) if deg(p) <= n and deg(p) % 2 == 1 and is_irred2(p)]

def ordp(a, p):
    k = 0
    while a:
        q, r = divmod2(a, p)
        if r:
            break
        a = q; k += 1
    return k

def run(basis, B, Ts, out, stats, check_primes=False, rng=None):
    d = len(basis); Q = 1 << d
    V = span(basis); co = coords(basis)
    cB = co[B]
    forms = [l for l in range(1, Q) if parity(l & cB) == 1]
    s = len(forms); assert s == Q // 2
    sides = {}
    for l in forms:
        sides[l] = [[v for v in V if parity(l & co[v]) == c] for c in (0, 1)]
    sideset = {(l, c): set(sides[l][c]) for l in forms for c in (0, 1)}
    # pairs
    pairs = []
    seen = set()
    for v in V:
        if v in seen:
            continue
        seen.add(v); seen.add(v ^ B)
        pairs.append(v)
    oddp = odd_primes_upto(8)
    for T in Ts:
        content = T
        for v in basis:
            content = gcd2(content, v)
        if content != 1:
            continue
        Y = (T, B)
        lin = {v: add4(Y, (v, 0)) for v in V}
        Fp = {(l, c): prod4([lin[v] for v in sides[l][c]]) for l in forms for c in (0, 1)}
        degF = d4(Fp[(forms[0], 0)])
        assert all(d4(x) == degF for x in Fp.values())
        lcs = {k: lc4(x) for k, x in Fp.items()}
        mon = {k: scal(x, INV[lcs[k]]) for k, x in Fp.items()}
        keys = list(Fp)
        ddeg = {}
        for a in keys:
            for b in keys:
                if a < b:
                    ddeg[(a, b)] = d4(add4(mon[a], mon[b]))
        h = {}
        for v in pairs:
            g = gcd2(T ^ v, B)
            h[v] = max(deg(T ^ v), deg(B)) - deg(g)
        # DFS over sides
        best = None
        def dfs(i, G, chosen, maxd, mind=10**9):
            nonlocal best
            if i == s:
                Dp = degF - d4(G)
                crit = (maxd - d4(G)) <= Dp // 2 and mind >= 0
                # variable pairs
                var = 0
                for v in pairs:
                    ins = [ (v in sideset[k]) for k in chosen ]
                    if any(ins) and not all(ins):
                        var += h[v]
                stats['inst'] += 1
                if 2 * Dp < var:
                    stats['I1_fail'] += 1
                    out.append(('I1_fail', T, B, chosen, Dp, var))
                if Dp < (Q / 8 - 0.5) * deg(B):
                    stats['I2_fail'] += 1
                stats['I1_tight'] = max(stats.get('I1_tight', 0), var / (2 * Dp) if Dp else 0)
                if crit:
                    stats['crit'] += 1
                    if best is None or Dp < best[0]:
                        best = (Dp, T, tuple(chosen), d4(G))
                    critical_checks(chosen, G, Dp)
                return
            l = forms[i]
            for c in (0, 1):
                k = (l, c)
                G2 = Fp[k] if G is None else gcd4(G, Fp[k])
                md = maxd; mi = mind
                for k2 in chosen:
                    a, b = (k, k2) if k < k2 else (k2, k)
                    md = max(md, ddeg[(a, b)])
                    mi = min(mi, ddeg[(a, b)])
                # prune: if even the best possible further gcd cannot rescue... (no pruning, full enumeration)
                dfs(i + 1, G2, chosen + [k], md, mi)
        def critical_checks(chosen, G, Dp):
            # common-unit subfamily
            byu = {}
            for k in chosen:
                byu.setdefault(lcs[k], []).append(k)
            sub = max(byu.values(), key=len)
            sp = len(sub)
            sumdelta = 0
            for a, b in itertools.combinations(sub, 2):
                Fa, Fb = Fp[a], Fp[b]
                diff = add4(Fa, Fb)
                g = gcd4(Fa, Fb)
                delta = d4(diff) - d4(g)
                sumdelta += delta
                # identity F_i - F_j = P_J(B) P_J(Y+a0)
                la, lb = a[0], b[0]
                J = [v for v in V if parity(la & co[v]) == 0 and parity(lb & co[v]) == 0]
                Sa = set(sides[la][a[1]]); Sb = set(sides[lb][b[1]])
                a0 = next(iter(Sa & Sb))
                PJB = 1
                for v in J:
                    PJB = clmul(PJB, B ^ v)
                PJY = prod4([add4(Y, (a0 ^ v, 0)) for v in J])
                if mul4((PJB, 0), PJY) != diff:
                    stats['ident_fail'] += 1
                iner = sum(ordp(PJB, p) * deg(p) for p in oddp if mod2(B, p) != 0)
                if delta < iner:
                    stats['inert_fail'] += 1
            if 4 * sumdelta > sp * Dp:
                stats['margin_fail'] += 1
            stats['margin_ratio'] = max(stats.get('margin_ratio', 0), 4 * sumdelta / (sp * Dp))
        dfs(0, None, [], -1)
        if best is not None:
            key = (tuple(basis), B)
            if key not in stats['best'] or best[0] < stats['best'][key][0]:
                stats['best'][key] = best

if __name__ == '__main__':
    which = sys.argv[1]
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    out = []
    stats = dict(inst=0, crit=0, I1_fail=0, I2_fail=0, ident_fail=0, inert_fail=0, margin_fail=0, best={})
    if which == 'd3':
        basis = [1, 2, 4]
        V = span(basis)
        Ts = list(range(0, 1 << 10))
        for B in [v for v in V if v]:
            run(basis, B, Ts, out, stats)
    elif which == 'd4':
        basis = [1, 2, 4, 8]
        Bs = [1, 2, 3, 4, 7, 6, 8, 11, 13, 15]
        Ts = []
        for L in range(6, 21):
            for _ in range(10):
                Ts.append((1 << L) ^ rng.randrange(1 << min(L, 7)))
        Ts += [1 << L for L in range(4, 21)]
        for B in Bs:
            run(basis, B, Ts, out, stats)
    elif which == 'd4B':
        basis = [1, 2, 4, 8]
        B = int(sys.argv[3]); Lmin, Lmax, low = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
        Ts = [(1 << L) ^ x for L in range(Lmin, Lmax + 1) for x in range(1 << low)]
        run(basis, B, Ts, out, stats)
    elif which == 'd4rand':
        # random 4-dim V in F_2[t]_{<7}, B random nonconstant element
        for trial in range(int(sys.argv[3])):
            while True:
                basis = [rng.randrange(1, 1 << 7) for _ in range(4)]
                if len(span(basis)) == 16 and len(set(span(basis))) == 16:
                    break
            V = span(basis)
            B = rng.choice([v for v in V if v])
            Ts = [(1 << L) ^ rng.randrange(1 << 7) for L in range(8, 22, 2)] + [rng.randrange(1 << 9) for _ in range(5)]
            run(basis, B, Ts, out, stats)
    best = {str(k): v for k, v in stats.pop('best').items()}
    print(json.dumps(stats))
    for k, v in sorted(best.items()):
        print("min critical D' for (basis,B)", k, "->", v[0], "T=", bin(v[1]), "gcd deg", v[3])
    for o in out[:10]:
        print(o)
