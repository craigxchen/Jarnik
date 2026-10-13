"""General closure (Def. 2.5 sub-closure) for relations of any weight, Pi an arbitrary sublattice of Z^r.
Certificates searched each round, relative to the current Pi (canonical coset representatives):
  * k=2, n=(1,1): sigma v + sigma' v' in Pi, sin^2 d + sin^2 d' > 1;  self n=2: 2 sigma v in Pi, sin^2 d > 1/2;
  * k=3, n_i in {1,2}, two of the three terms sharing a column (third looked up by coset), arc < pi.
Every certificate is re-checked from scratch (exact tightness, exact sigma, 50-digit arc sum,
lattice membership in the lattice of previously added vectors).
usage: python3 gclosure.py primes.json|random:M:b:K:seed M [w3]
"""
import sys, json, math, itertools
from collections import defaultdict
import mpmath
from rcore import gauss, Lprime, tight_list, gmul, Ds
from rclosure import Lattice, A_of
from eval_set import triple_list
from farkas import load
mpmath.mp.dps = 50

def canon(lat, vec):
    v = {k: x for k, x in vec.items() if x}
    for c in sorted(lat.piv):
        if c not in v: continue
        P = lat.piv[c]; p = P[c]; q = v[c] // p
        if q:
            for k, x in P.items():
                y = v.get(k, 0) - q * x
                if y: v[k] = y
                else: v.pop(k, None)
    return tuple(sorted(v.items()))

def closure(r, rels, k3=True):
    lat = Lattice(); log = []
    s2f = [float(rel[3]) for rel in rels]
    arc = [math.pi - 2 * math.asin(math.sqrt(x)) for x in s2f]
    def vecof(i, sgn):
        v, s, sg, s2 = rels[i]; return {c: sgn * sg * co for c, co in v}
    changed = True
    while changed:
        changed = False
        plus = {}; minus = defaultdict(list); fire = []
        for i in range(len(rels)):
            cp = canon(lat, vecof(i, 1))
            if not cp: continue                       # already in Pi
            cm = canon(lat, vecof(i, -1))
            if cp == cm and rels[i][3] * 2 > 1:       # self certificate
                fire.append(([(i, 2)], [i])); continue
            plus.setdefault(cp, []).append(i); minus[cm].append(i)
        for key, lst in plus.items():                 # sigma v_i + sigma' v_j in Pi  <=> canon(+i) = canon(-j)
            js = minus.get(key, ())
            if not js: continue
            bi = max(lst, key=lambda i: rels[i][3]); bj = max(js, key=lambda j: rels[j][3])
            if bi != bj and rels[bi][3] + rels[bj][3] > 1:
                fire.append(([(bi, 1), (bj, 1)], [bi, bj]))
            elif bi == bj:
                for i in lst:
                    for j in js:
                        if i != j and rels[i][3] + rels[j][3] > 1:
                            fire.append(([(i, 1), (j, 1)], [i, j])); break
                    else: continue
                    break
        if not fire and k3:
            live = [i for i in range(len(rels)) if canon(lat, vecof(i, 1))]
            T = [(i, n) for i in live for n in (1, 2) if n * arc[i] < math.pi]
            bycol = defaultdict(list); coset = defaultdict(list)
            for (i, n) in T:
                vec = {c: n * x for c, x in vecof(i, 1).items()}
                coset[canon(lat, {c: -x for c, x in vec.items()})].append((i, n))
                for c in vec: bycol[c].append((i, n))
            for c, lst in bycol.items():
                for (a, na), (b, nb) in itertools.combinations(lst, 2):
                    if a == b or na * arc[a] + nb * arc[b] >= math.pi: continue
                    u = defaultdict(int)
                    for cc, x in vecof(a, na).items(): u[cc] += x
                    for cc, x in vecof(b, nb).items(): u[cc] += x
                    for (t, nt) in coset.get(canon(lat, u), ()):
                        if t not in (a, b) and na * arc[a] + nb * arc[b] + nt * arc[t] < math.pi - 1e-9:
                            fire.append(([(a, na), (b, nb), (t, nt)], [a, b, t])); break
                if len(fire) > 200: break
        changed = False
        for terms, idxs in fire:
            if all(lat.member(vecof(i, 1)) for i in idxs): continue
            log.append((terms, [rels[i][0] for i in idxs]))
            for i in idxs: lat.insert(vecof(i, 1))
            changed = True
    return lat, log

def check(gp, L, rels, log):
    lat = Lattice(); ok = 0
    for terms, added in log:
        u = defaultdict(int); tot = mpmath.mpf(0)
        for idx, n in terms:
            v, s, _, _ = rels[idx]; A = A_of(gp, v); D = Ds(A)[s]
            assert D != 0 and L % abs(D) == 0
            B = gmul(A, A)
            for _ in range(s): B = gmul(B, (0, -1))
            sig = (B[1] > 0) - (B[1] < 0)
            Cz = A
            for _ in range(s): Cz = gmul(Cz, (1, -1))
            s2 = mpmath.mpf(Cz[1] ** 2) / (2 ** s * (A[0] ** 2 + A[1] ** 2))
            tot += n * (mpmath.pi - 2 * mpmath.asin(mpmath.sqrt(s2)))
            for c, co in v: u[c] += n * sig * co
        assert tot < mpmath.pi - mpmath.mpf(10) ** -30, 'arc'
        assert lat.member(dict(u)), 'membership'
        for v in added:
            lat.insert({c: co for c, co in v})
        ok += 1
    r = len(gp)
    return ok, sum(lat.member({k: 1}) for k in range(r)), sum(lat.member({k: 2}) for k in range(r))

if __name__ == '__main__':
    ps = load(sys.argv[1]); M = int(sys.argv[2]); w3 = len(sys.argv) > 3
    gp = [gauss(p) for p in ps]; L = Lprime(2 * M, set(ps))
    rels = tight_list(gp, L) + (triple_list(gp, L) if w3 else [])
    lat, log = closure(len(gp), rels)
    ok, ex, hf = check(gp, L, rels, log)
    print(json.dumps(dict(set=sys.argv[1], r=len(gp), tight=len(rels), w3=w3, certs=len(log), certs_checked=ok,
                          k3_certs=sum(1 for t, _ in log if len(t) == 3), exact_pinned=ex, half_pinned=hf,
                          free=len(gp) - hf)), flush=True)
