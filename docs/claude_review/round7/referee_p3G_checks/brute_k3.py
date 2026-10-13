"""Brute force: all level-0 certificates (Lemma 2.4 relative to Pi_0 = {0}) with at most 3 terms,
multiplicities n_i in {1, 2}, among ALL tight relations of weight <= 2 (singles, doubles, pairs),
and also with tight weight-3 relations if 'w3' is given.  Arc sums are multiples of pi/2 (Lemma R.1),
so 'arc < pi' is tested as 'arc < 3pi/4' (robust).
usage: python3 brute_k3.py primes.json|random:M:b:K:seed M [w3]
"""
import sys, json, math, itertools
from collections import defaultdict
from rcore import gauss, Lprime, tight_list
from eval_set import triple_list
from farkas import load

def terms_of(rels):
    T = []
    for idx, (v, s, sg, s2) in enumerate(rels):
        w = math.pi - 2 * math.asin(math.sqrt(float(s2)))
        for n in (1, 2):
            vec = tuple(sorted((c, n * sg * co) for c, co in v))
            T.append((vec, n * w, idx, n))
    return T

def neg(vec): return tuple((c, -x) for c, x in vec)

def add(a, b):
    d = defaultdict(int)
    for c, x in a: d[c] += x
    for c, x in b: d[c] += x
    return tuple(sorted((c, x) for c, x in d.items() if x))

def search(rels):
    T = terms_of(rels); byvec = defaultdict(list); bycol = defaultdict(list)
    for t in T:
        byvec[t[0]].append(t)
        for c, _ in t[0]: bycol[c].append(t)
    found = []
    for t in T:                                       # k = 1 impossible; k = 2
        for u in byvec.get(neg(t[0]), ()):
            if u[2] != t[2] and t[1] + u[1] < 0.75 * math.pi: found.append((t[2:], u[2:]))
    seen = set()
    for c, lst in bycol.items():                      # k = 3: two terms share column c
        for a, b in itertools.combinations(lst, 2):
            if a[2] == b[2] or a[1] + b[1] >= 0.75 * math.pi: continue
            s = add(a[0], b[0])
            if not s: continue
            for u in byvec.get(neg(s), ()):
                if u[2] in (a[2], b[2]): continue
                if a[1] + b[1] + u[1] < 0.75 * math.pi:
                    key = tuple(sorted((a[2:], b[2:], u[2:])))
                    if key not in seen: seen.add(key); found.append(key)
    return found

if __name__ == '__main__':
    ps = load(sys.argv[1]); M = int(sys.argv[2]); w3 = len(sys.argv) > 3
    gp = [gauss(p) for p in ps]; L = Lprime(2 * M, set(ps))
    rels = tight_list(gp, L) + (triple_list(gp, L) if w3 else [])
    f = search(rels)
    print(json.dumps(dict(set=sys.argv[1], r=len(ps), tight=len(rels), w3=w3, level0_certs_k_le_3=len(f),
                          example=[[(rels[i][0], rels[i][1], n) for i, n in c] for c in f[:1]])), flush=True)
