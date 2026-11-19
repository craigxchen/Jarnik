# (b) Count fibres of the Moore-minor map over F_32 exhaustively (x_1 = 1 normalisation),
# and test the dimension-4 variant (5 hyperplanes of span(1,x2,x3,x4) avoiding 1) by image counting.
import itertools, sys
from gf import GF2k
F = GF2k(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
q = F.q
# log tables
exp = [0] * (2 * q); log = [0] * q
g = None
for cand in range(2, q):
    x = 1; seen = set()
    for i in range(q - 1):
        seen.add(x); x = F.mul(x, cand)
    if len(seen) == q - 1: g = cand; break
x = 1
for i in range(q - 1):
    exp[i] = x; exp[i + q - 1] = x; log[x] = i; x = F.mul(x, g)
def mul(a, b): return 0 if a == 0 or b == 0 else exp[log[a] + log[b]]
def inv(a): return exp[(q - 1 - log[a]) % (q - 1)]
def span(vs):
    out = [0]
    for v in vs: out = out + [u ^ v for u in out]
    return out
def prodnz(vs):
    r = 1
    for u in span(vs)[1:]: r = mul(r, u)
    return r

mode = sys.argv[2] if len(sys.argv) > 2 else "dim5"
if mode == "dim5":
    images = {}; dom = 0
    for x2, x3, x4, x5 in itertools.product(range(q), repeat=4):
        xs = (1, x2, x3, x4, x5)
        if len(set(span(xs))) != 32: continue  # need F_2-independence (only possible if q >= 32)
        dom += 1
        ys = [prodnz([xs[k] for k in range(5) if k != j]) for j in range(5)]
        i0 = inv(ys[0]); key = tuple(mul(i0, y) for y in ys)
        images[key] = images.get(key, 0) + 1
    print(f"F_{q}: projective source points (x1=1, independent): {dom}; distinct images [y]: {len(images)};"
          f" max fibre size: {max(images.values())}")
else:
    # dimension 4: V = span(1,x2,x3,x4); functionals l with l(1)=1 given by values on x2,x3,x4
    import random
    sel = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)] if mode == "dim4a" else \
          [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1)]
    images = set(); dom = 0
    for x2, x3, x4 in itertools.product(range(q), repeat=3):
        V = span((1, x2, x3, x4))
        if len(set(V)) != 16: continue
        Hs = []
        for lv in sel:
            # kernel of l: basis of H: elements v in V with l(v)=0; l(1)=1, l(x2)=lv0...
            H = []
            for c in itertools.product([0, 1], repeat=4):
                v = (c[0] and 1) ^ (x2 if c[1] else 0) ^ (x3 if c[2] else 0) ^ (x4 if c[3] else 0)
                val = c[0] ^ (c[1] & lv[0]) ^ (c[2] & lv[1]) ^ (c[3] & lv[2])
                if val == 0 and v != 0: H.append(v)
            Hs.append(H)
        Vset = set(V)
        for Y in range(q):
            if Y in Vset: continue
            dom += 1
            vals = []
            for H in Hs:
                r = 1
                for v in H: r = mul(r, Y ^ v)
                vals.append(r)
            i0 = inv(vals[0]); images.add(tuple(mul(i0, y) for y in vals[1:]))
    print(f"F_{q} {mode}: domain points (x2,x3,x4,Y): {dom}; distinct images in T^4: {len(images)};"
          f" q^3 = {q**3}, q^4 = {q**4}")
