"""Evidence: does adding pullbacks from ALL pairs and triples of 'cheap' directions (columns and class-combination
directions m in {-1,0,1}^8 with at most 8 values on Y_3(2)) close the certified degree-7 deficit of J(Y_3(2))?
Generators exact (leading forms of 2D/3D images over Q, mapped to y-coordinates with integer coefficients); dims mod p
(so a remaining deficit is evidence only; closure would be certified by the sandwich)."""
import itertools, sys
from math import gcd
sys.path.insert(0, '.')
from lfi import monos, leading_forms, ideal_dim_modp, hilbert
from circuits_loc import section
from fano_block_loc import chi
import flint
j, D = 2, 7
F = list(range(1, 8)); n = 7; cap = 4
P = section(3, F, j)
HP = [1, 7, 28, 84, 210, 455, 812, 756]   # exact (fano_block_cert_out.txt)
dirs, seen = [], set()
for q in range(7):
    v = [0] * 7; v[q] = 1; dirs.append(tuple(v)); seen.add(tuple(v))
for m in itertools.product((-1, 0, 1), repeat=8):
    coef = [-sum(m[z] * chi(a, z) for z in range(8)) for a in F]
    if all(c == 0 for c in coef):
        continue
    g = 0
    for c in coef:
        g = gcd(g, abs(c))
    prim = tuple(c // g for c in coef)
    if prim[next(i for i, c in enumerate(prim) if c)] < 0:
        prim = tuple(-c for c in prim)
    if prim in seen:
        continue
    seen.add(prim)
    if len(set(sum(c * y for c, y in zip(prim, pt)) for pt in P)) <= 8:
        dirs.append(prim)
print('cheap directions:', len(dirs)); sys.stdout.flush()

def lin_pow_expand(form_coeffs_list, e):
    """product prod_q (sum_a V[q][a] y_a)^{e_q} as dict"""
    poly = {tuple([0] * n): 1}
    for V, k in form_coeffs_list:
        for _ in range(k):
            new = {}
            for ex, c in poly.items():
                for a, va in enumerate(V):
                    if va:
                        f = list(ex); f[a] += 1; f = tuple(f)
                        new[f] = new.get(f, 0) + c * va
            poly = new
    return poly

def family_gens(G):
    img = sorted(set(tuple(sum(c * y for c, y in zip(v, pt)) for v in G) for pt in P))
    gens = []
    b = len(G)
    for i in range(2, D + 1):
        top, basis = leading_forms(img, i)
        for vec in basis:
            poly = {}
            for e, c in zip(top, vec):
                if c:
                    t2 = {tuple([0] * n): 1}
                    for q in range(b):
                        t2 = mul(t2, power(G[q], e[q]))
                    for ex, cc in t2.items():
                        poly[ex] = poly.get(ex, 0) + c * cc
            poly = {ex: c for ex, c in poly.items() if c and max(ex) <= cap}
            if poly:
                gens.append(poly)
    return gens

_pc = {}
def power(V, k):
    key = (V, k)
    if key not in _pc:
        poly = {tuple([0] * n): 1}
        for _ in range(k):
            new = {}
            for ex, c in poly.items():
                for a, va in enumerate(V):
                    if va:
                        f = list(ex); f[a] += 1; f = tuple(f)
                        new[f] = new.get(f, 0) + c * va
            poly = new
        _pc[key] = poly
    return _pc[key]

def mul(p1, p2):
    out = {}
    for e1, c1 in p1.items():
        for e2, c2 in p2.items():
            f = tuple(a + b for a, b in zip(e1, e2))
            out[f] = out.get(f, 0) + c1 * c2
    return out

for b in (1, 2, 3):
    gens = []
    for G in itertools.combinations(dirs, b):
        gens += family_gens(list(G))
    # degree 6,7 comparison mod p (reduced monomials; grid implicit)
    out = []
    for i in (6, 7):
        Ni = len(monos(n, i, cap))
        out.append((i, Ni - HP[i] - ideal_dim_modp(i, n, gens, cap=cap)))
    print(f'families of <= {b} cheap directions: deficits (degree, dim J - dim_p L) = {out}'); sys.stdout.flush()
