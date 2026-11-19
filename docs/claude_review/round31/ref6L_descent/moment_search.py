"""Search for affine pencils G_e = T + c_e (T real integer, c_e = x_e + i y_e small Gaussian
integers) with strict distinctness whose ten cofactor polynomials f_ij(T) have low degree.
deg f_ij <= 5-k for all pairs  <=>  Im star power sums of the roots equal for orders 1..k (II.8).
On such a template Delta ~ T^{-(k+1)} ~ P^{-(k+1)}, i.e. nu = 5-k.
  k=4 would give nu=1 (exactly Astra's boundary, exponent 4);  k=5 is the open strict linear system.
We enumerate y in the integer kernel of E (star sums of y equal), x in ker E diag(y) (k=2),
then test k=3,4,5 exactly.
"""
import itertools, sys
import numpy as np
import sympy as sp
from itertools import combinations

edges = list(combinations(range(5), 2))
Inc = np.array([[int(i in e) for e in edges] for i in range(5)])
E = Inc[:4] - Inc[4]

def intkernel(Mat):
    Ms = sp.Matrix(Mat)
    ns = Ms.nullspace()
    vecs = []
    for v in ns:
        den = sp.ilcm(*[sp.fraction(x)[1] for x in v])
        vecs.append([int(x*den) for x in v])
    # LLL-ish reduction through sympy is unavailable; return HNF-free basis
    return np.array(vecs, dtype=np.int64)

Ymax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
Xmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3

def strict(x, y):
    c = list(zip(x, y))
    s = set(c) | set((a, -b) for a, b in c)
    return all(b != 0 for b in y) and len(s) == 20

def imsum(x, y, k):
    z = (np.array(x, dtype=object) + 1j*0)  # placeholder unused
    vals = []
    for a, b in zip(x, y):
        w = (sp.Integer(a) + sp.I*sp.Integer(b))**k
        vals.append(int(sp.im(sp.expand(w))))
    return np.array(vals, dtype=object)

def impow(x, y, k):
    # integer Im((x+iy)^k) via binomial expansion
    out = []
    for a, b in zip(x, y):
        s = 0
        for j in range(1, k+1, 2):
            s += sp.binomial(k, j) * a**(k-j) * b**j * (-1)**((j-1)//2)
        out.append(int(s))
    return np.array(out, dtype=object)

def ok(x, y, k):
    v = impow(x, y, k)
    return all(sum(v[e] for e in range(10) if i in edges[e]) == sum(v[e] for e in range(10) if 4 in edges[e]) for i in range(4))

# enumerate y: vectors with entries in [-Ymax..Ymax]\{0} and equal star sums
best = {}
count = {2: 0, 3: 0, 4: 0, 5: 0}
seen_y = 0
rng = [v for v in range(-Ymax, Ymax+1) if v != 0]
# use structure: choose y on 6 free edges, solve rest via E y = 0 (4 equations); keep integral ones
Ys = []
free = None
for cols in combinations(range(10), 4):
    if abs(round(np.linalg.det(E[:, cols]))) == 1:
        free = [c for c in range(10) if c not in cols]
        piv = list(cols)
        break
Ep = E[:, piv].astype(float); Ef = E[:, free].astype(float)
Einv = np.linalg.inv(Ep)
for yf in itertools.product(rng, repeat=6):
    yp = -Einv @ (Ef @ np.array(yf))
    ypr = np.rint(yp)
    if np.max(np.abs(yp-ypr)) > 1e-9:
        continue
    if np.any(ypr == 0) or np.max(np.abs(ypr)) > Ymax:
        continue
    y = np.zeros(10, dtype=np.int64); y[free] = yf; y[piv] = ypr.astype(np.int64)
    # normalize: skip y and -y duplicates and global sign
    Ys.append(tuple(y))
print('number of y with equal star sums:', len(Ys))
res3 = []
for y in Ys:
    Ey = E * np.array(y)[None, :]
    K = sp.Matrix(Ey).nullspace()
    # integral basis
    B = []
    for v in K:
        den = sp.ilcm(*[sp.fraction(t)[1] for t in v])
        B.append([int(t*den) for t in v])
    B = np.array(B, dtype=np.int64)
    # enumerate small combinations; x modulo constant: force x[0]=0 via filter
    nb = B.shape[0]
    for coef in itertools.product(range(-2, 3), repeat=nb):
        x = np.array(coef) @ B
        if x[0] != 0 or np.max(np.abs(x)) > Xmax:
            continue
        if not strict(list(x), list(y)):
            continue
        count[2] += 1
        if ok(list(x), list(y), 3):
            count[3] += 1
            res3.append((tuple(int(a) for a in x), y))
            if ok(list(x), list(y), 4):
                count[4] += 1
                print('FOUR-MOMENT', x, y)
                if ok(list(x), list(y), 5):
                    count[5] += 1
                    print('FIVE-MOMENT (constant cofactors!)', x, y)
print('counts (strict, by highest Im-moment order satisfied):', count)
for r in res3[:10]:
    print('k=3 example x=', r[0], 'y=', r[1])
