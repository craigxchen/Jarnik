# Judge's independent check of Part VI.4 (ternary construction), own F_729 arithmetic.
# Usage: python3 ternary.py FORMS L [signs]
#   FORMS like "10,01,11,12" (coefficients (a0,a1) of lambda = a0 x0 + a1 x1 over F_3), q=27, k=2.
import sys, itertools
import numpy as np

P = 3
DEG = 6
QQ = P**DEG  # 729

def polymulmod(a, b, mod):  # base-3 digit lists, mod monic of degree 6 (list of 7 coeffs, low first)
    res = [0]*(2*DEG-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                res[i+j] = (res[i+j] + x*y) % P
    for d in range(len(res)-1, DEG-1, -1):
        c = res[d]
        if c:
            for k in range(DEG+1):
                res[d-DEG+k] = (res[d-DEG+k] - c*mod[k]) % P
    return res[:DEG]

def toint(a): return sum(c*P**i for i, c in enumerate(a))
def todig(x): return [(x//P**i) % P for i in range(DEG)]

# find a primitive modulus
def find_prim():
    for tail in itertools.product(range(P), repeat=DEG):
        if tail[0] == 0: continue
        mod = list(tail)+[1]
        x = [0, 1, 0, 0, 0, 0]
        cur = [1, 0, 0, 0, 0, 0]
        order = None
        for e in range(1, QQ):
            cur = polymulmod(cur, x, mod)
            if cur == [1, 0, 0, 0, 0, 0]:
                order = e; break
        if order == QQ-1:
            return mod
mod = find_prim()
EXP = np.zeros(2*(QQ-1), dtype=np.int64)
LOG = np.full(QQ, -1, dtype=np.int64)
cur = [1, 0, 0, 0, 0, 0]
for e in range(QQ-1):
    v = toint(cur); EXP[e] = v; EXP[e+QQ-1] = v; LOG[v] = e
    cur = polymulmod(cur, [0, 1, 0, 0, 0, 0], mod)
DIG = np.array([todig(x) for x in range(QQ)], dtype=np.int64)
W3 = np.array([P**i for i in range(DEG)], dtype=np.int64)
ADD = np.zeros((QQ, QQ), dtype=np.int64)
for a in range(QQ):
    ADD[a] = ((DIG[a][None, :] + DIG) % P) @ W3
NEG = ((-DIG) % P) @ W3

def mul_s(s, arr):  # scalar times array
    if s == 0: return np.zeros_like(arr)
    out = np.zeros_like(arr)
    nz = arr != 0
    out[nz] = EXP[LOG[arr[nz]] + LOG[s]]
    return out
def mul_el(a, b):
    if a == 0 or b == 0: return 0
    return int(EXP[LOG[a]+LOG[b]])
def pw(a, e):
    if a == 0: return 0
    return int(EXP[(LOG[a]*e) % (QQ-1)])
def add_arr(a, b): return ADD[a, b]

g = int(EXP[1])
I = pw(g, 182)
assert mul_el(I, I) == NEG[1]
F27 = [0] + [pw(g, 28*k) for k in range(26)]
CHI = {0: 0}
for k in range(26): CHI[pw(g, 28*k)] = 1 if k % 2 == 0 else -1
F9 = set([0] + [pw(g, 91*k) for k in range(8)])
assert all(pw(x, 27) == x for x in F27)

forms = [tuple(int(c) for c in f) for f in sys.argv[1].split(',')]
L = int(sys.argv[2])
signs = [int(s) for s in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1]*len(forms)
m = len(forms)

def smul3(a, x):  # a in {0,1,2} times x in field
    return 0 if a == 0 else (x if a == 1 else int(NEG[x]))

V = []
for x0 in F27:
    for x1 in F27:
        V.append((x0, x1, mul_el(I, x0), mul_el(I, x1)))  # v = c0 + c1 t
n = len(V); h = n//27
E = []
for f, sg in zip(forms, signs):
    ex = []
    for (x0, x1, c0, c1) in V:
        lam = int(ADD[smul3(f[0], x0), smul3(f[1], x1)])
        if sg == -1: lam = int(NEG[lam])
        ex.append(1 + CHI[lam])
    E.append(ex)
gmin = [min(E[j][a] for j in range(m)) for a in range(n)]
gdeg = sum(gmin)
print("n=%d h=%d m=%d formal g=%d  sum exps=%s" % (n, h, m, gdeg, [sum(e) for e in E]))

def prodpoly(exps):
    deg = sum(exps)*L
    P_ = np.zeros(deg+1, dtype=np.int64); P_[0] = 1; cur = 0
    for (x0, x1, c0, c1), e in zip(V, exps):
        for _ in range(e):
            # multiply by t^L - c0 - c1 t
            new = np.zeros(deg+1, dtype=np.int64)
            new[L:cur+L+1] = P_[:cur+1]
            new[:cur+1] = add_arr(new[:cur+1], mul_s(int(NEG[c0]), P_[:cur+1]))
            new[1:cur+2] = add_arr(new[1:cur+2], mul_s(int(NEG[c1]), P_[:cur+1]))
            P_ = new; cur += L
    return P_

Z = [prodpoly([E[j][a]-gmin[a] for a in range(n)]) for j in range(m)]
D = len(Z[0])-1
assert all(len(z)-1 == D for z in Z) and all(z[-1] == 1 for z in Z)
# coefficients in F_9?
inF9 = all(set(np.unique(z).tolist()) <= F9 for z in Z)
print("D=%d (=(n-g)L=%d), monic, coefficients in F_9: %s" % (D, (n-gdeg)*L, inF9))
# F_9 = a + b i
F9map = {}
for a in range(3):
    for b in range(3):
        F9map[int(ADD[a, smul3(b, I)])] = (a, b)
def split(z):
    A = np.array([F9map[int(c)][0] for c in z], dtype=np.int64)
    B = np.array([F9map[int(c)][1] for c in z], dtype=np.int64)
    return A, B
norms = []
for z in Z:
    A, B = split(z)      # z = A + iB, conj = A - iB, norm = A^2 + B^2
    Nn = (np.convolve(A, A) + np.convolve(B, B)) % 3
    norms.append(Nn)
# also direct check: conj via 27-power on F_729 coefficients equals A - iB
z0c = np.array([pw(int(c), 27) for c in Z[0]])
A0, B0 = split(Z[0]); assert all(F9map[int(c)] == (a, (-b) % 3) for c, a, b in zip(z0c, A0, B0))
print("norms equal:", all(np.array_equal(norms[0], x) for x in norms), " deg N =", len(norms[0])-1, " in F_3[t] by construction (A^2+B^2)")
maxd = 0
for i in range(m):
    for j in range(i+1, m):
        d = add_arr(Z[i], NEG[Z[j]])
        nzi = np.nonzero(d)[0]
        dd = int(nzi[-1]) if len(nzi) else -1
        maxd = max(maxd, dd)
        print("pair", i, j, "deg diff (formal G only) =", dd, " <= D/2 =", D/2, dd <= D/2)
print("distinct:", all(not np.array_equal(Z[i], Z[j]) for i in range(m) for j in range(i+1, m)))
# no point is a unit multiple of another's conjugate: check conj(Z_i) != Z_j (all monic)
cz = [np.array([pw(int(c), 27) for c in z]) for z in Z]
print("no conj coincidences:", all(not np.array_equal(cz[i], Z[j]) for i in range(m) for j in range(m)))

# full gcd of the reduced points (over F_729 tables; inputs are in F_9)
def trim(a):
    nz = np.nonzero(a)[0]
    return a[:nz[-1]+1] if len(nz) else a[:0]
def pgcd(a, b):
    a = trim(a.copy()); b = trim(b.copy())
    while len(b):
        inv_lb = int(EXP[(QQ-1-LOG[b[-1]]) % (QQ-1)])
        db = len(b)-1
        while len(a)-1 >= db and len(a):
            c = mul_el(int(a[-1]), inv_lb)
            s = len(a)-1-db
            a[s:] = add_arr(a[s:], NEG[mul_s(c, b)])
            a = trim(a)
        a, b = b, a
    return a
import os
if os.environ.get("NOGCD"): sys.exit(0)
G = Z[0]
for z in Z[1:]:
    G = pgcd(G, z)
print("residual gcd degree of reduced points:", len(G)-1)
print("full gcd of raw marks = formal g*L + residual =", gdeg*L + len(G)-1, "; actual D =", D-(len(G)-1))
