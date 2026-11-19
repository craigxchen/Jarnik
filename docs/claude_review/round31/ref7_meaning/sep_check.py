# Own small check of the VII.1 separation device and of the
# "one parent, many same-degree right composition factors" identity,
# constant model, p = q = 7 (7 = 3 mod 4), k = 2, F_49 = F_7(i).
# Polynomials in t over F_49: pair (re, im) of int64 numpy arrays, low degree first.
import sys, itertools
import numpy as np
P = 7

def tr(a):
    re, im = a
    n = max(len(re), len(im))
    re = np.pad(re % P, (0, n - len(re))); im = np.pad(im % P, (0, n - len(im)))
    nz = np.nonzero(re | im)[0]
    if len(nz) == 0:
        return (np.zeros(1, np.int64), np.zeros(1, np.int64))
    d = nz[-1] + 1
    return (re[:d].copy(), im[:d].copy())

def deg(a):
    a = tr(a)
    return -1 if (len(a[0]) == 1 and a[0][0] == 0 and a[1][0] == 0) else len(a[0]) - 1

def add(a, b):
    n = max(len(a[0]), len(b[0]))
    return tr((np.pad(a[0], (0, n - len(a[0]))) + np.pad(b[0], (0, n - len(b[0]))),
               np.pad(a[1], (0, n - len(a[1]))) + np.pad(b[1], (0, n - len(b[1])))))

def neg(a): return tr((-a[0], -a[1]))
def sub(a, b): return add(a, neg(b))
def conj(a): return tr((a[0], -a[1]))

def mul(a, b):
    ar, ai = a; br, bi = b
    c = lambda x, y: np.convolve(x, y) % P
    return tr((c(ar, br) - c(ai, bi), c(ar, bi) + c(ai, br)))

def const(re, im=0): return tr((np.array([re], np.int64), np.array([im], np.int64)))
def mono(d, re=1, im=0):
    r = np.zeros(d + 1, np.int64); i = np.zeros(d + 1, np.int64); r[d] = re; i[d] = im
    return tr((r, i))

def pw(a, e):
    r = const(1)
    while e:
        if e & 1: r = mul(r, a)
        a = mul(a, a); e >>= 1
    return r

def inv49(re, im):
    nrm = (re * re + im * im) % P
    ni = pow(int(nrm), P - 2, P)
    return (re * ni) % P, (-im * ni) % P

def divmod_(a, b):
    a = tr(a); b = tr(b)
    db = deg(b); lr, li = inv49(int(b[0][db]), int(b[1][db]))
    q = [np.zeros(max(deg(a) - db + 1, 1), np.int64), np.zeros(max(deg(a) - db + 1, 1), np.int64)]
    r = [a[0].copy(), a[1].copy()]
    for k in range(deg(a) - db, -1, -1):
        cr, ci = int(r[0][k + db]) % P, int(r[1][k + db]) % P
        if cr == 0 and ci == 0: continue
        fr, fi = (cr * lr - ci * li) % P, (cr * li + ci * lr) % P
        q[0][k] = fr; q[1][k] = fi
        br_, bi_ = b
        r[0][k:k + db + 1] = (r[0][k:k + db + 1] - (fr * br_ - fi * bi_)) % P
        r[1][k:k + db + 1] = (r[1][k:k + db + 1] - (fr * bi_ + fi * br_)) % P
    return tr(tuple(q)), tr(tuple(r))

def gcd(a, b):
    a = tr(a); b = tr(b)
    while deg(b) >= 0:
        a, b = b, divmod_(a, b)[1]
    return deg(a)

def deriv(a):
    a = tr(a); n = len(a[0])
    if n == 1: return const(0)
    k = np.arange(1, n)
    return tr((a[0][1:] * k, a[1][1:] * k))

squares = {(x * x) % P for x in range(1, P)}
def chi(x):
    x %= P
    return 0 if x == 0 else (1 if x in squares else -1)

# ----- set-up
s = int(sys.argv[1]) if len(sys.argv) > 1 else 300
forms = [(1, 0), (0, 1), (1, 1)] if len(sys.argv) <= 2 else [tuple(map(int, f.split(','))) for f in sys.argv[2:]]
C = sub(mono(7), mono(1))                 # t^7 - t  (Delta = -1 constant model)
C7 = pw(C, 7)
assert deg(sub(C7, sub(mono(49), mono(7)))) == -1
labels = list(itertools.product(range(P), repeat=2))
def label(x):
    return mul(mul(const(0, 1), C7), add(const(x[0]), mono(7, x[1])))   # i C^7 (x0 + x1 t^7)
V = {x: label(x) for x in labels}
assert all(deg(deriv(v)) == -1 for v in V.values())        # dv = 0
if s % P == 0:
    T = add(const(1), mul(C7, add(mono(s), mono(1))))
else:
    T = add(const(1), mul(C7, mono(s)))
L = deg(T); b = max(deg(v) for v in V.values())
print("q=p=7 k=2 n=49 h=7; b=%d L=%d" % (b, L))
dT = deriv(T)
# simple zeros and disjointness
ok1 = all(gcd(sub(T, v), dT) == 0 for v in V.values())
ok2 = all(gcd(sub(T, V[x]), sub(T, V[y])) == 0 for x, y in itertools.combinations(labels, 2))
print("each T-v squarefree:", ok1, " pairwise coprime:", ok2)
R = const(1)
for v in V.values(): R = mul(R, sub(T, v))
print("deg R =", deg(R), " R in F_7[t]:", not np.any(R[1]))

# ----- composition identity P_V(X) = A_j(X)^q - c_j^{q-1} A_j(X), tested at X = T
PV_T = R
for (a1, a2) in forms:
    ker = [x for x in labels if (a1 * x[0] + a2 * x[1]) % P == 0]
    A_T = const(1)
    for x in ker: A_T = mul(A_T, sub(T, V[x]))
    x0 = next(x for x in labels if (a1 * x[0] + a2 * x[1]) % P == 1)
    c = const(1)
    for x in ker: c = mul(c, sub(V[x0], V[x]))
    rhs = sub(pw(A_T, 7), mul(pw(c, 6), A_T))
    print("form", (a1, a2), ": P_V(T) == A^7 - c^6 A :", deg(sub(PV_T, rhs)) == -1, " deg A(T) =", deg(A_T))

# ----- marks
ex = {}
for j, (a1, a2) in enumerate(forms):
    for x in labels: ex[(j, x)] = 1 + chi(a1 * x[0] + a2 * x[1])
mn = {x: min(ex[(j, x)] for j in range(len(forms))) for x in labels}
g = sum(mn.values())
Tv = {x: sub(T, V[x]) for x in labels}
Z = []
for j in range(len(forms)):
    z = const(1)
    for x in labels:
        e = ex[(j, x)] - mn[x]
        if e: z = mul(z, pw(Tv[x], e))
    Z.append(z)
D = deg(Z[0])
norms = [mul(z, conj(z)) for z in Z]
print("g =", g, " D =", D, " equal degrees:", all(deg(z) == D for z in Z),
      " one norm:", all(deg(sub(nn, norms[0])) == -1 for nn in norms),
      " norm in F_7[t]:", not np.any(norms[0][1]))
# parent exponents: label x contributes (T-v_x) with exponent e_j(x)+e_j(-x)
capexp = {}
for x in labels:
    mx = ((-x[0]) % P, (-x[1]) % P)
    capexp[x] = (ex[(0, x)] - mn[x]) + (ex[(0, mx)] - mn[mx])
print("parent exponent histogram on the T-v:", {e: sum(1 for x in labels if capexp[x] == e) for e in set(capexp.values())})
md = max(deg(sub(Z[i], Z[j])) for i, j in itertools.combinations(range(len(Z)), 2))
print("max deg(z_i-z_j) =", md, " D/2 =", D / 2, " critical:", 2 * md <= D, " slack D-2md =", D - 2 * md)
# exponent-one share of deg N (unpadded)
one = sum(L for x in labels if capexp[x] == 1); tot = sum(L * capexp[x] for x in labels)
print("fraction of deg N on exponent-one primes:", one / tot if tot else None)
# padding by a common NONREAL factor of degree <= slack: exponent-one primes, cluster preserved
sl = D - 2 * md
if sl > 0:
    # G = prod (t - a_r) with a_r in F_49 \ F_7, distinct, conjugates not used; coprime to R (R has no roots there? check)
    import random
    random.seed(1)
    F49=[(a,bb) for a in range(P) for bb in range(P)]
    def ev(f,al):
        r=(0,0)
        for k in range(deg(f),-1,-1):
            r=((r[0]*al[0]-r[1]*al[1]+int(f[0][k]))%P,(r[0]*al[1]+r[1]*al[0]+int(f[1][k]))%P)
        return r
    G = const(1); used = 0; picked = []
    while used + 2 <= sl and len(picked) < 60:
        c0=(random.randrange(P),random.randrange(P)); c1=(random.randrange(P),random.randrange(1,P))
        f = tr((np.array([c0[0],c1[0],1],np.int64),np.array([c0[1],c1[1],0],np.int64)))
        if any(ev(f,al)==(0,0) for al in F49): continue
        key=(c0,c1); ck=((c0[0],(-c0[1])%P),(c1[0],(-c1[1])%P))
        if key in picked or ck in picked: continue
        if gcd(f,R)>0 or gcd(f,G)>0 or gcd(conj(f),G)>0: continue
        picked.append(key); G=mul(G,f); used+=2
    Zp = [mul(G, z) for z in Z]
    Dp = deg(Zp[0]); mdp = max(deg(sub(Zp[i], Zp[j])) for i, j in itertools.combinations(range(len(Zp)), 2))
    print("padded by nonreal common G of degree", used, ": D' =", Dp, " max diff =", mdp, " critical:", 2 * mdp <= Dp,
          " exponent-one share of deg N' =", (one + 2 * used) / (tot + 2 * used))
