"""p = 7 (form A^2+B^2 over F_7(t), F_49 = F_7(i)): VI.4 with q = p = 7, k = 2, direct products, full gcd.
Pure-python F_49[t] arithmetic (pairs (a,b) = a + b i)."""
import itertools, sys
p = 7
def add(x, y): return ((x[0]+y[0]) % p, (x[1]+y[1]) % p)
def mul(x, y): return ((x[0]*y[0]-x[1]*y[1]) % p, (x[0]*y[1]+x[1]*y[0]) % p)
def inv(x):
    for a in range(p):
        for b in range(p):
            if mul(x, (a, b)) == (1, 0): return (a, b)
Z = (0, 0)
def trim(f):
    while f and f[-1] == Z: f.pop()
    return f
def pmul(f, g):
    r = [Z]*(len(f)+len(g)-1)
    for i, a in enumerate(f):
        if a == Z: continue
        for j, b in enumerate(g):
            if b != Z: r[i+j] = add(r[i+j], mul(a, b))
    return trim(r)
def psub(f, g):
    n = max(len(f), len(g)); f = f + [Z]*(n-len(f)); g = g + [Z]*(n-len(g))
    return trim([add(a, ((-b[0]) % p, (-b[1]) % p)) for a, b in zip(f, g)])
def pmod(f, g):
    f = f[:]; il = inv(g[-1])
    while len(f) >= len(g):
        c = mul(f[-1], il); s = len(f) - len(g)
        for j, b in enumerate(g):
            m = mul(c, b); f[s+j] = add(f[s+j], ((-m[0]) % p, (-m[1]) % p))
        trim(f)
    return f
def pdiv(f, g):
    f = f[:]; il = inv(g[-1]); q = [Z]*(len(f)-len(g)+1)
    while len(f) >= len(g):
        c = mul(f[-1], il); s = len(f) - len(g); q[s] = c
        for j, b in enumerate(g):
            m = mul(c, b); f[s+j] = add(f[s+j], ((-m[0]) % p, (-m[1]) % p))
        trim(f)
    assert not f
    return q
def gcd(f, g):
    while g: f, g = g, pmod(f, g)
    il = inv(f[-1]); return [mul(il, a) for a in f]
def conj(f): return [(a, (-b) % p) for a, b in f]
sq = {(x*x) % p for x in range(1, p)}
chi = lambda y: 0 if y % p == 0 else (1 if y % p in sq else -1)
assert chi(-1) == -1
forms = [(1, 0)] + [(a, 1) for a in range(p)]
def mark(lam, eps, L):
    f = [(1, 0)]
    for x in itertools.product(range(p), repeat=2):
        e = 1 + chi(eps*(lam[0]*x[0] + lam[1]*x[1]))
        fac = [Z]*(L+1); fac[L] = (1, 0)
        fac[0] = add(fac[0], (0, (-x[0]) % p)); fac[1] = add(fac[1], (0, (-x[1]) % p))
        for _ in range(e): f = pmul(f, fac)
    return f
best = {}
for L in range(2, int(sys.argv[1])+1):
    M = {(lam, e): mark(lam, e, L) for lam in forms for e in (1, -1)}
    keys = list(M)
    N0 = pmul(M[keys[0]], conj(M[keys[0]]))
    for kk in keys: assert pmul(M[kk], conj(M[kk])) == N0 and all(b == 0 for a, b in N0)
    for m in (3, 4):
        for S in itertools.combinations(keys, m):
            if len({s[0] for s in S}) < m: continue
            G = M[S[0]]
            for s in S[1:]: G = gcd(G, M[s])
            red = [pdiv(M[s], G) for s in S]
            D = len(red[0]) - 1
            ds = [len(psub(a, b)) - 1 for a, b in itertools.combinations(red, 2)]
            if min(ds) >= 0 and 2*max(ds) <= D:
                if m not in best or D < best[m][0]: best[m] = (D, L, S)
    print('L', L, 'best so far', {m: best[m][:2] for m in best}, flush=True)
print(best)
