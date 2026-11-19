"""Independent reconstruction of the five-mark integer star model (Astra III.11).
Usage: python3 five_model.py MODE   MODE = small (K~1e7, full factorisation) or crt (Astra's recipe K)."""
import sys, random
from math import gcd, prod, isqrt
from itertools import combinations
from fractions import Fraction
import sympy
from sympy import Matrix, factorint, isprime
from mpmath import mp, mpf, atan2, log as mlog
mp.dps = 120

s = [3, 7, 11, 19, 23]
u = [x*x for x in s]
pairs = list(combinations(range(5), 2))
F = [(0,1),(1,2),(2,3),(3,4),(0,4)]
Fset = set(F)

# --- 1. polynomial row (f,h): kernel of the 5 linear conditions, computed independently
z = sympy.symbols('z')
cols = []
rows = []
for (i,j) in F:
    w = -u[i]*u[j]
    c = u[i]+u[j]
    rows.append([1, w, w*w, c, c*w, c*w*w])
Mk = Matrix(rows)
ker = Mk.nullspace()
print("kernel dim", len(ker))
kv = ker[0]
den = sympy.ilcm(*[sympy.fraction(x)[1] for x in kv])
kv = [int(x*den) for x in kv]
g0 = 0
for x in kv: g0 = gcd(g0, x)
kv = [x//g0 for x in kv]
fA = (3566328745136157,-212247537740,-232553177)
hA = (-68417652775626,-13737560909,189877)
astra = list(fA)+list(hA)
ratio = Fraction(astra[0], kv[0])
print("my primitive kernel vector:", kv)
print("Astra's (f,h) proportional to mine:", all(Fraction(a)==ratio*b for a,b in zip(astra,kv)), "ratio", ratio)
f = kv[:3]; h = kv[3:]
fp = f[0]+f[1]*z+f[2]*z**2; hp = h[0]+h[1]*z+h[2]*z**2
res = sympy.resultant(fp, hp, z)
print("resultant(f,h) =", res, "factored", factorint(abs(int(res))) if res != 0 else 0)

def choose_K(mode):
    if mode == 'small':
        rng = random.Random(31)
        base = 2*5*13*17*29*31*37*47
        from sympy.ntheory.modular import crt as _crt
        t0, md = _crt(s, [pow(base, -1, p) for p in s]); t0, md = int(t0), int(md)
        while True:
            K = base*(t0 + md*rng.randrange(1, 300))
            if any(K % p != 1 for p in s):
                continue
            J = {e: K*K+u[e[0]]*u[e[1]] for e in pairs}
            if any(gcd(J[a], J[b]) != 1 for a,b in combinations(pairs,2)):
                continue
            if any(gcd(J[e], abs(u[a]-u[b])) != 1 for e in pairs for a,b in pairs):
                continue
            fac = {e: factorint(J[e]) for e in pairs}
            if all(any(ex == 1 for ex in fac[e].values()) for e in pairs):
                return K, fac
    else:  # Astra's recipe: K even, divisible by primes of differences of products (not s), K = 1 mod s_i,
           # plus a distinguished split prime q_e with q_e || J_e
        prods = [u[i]*u[j] for i,j in pairs]
        assert len(set(prods)) == 10
        dp = set()
        for a,b in combinations(prods,2):
            dp |= set(factorint(abs(a-b)).keys())
        dp -= set(s)
        base = 2*prod(dp)
        qs = []
        q = max(prods)+1
        while len(qs) < 10:
            if isprime(q) and q % 4 == 1 and q not in dp: qs.append(q)
            q += 1
        # CRT: K = 0 mod base, K = 1 mod s_i, K = r_e mod q_e^2 with r_e^2 = -u_iu_j mod q_e, not mod q_e^2
        mods = [base] + s[:]; vals = [0] + [1]*5
        for e,q in zip(pairs, qs):
            c = (-u[e[0]]*u[e[1]]) % q
            r = sympy.sqrt_mod(c, q)
            # lift r to a residue mod q^2 with r^2 + u u' = 0 mod q but not mod q^2
            for t in range(q):
                R = r + t*q
                if (R*R + u[e[0]]*u[e[1]]) % (q*q) != 0: break
            mods.append(q*q); vals.append(R)
        from sympy.ntheory.modular import crt
        K, Mod = crt(mods, vals)
        K = int(K)
        if K < 10**30: K += int(Mod)
        J = {e: K*K+u[e[0]]*u[e[1]] for e in pairs}
        fac = {e: {q: (1 if J[e] % q == 0 and J[e] % (q*q) else None)} for e,q in zip(pairs, qs)}
        return K, fac

mode = sys.argv[1] if len(sys.argv) > 1 else 'small'
K, fac = choose_K(mode)
J = {e: K*K+u[e[0]]*u[e[1]] for e in pairs}
print("\nMODE", mode, " K =", K, " (digits %d)" % len(str(K)))
# block checks
assert K % 2 == 0 and all(K % p == 1 for p in s)
assert all(J[e] % 2 == 1 for e in pairs)
assert all(gcd(J[a],J[b]) == 1 for a,b in combinations(pairs,2)), "blocks not coprime"
assert all(gcd(K, s[i]*s[j]) == 1 for i,j in pairs)  # => J = K^2+(s_i s_j)^2 primitive sum of 2 squares => odd primes split
for e in pairs:
    for d in set(abs(u[a]-u[b]) for a,b in pairs):
        assert gcd(J[e], d) == 1
if mode == 'small':
    for e in pairs:
        assert all(p % 4 == 1 for p in fac[e]), "non-split prime"
    print("all ten blocks fully factored; all primes = 1 mod 4; each block has a prime to exponent one")
    dist = {e: min(p for p,ex in fac[e].items() if ex == 1) for e in pairs}
    print("squarefree blocks:", sum(all(x==1 for x in fac[e].values()) for e in pairs), "/10")
else:
    for e in pairs:
        assert list(fac[e].values())[0] == 1
    dist = {e: list(fac[e].keys())[0] for e in pairs}
    print("each block contains its distinguished prime to exponent exactly one")

N = prod(J.values()); mp.dps = len(str(N)) + 60
# --- profile after complementing the cycle cuts
marked = {}
for e in pairs:
    marked[e] = set(e) if e not in Fset else set(range(5)) - set(e)
n = [prod(J[e] for e in pairs if i in marked[e]) for i in range(5)]
g = {(i,j): prod(J[e] for e in pairs if i in marked[e] and j in marked[e]) for i,j in pairs}
for i in range(5): assert sum(1 for e in pairs if i in marked[e]) == 5
for p_ in pairs: assert sum(1 for e in pairs if set(p_) <= marked[e]) == 2
D = prod(J[e] for e in F)
A = [prod(J[e] for e in F if i in e) for i in range(5)]
for i,j in pairs:
    assert g[i,j] == D*J[i,j]//(A[i]*A[j]) and D*J[i,j] % (A[i]*A[j]) == 0

# --- index one: valuation columns at the distinguished primes (anchored at mark 0)
colmat = []
for e in pairs:
    col = [1 if i in marked[e] else 0 for i in range(5)]
    colmat.append([col[i]-col[0] for i in range(1,5)])
from sympy.matrices.normalforms import smith_normal_form
snf = smith_normal_form(Matrix(colmat).T, domain=sympy.ZZ)
print("Smith form of anchored 4x10 valuation matrix (distinguished primes only):", [snf[i,i] for i in range(4)])

# --- lattice L = {(a,b): a + c b = 0 mod D}, c = u_i+u_j mod J_ij on F
from sympy.ntheory.modular import crt
c, _ = crt([J[e] for e in F], [(u[e[0]]+u[e[1]]) % J[e] for e in F])
c = int(c)
def inL(a,b): return (a + c*b) % D == 0
fK = f[0]+f[1]*K**2+f[2]*K**4; hK = h[0]+h[1]*K**2+h[2]*K**4
assert inL(fK, hK)
alpha, beta = (fK + c*hK)//D, hK      # coordinates in basis (D,0), (-c,1)
gg = gcd(alpha, beta)
print("lattice content of (f(K^2),h(K^2)):", gg, "  plain gcd:", gcd(fK,hK))
alpha//=gg; beta//=gg
v = (alpha*D - c*beta, beta)
assert inL(*v)
# complete basis: alpha*delta - beta*gamma = 1
x1, y1, gq = sympy.gcdex(alpha, -beta)
gq, x1, y1 = int(gq), int(x1), int(y1)   # alpha*x1 + (-beta)*y1 = gq
assert gq in (1,-1)
delta, gamma = x1*gq, y1*gq
assert alpha*delta - beta*gamma == 1
w = (gamma*D - c*delta, delta)
assert inL(*w)
assert abs(w[0]*v[1]-w[1]*v[0]) == D
# reduce w against v
vv = v[0]**2+v[1]**2
t = round(Fraction(w[0]*v[0]+w[1]*v[1], vv))
w = (w[0]-t*v[0], w[1]-t*v[1])
print("log_K |v| = %.4f, log_K |v_2| = %.4f, log_K |w| = %.4f" % (
    float(mlog(mpf(vv))/2/mlog(K)), float(mlog(abs(mpf(v[1])))/mlog(K)), float(mlog(mpf(w[0]**2+w[1]**2))/2/mlog(K))))
j0 = (K**7 - w[1])//v[1]
U = [(u[i], u[i]**2 - K*K) for i in range(5)]
for i,j in pairs:
    assert U[i][0]*U[j][1]-U[i][1]*U[j][0] == J[i,j]*(u[j]-u[i])
    assert gcd(*U[i]) == 1
chosen = None
for jj in range(j0-1, j0+3):
    r1 = (w[0]+jj*v[0], w[1]+jj*v[1])
    Mat = (r1, v)
    assert inL(*r1)
    V = []
    for i in range(5):
        a = r1[0]*U[i][0]+r1[1]*U[i][1]; b = v[0]*U[i][0]+v[1]*U[i][1]
        assert a % A[i] == 0 and b % A[i] == 0
        V.append((a//A[i], b//A[i]))
    norms = [x*x+y*y for x,y in V]
    eq = [norms[i] == n[i] for i in range(5)]
    print(" j=j0%+d: norms equal to n_i? %s ; second entry of row1 / K^7 = %.6f" % (jj-j0, eq, float(Fraction(r1[1], K**7))))
    if chosen is None and not all(eq):
        chosen = (jj, r1, V, norms)
jj, r1, V, norms = chosen
print("\nUsing j = j0%+d" % (jj-j0))
# primitivity
print("columns primitive:", all(gcd(x,y) == 1 for x,y in V))
# divisibility & cofactors
cof = {}
for i,j in pairs:
    dt = V[i][0]*V[j][1]-V[i][1]*V[j][0]
    assert dt % g[i,j] == 0
    cof[i,j] = dt//g[i,j]
print("all ten g_ij | det(V_i,V_j); cofactors:", {k: cof[k] for k in pairs})
print("cofactor == +-(u_j-u_i) for all pairs:", all(abs(cof[i,j]) == abs(u[j]-u[i]) for i,j in pairs))
# norms
for i in range(5):
    rel = Fraction(norms[i]-n[i], n[i])
    rel2 = Fraction(norms[i], isqrt(N))  # vs sqrt N
    print(" i=%d |V|^2/n_i - 1 = %.3e   K^2*(|V|^2/n_i-1) = %.4f   |V|^4/N - 1 = %.3e  n_i | |V|^2 ? %s" % (
        i, float(rel), float(rel*K*K), float(Fraction(norms[i]**2, N)-1), norms[i] % n[i] == 0))
# angles
ang = [2*atan2(mpf(y), mpf(x)) for x,y in V]
diam = max(ang)-min(ang)
print("diam of doubled args = %s" % mp.nstr(diam, 8))
print("diam * N^(3/10) = %s ;  diam * N^(1/4) = %s  (shortened convention needs <= 1/2)" % (
    mp.nstr(diam*mpf(N)**(mpf(3)/10), 8), mp.nstr(diam*mpf(N)**(mpf(1)/4), 8)))
print("min pair doubled-angle sep * N^(3/10) = %s" % mp.nstr(min(abs(ang[i]-ang[j]) for i,j in pairs)*mpf(N)**(mpf(3)/10), 8))
print("V_0 =", V[0])
# Gaussian isotropy at distinguished primes: is V_i mod p on an isotropic line?
iso = 0; tot = 0
for e in pairs:
    p = dist[e]
    for i in marked[e]:
        tot += 1
        if (V[i][0]**2+V[i][1]**2) % p == 0: iso += 1
print("isotropic reductions at distinguished primes: %d of %d marked (column,prime) incidences" % (iso, tot))
print("projective invariance: cross-ratios of V equal those of U:",
      all(Fraction((V[a][0]*V[b][1]-V[a][1]*V[b][0])*(V[cc][0]*V[d][1]-V[cc][1]*V[d][0]),
                   (V[a][0]*V[d][1]-V[a][1]*V[d][0])*(V[cc][0]*V[b][1]-V[cc][1]*V[b][0])) ==
          Fraction((U[a][0]*U[b][1]-U[a][1]*U[b][0])*(U[cc][0]*U[d][1]-U[cc][1]*U[d][0]),
                   (U[a][0]*U[d][1]-U[a][1]*U[d][0])*(U[cc][0]*U[b][1]-U[cc][1]*U[b][0]))
          for a,b,cc,d in [(0,1,2,3),(0,1,2,4),(1,2,3,4),(0,2,3,4)]))

# --- (L_5) datum: anchor column 0 at infinity, P_j = (x_j, 1) with x_j = common/(t_j - t_0)
# t_i = U_i[1]/U_i[0]
tt = [Fraction(U[i][1], U[i][0]) for i in range(5)]
xs = [1/(tt[i]-tt[0]) for i in range(1,5)]
den = 1
for x in xs: den = den*x.denominator//gcd(den, x.denominator)
X = [int(x*den) for x in xs]
gx = 0
for i in range(4):
    for j in range(i+1,4): gx = gcd(gx, X[i]-X[j])
cont = 0
for x in X: cont = gcd(cont, x)
print("\n(L_5) chart: x_j (j=1..4) integers; content", cont)
# blocks in anchored L-form: T subset of {1..4}: pair cut {a,b} -> T = {a,b} if 0 not in it else complement in {1..4}
H = 0
for a in range(1,5):
    for b in range(a+1,5):
        blk = 1
        for e in pairs:
            T = set(e) if 0 not in e else set(range(1,5)) - set(e)
            if a in T and b in T: blk *= J[e]
        dif = X[a-1]-X[b-1]
        assert dif % blk == 0, (a,b)
        H += abs(dif//blk)
        print("  t_%d%d = %s" % (a, b, dif//blk))
print("  (with P_0=(1,0), P_j=(x_j,1): Y_j = 1; cofactors above; chain condition trivial: blocks pairwise coprime)")
