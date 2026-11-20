"""Function-field check of the projection argument on the five-mark ternary witness (D=12).
Referee's own code; uses python-flint only as an F_9[t] arithmetic library. The five encoded marks
are copied as data from the Part VII witness list quoted in Astra's script."""
from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx
from itertools import combinations, product
P3 = fmpz_mod_poly_ctx(3)
F9 = fq_default_ctx(modulus=P3([1, 0, 1]), var="i")
R9 = fq_default_poly_ctx(F9)
I = F9.gen()
def decode(index):
    cs = []
    while index:
        index, d = divmod(index, 9)
        hi, lo = divmod(d, 3)
        cs.append(F9(lo) + I*hi)
    return R9(cs)
def bar(z):
    return R9([c**3 for c in z.coeffs()])
codes = [391358742215, 391358571511, 391360438826, 391357359206, 391358930855]
Z = [decode(c) for c in codes]
N = Z[0]*bar(Z[0])
D = Z[0].degree()
print("degrees", [z.degree() for z in Z], "deg N", N.degree(), "monic", all(z.coeffs()[-1] == F9(1) for z in Z))
assert all(z*bar(z) == N for z in Z)
print("pair difference degrees", sorted((z-w).degree() for z, w in combinations(Z, 2)), "D/2 =", D/2)
# all elements of F9 and of F9[t]/lambda
F9el = [F9(a) + I*b for a in range(3) for b in range(3)]
def is_in_F3(p):
    return all(c**3 == c for c in p.coeffs())
hist = {}
for zi, zj in combinations(Z, 2):
    prod_ = zi*zj
    G = prod_.gcd(N)
    b0 = prod_ // G
    # find a constant eta with b = eta*b0 satisfying prod*bar(b) == N*b
    sols = [eta for eta in F9el if eta != F9(0) and prod_*bar(R9([eta])*b0) == N*(R9([eta])*b0)]
    assert sols, "no b"
    b = R9([sols[0]])*b0
    traces = [z*bar(b) + bar(z)*b for z in Z]
    assert all(is_in_F3(t) for t in traces)
    maxdeg = max((t-u).degree() for t, u in combinations(traces, 2) if t != u)
    ndist = len(set(str(t) for t in traces))
    T = b*bar(b)
    hist.setdefault((b.degree(), maxdeg, ndist), 0)
    hist[(b.degree(), maxdeg, ndist)] += 1
    assert maxdeg <= b.degree()
print("(deg b, max deg trace difference, #distinct traces): count", hist)
# general sum-two certificates c in {-1,..,2}^5: b/bbar = prod z^c / N; check deg of trace differences <= deg b
cnt = 0; worst = 0
for c in product(range(-1, 3), repeat=5):
    if sum(c) != 2: continue
    num = R9([1]); den = N
    for ci, z in zip(c, Z):
        if ci > 0: num = num*z**ci
        elif ci < 0: num = num*bar(z)**(-ci); den = den*N**(-ci)
    g = num.gcd(den)
    a_, d_ = num//g, den//g
    # a_/d_ = b/bbar in lowest terms: b ~ a_, bbar ~ d_ up to constants
    sols = [eta for eta in F9el if eta != F9(0) and num*bar(R9([eta])*a_) == den*(R9([eta])*a_)]
    if not sols: print("no b for", c); continue
    b = R9([sols[0]])*a_
    traces = [z*bar(b) + bar(z)*b for z in Z]
    md = max([(t-u).degree() for t, u in combinations(traces, 2) if t != u] + [-1])
    assert md <= b.degree(), (c, md, b.degree())
    cnt += 1; worst = max(worst, md - b.degree())
print("sum-two certificates with entries in [-1,2]:", cnt, "all satisfy deg(x_i-x_j) <= deg b; max(maxdeg - deg b) =", worst)
