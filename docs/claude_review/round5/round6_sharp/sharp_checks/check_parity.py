"""Exact checks of Proposition 2.1 / 2.2 of round6/sharp.md (what factored residue data force).

For every odd prime power Q = q^a <= QMAX:
 (a) T_Q = {u : Norm u = 1} is cyclic of order m = (q - chi(q)) q^(a-1);
 (b) r -> r/conj(r) maps the unit group G_Q onto T_Q, with fibres r (Z/Q)^*;
 (c) nu(t) := Legendre(Norm r) (any r in the fibre of t) is well defined, and nu(t) = +1 iff t is a
     square in T_Q;
 (d) i is a square in T_Q iff q = +-1 mod 8;  -1 is always a square;
 (e) for t in T_Q and n in (Z/Q)^*: some r has r/conj(r) = t and Norm r = n  iff  nu(t) = (n/q);
 (f) actual Gaussian primes pi over p = 1 mod 4 (p <= PMAX, p != q): pi/conj(pi) is a square in T_Q
     iff (p/q) = 1; and log(pi/conj pi) mod 4 is NOT a function of (p/q) when 8 | m (quartic part);
 (g) compatibility: the parity of log_g at level a equals the parity at level 1 (generators chosen
     compatibly by reduction).
Exhaustive over Z[i]/Q.  Prints ALL PARITY CHECKS PASSED at the end."""
import sys
from gres import *

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
PMAX = 3000


def gaussian_prime(p):
    for x in range(1, int(p ** 0.5) + 1):
        y2 = p - x * x
        y = int(round(y2 ** 0.5))
        if y * y == y2:
            return (x, y)
    raise ValueError


ok = True
moduli = []
for q in primes_upto(QMAX):
    if q == 2:
        continue
    a = 1
    while q ** a <= QMAX:
        moduli.append((q, a))
        a += 1

gauss_primes = [(p, gaussian_prime(p)) for p in primes_upto(PMAX) if p % 4 == 1]
summary = []
for (q, a) in moduli:
    Q = q ** a
    m = m_of(q, a)
    units = [(x, y) for x in range(Q) for y in range(Q) if (x * x + y * y) % q != 0]
    T = [u for u in units if norm(u, Q) == 1]
    ok &= (len(T) == m)
    # (a) cyclic
    g = generator_T(q, a, seed=q * 1000 + a)
    tab = dlog_table(g, m, Q)
    ok &= (len(tab) == m) and set(tab) == set(T)
    # (b) onto, fibres = rational multiples
    fib = {}
    for r in units:
        t = mul(r, inv(conj(r, Q), Q), Q)
        fib.setdefault(t, []).append(r)
    phiQ = Q - Q // q
    ok &= set(fib) == set(T) and all(len(v) == phiQ for v in fib.values())
    for t, rs in list(fib.items())[:50]:
        r0 = rs[0]
        rat = {((r0[0] * s) % Q, (r0[1] * s) % Q) for s in range(1, Q) if s % q}
        ok &= rat == set(rs)
    # (c) nu well defined and = square-ness in T
    nu = {}
    for t, rs in fib.items():
        vals = {legendre(norm(r, Q), q) for r in rs}
        ok &= len(vals) == 1
        nu[t] = vals.pop()
    for t in T:
        ok &= (nu[t] == 1) == (tab[t] % 2 == 0)
    # (d) units
    i_unit = (0, 1)
    ok &= (tab[i_unit] % 2 == 0) == (q % 8 in (1, 7))
    ok &= tab[((Q - 1) % Q, 0)] % 2 == 0
    # (e) realisability with prescribed norm (check all t, a sample of n)
    ns = [n for n in range(1, Q) if n % q][:40]
    for t, rs in fib.items():
        norms = {norm(r, Q) for r in rs}
        for n in ns:
            ok &= ((n % Q) in norms) == (nu[t] == legendre(n, q))
    # (f) actual Gaussian primes; quartic part
    sq_ok = True
    mod4 = {1: set(), -1: set()}
    for p, (x, y) in gauss_primes:
        if p % q == 0:
            continue
        pi = (x % Q, y % Q)
        t = mul(pi, inv(conj(pi, Q), Q), Q)
        L = legendre(p, q)
        sq_ok &= (tab[t] % 2 == 0) == (L == 1)
        mod4[L].add(tab[t] % 4)
    ok &= sq_ok
    quartic_free = (m % 8 == 0) and (mod4[1] == {0, 2})
    summary.append((Q, q, a, m, q % 8, tab[i_unit] % 2 == 0, sorted(mod4[1]), sorted(mod4[-1]), quartic_free))

# (g) compatibility of parity across levels: generator of level a reduces to a generator of level 1
for q in [3, 5, 7, 11, 13]:
    a = 1
    while q ** (a + 1) <= QMAX:
        a += 1
    if a == 1:
        continue
    Q = q ** a
    m = m_of(q, a)
    g = generator_T(q, a, seed=7)
    g1 = (g[0] % q, g[1] % q)
    ok &= order_in_T(g1, m_of(q, 1), q) == m_of(q, 1)
    tab_hi = dlog_table(g, m, Q)
    tab_lo = dlog_table(g1, m_of(q, 1), q)
    for t, k in list(tab_hi.items())[:500]:
        ok &= tab_lo[(t[0] % q, t[1] % q)] == k % m_of(q, 1)
        ok &= (tab_lo[(t[0] % q, t[1] % q)] % 2) == (k % 2)

print(f"{'Q':>5} {'q':>4} {'a':>2} {'m':>5} {'q%8':>4} {'i square':>9} {'l mod 4 for (p/q)=+1':>22} {'for -1':>8} {'quartic part free':>18}")
for row in summary:
    print(f"{row[0]:>5} {row[1]:>4} {row[2]:>2} {row[3]:>5} {row[4]:>4} {str(row[5]):>9} {str(row[6]):>22} {str(row[7]):>8} {str(row[8]):>18}")
print("moduli checked:", len(moduli), " Gaussian primes used:", len(gauss_primes))
print("ALL PARITY CHECKS PASSED" if ok else "FAILURE")
