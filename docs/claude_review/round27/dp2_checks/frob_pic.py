"""Characteristic polynomial of Frob_p on NS(S_{Fbar_p}) (rank 8) from point counts over F_{p^n}, n=1..nmax.
chi(x) = (x-1) psi(x) (K is Frobenius-invariant); psi recovered from power sums s_n - 1 (n=1..7) by Newton.
If psi(1) != 0 then rank NS(S)^{Gal(Qbar/Q)} = 1 (the Galois-invariant part is contained in the Frob-invariant part),
so every curve on S defined over Q has class -nK."""
import sys, time
import sympy as sp
from ffcount import count_S
from search_f import candidate
from geom import default_f
x, y, z, X = sp.symbols('x y z X')

def charpoly_from_traces(tr):
    # power sums of psi's roots
    s = [None] + [tr[n] - 1 for n in range(1, 8)]
    e = [1] + [0]*7
    for k in range(1, 8):
        acc = 0
        for i in range(1, k + 1):
            acc += (-1)**(i - 1) * e[k - i] * s[i]
        e[k] = sp.Rational(acc, k)
    psi = sum((-1)**k * e[k] * X**(7 - k) for k in range(8))
    return sp.expand(psi)

if __name__ == '__main__':
    p = int(sys.argv[1]); nmax = int(sys.argv[2]); seeds = [int(a) for a in sys.argv[3:]]
    for seed in seeds:
        f = candidate(seed) if seed >= 0 else default_f()
        monos = [((i, j, k), int(co)) for (i, j, k), co in sp.Poly(f, x, y, z).terms()]
        tr = {}
        for n in range(1, nmax + 1):
            t0 = time.time()
            N = count_S(monos, p, n)
            q = p**n
            assert (N - q*q - 1) % q == 0
            tr[n] = (N - q*q - 1)//q
        psi = charpoly_from_traces(tr)
        fac = sp.factor_list(psi)
        print('seed', seed, 'traces', [tr[n] for n in range(1, nmax+1)], '| psi =', sp.factor(psi), '| psi(1) =', psi.subs(X, 1), flush=True)
        if nmax >= 8:
            # consistency check: power sum s_8 from psi roots
            roots = sp.Poly(psi, X).nroots()
            s8 = 1 + sum(r**8 for r in roots)
            print('   check n=8: predicted trace', sp.N(s8, 8), 'counted', tr[8])
