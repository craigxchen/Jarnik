"""Theorem A checks (exact): the regular model of the conic u v = N t^2 over Z[i].

For each normalised cluster:
 (A1) at pi | p | N (p^e || N): v_pi(z_j) + v_pibar(z_j) = e; the section P_j=[z_j:zbar_j:1] reduces to
      component C_{a} of the A_{e-1} chain with a = v_pi(u) = v_pi(z_j); O_u=[1:0:0] -> C_0, O_v=[0:1:0] -> C_e.
 (A2) model local intersection i_pi(P_i,P_j) := 0 if a_i != a_j, else v_pi((z_i - z_j)/pi^a) (chart x = u/pi^a);
      at Gaussian primes varpi not dividing N: v_varpi(z_i - z_j).  Check prod Norm(varpi)^{i} == Norm(c_ij),
      c_ij = (z_i - z_j)/gcd(z_i,z_j), computed independently by the Gaussian Euclidean algorithm.
 (A3) vertical correction: solve A x = -n on the chain (A = intersection matrix), -Phi^2 = -x^T A x;
      check == sum_tau S_tau(c)^2 for random zero-sum charges, and == cut slack Q_p for D0 = sum P_j - M/2 (O_u+O_v).
 (A4) the integer identity prod Norm(z_i-z_j)^4 = N^{M(M-2)} prod_p p^{4 Q_p} prod Norm(c_ij)^4 (all exact),
      with Q_p computed from the chain correction (A3), not from the cut formula.
"""
import random
from fractions import Fraction
from itertools import combinations
from sympy import factorint, Matrix
from gauss import sub, conj, norm, gcd_g, exact_div, divides, valuation, split_prime_factor
from clusters import standard_clusters

def gaussian_primes_over(q):
    if q == 2:
        return [((1, 1), 2)]
    if q % 4 == 3:
        return [((q, 0), q * q)]
    pi = split_prime_factor(q)
    return [(pi, q), (conj(pi), q)]

def gpow(z, k):
    out = (1, 0)
    for _ in range(k):
        out = (out[0]*z[0]-out[1]*z[1], out[0]*z[1]+out[1]*z[0])
    return out

def chain_matrix(e):
    n = e + 1
    A = [[0]*n for _ in range(n)]
    for k in range(n):
        A[k][k] = -1 if (k == 0 or k == e) else -2
        if k + 1 < n:
            A[k][k+1] = A[k+1][k] = 1
    if e == 0:
        A[0][0] = 0
    return A

def vertical_energy(e, charges_by_level):
    """-Phi^2 for the chain C_0..C_e with charges n_k (sum 0), by exact linear algebra."""
    n = [Fraction(x) for x in charges_by_level]
    assert sum(n) == 0
    A = chain_matrix(e)
    # solve A x = -n with x_0 = 0 (A has kernel (1,..,1)); use rows 0..e-1, unknowns x_1..x_e
    E = e
    if E == 0:
        return Fraction(0)
    Msub = Matrix([[Fraction(A[r][c]) for c in range(1, E+1)] for r in range(0, E)])
    rhs = Matrix([-n[r] for r in range(0, E)])
    sol = Msub.LUsolve(rhs)
    x = [Fraction(0)] + [Fraction(str(s)) for s in sol]
    # check the last row too
    last = sum(A[E][c]*x[c] for c in range(E+1))
    assert last == -n[E]
    Phi2 = sum(x[r]*A[r][c]*x[c] for r in range(E+1) for c in range(E+1))
    return -Phi2

def check_cluster(cl, rng, log):
    pts, N, M = cl['pts'], cl['N'], cl['M']
    fN = factorint(N)
    assert N % 2 == 1 and all(p % 4 == 1 for p in fN), 'not normalised'
    levels = {}
    ok = True
    # (A1)
    for p, e in fN.items():
        pi = split_prime_factor(p)
        pib = conj(pi)
        a = []
        for z in pts:
            ap = valuation(pi, z); bp = valuation(pib, z)
            ok &= (ap + bp == e)
            a.append(ap)
        ok &= (min(a) == 0 and max(a) == e)
        levels[p] = (pi, e, a)
    # (A2)
    nA2 = 0
    for i, j in combinations(range(M), 2):
        d = sub(pts[i], pts[j])
        c = exact_div(d, gcd_g(pts[i], pts[j]))
        prod = 1
        for q in factorint(norm(d)):
            for varpi, Nv in gaussian_primes_over(q):
                if q in fN:
                    pi, e, a = levels[q]
                    # place pi or pibar; levels at pibar are e - a
                    ai = valuation(varpi, pts[i]); aj = valuation(varpi, pts[j])
                    if ai != aj:
                        inter = 0
                    else:
                        x = exact_div(d, gpow(varpi, ai))
                        inter = valuation(varpi, x)
                else:
                    inter = valuation(varpi, d)
                prod *= Nv ** inter
        ok &= (prod == norm(c))
        nA2 += 1
    # (A3) random charges and D0
    Qp = {}
    for p, (pi, e, a) in levels.items():
        for _ in range(5):
            cs = [rng.randint(-3, 3) for _ in range(M)]
            cu = rng.randint(-3, 3)
            cv = -sum(cs) - cu
            nk = [0]*(e+1)
            for j in range(M):
                nk[a[j]] += cs[j]
            nk[0] += cu; nk[e] += cv
            S = [sum(nk[:t]) for t in range(1, e+1)]
            ok &= (vertical_energy(e, nk) == sum(s*s for s in S))
        # D0 (use 2*D0 to keep integers): charges 2 on P_j, -M on O_u and O_v
        nk = [0]*(e+1)
        for j in range(M):
            nk[a[j]] += 2
        nk[0] -= M; nk[e] -= M
        ve = vertical_energy(e, nk)          # = 4 * Q_p
        Scut = [sum(1 for j in range(M) if a[j] < t) for t in range(1, e+1)]
        ok &= (ve == sum((2*s - M)**2 for s in Scut))
        Qp[p] = ve                            # this is sum (2S-M)^2 = 4 Q_p
    # (A4) exact integer identity
    lhs = 1
    rhs_c = 1
    for i, j in combinations(range(M), 2):
        d = sub(pts[i], pts[j])
        lhs *= norm(d) ** 4
        c = exact_div(d, gcd_g(pts[i], pts[j]))
        rhs_c *= norm(c) ** 4
    rhs = N ** (M*(M-2)) * rhs_c
    for p, v in Qp.items():
        assert v.denominator == 1
        rhs *= p ** int(v)
    ok &= (lhs == rhs)
    log.write('%-24s M=%2d N=%s  pairs=%d  A1-A4 %s\n' % (cl['name'], M, dict(fN), nA2, 'PASS' if ok else 'FAIL'))
    return ok

if __name__ == '__main__':
    rng = random.Random(20)
    allok = True
    with open('dictionary_out.txt', 'w') as log:
        for cl in standard_clusters():
            allok &= check_cluster(cl, rng, log)
            log.flush()
        log.write('ALL PASS\n' if allok else 'SOME FAIL\n')
    print(open('dictionary_out.txt').read())
