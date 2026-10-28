"""Independent check of Theorem A (1)-(4) of arakelov.md."""
import random
from fractions import Fraction
from sympy import factorint, Matrix
from kgauss import *
from kclusters import configs, normalise

def chain_minus_phi2(e, charges):
    """charges: list n_0..n_e (sum 0). Solve (D+Phi).C_k = 0 on chain C_0..C_e
    (C_0^2=C_e^2=-1, middle -2, adjacent 1). Return -Phi^2 (exact)."""
    n = e + 1
    A = [[0]*n for _ in range(n)]
    for k in range(n):
        A[k][k] = -1 if k in (0, e) else -2
        if k+1 < n:
            A[k][k+1] = A[k+1][k] = 1
    if e == 0:
        return Fraction(0)
    # fix x_0 = 0, solve rows 1..e of A x = -n (rows 0..e-1 suffice since sum is consistent)
    # unknowns x_1..x_e; use equations k = 0..e-1
    Ared = Matrix([[A[k][j] for j in range(1, n)] for k in range(0, e)])
    b = Matrix([-charges[k] for k in range(0, e)])
    x = Ared.LUsolve(b)
    xs = [0] + [Fraction(int(v.p), int(v.q)) for v in x]
    # check all equations
    for k in range(n):
        assert sum(A[k][j]*xs[j] for j in range(n)) == -charges[k]
    phi2 = sum(xs[i]*A[i][j]*xs[j] for i in range(n) for j in range(n))
    return -phi2

def run():
    rng = random.Random(7)
    lines = []
    allok = True
    for name, n0, pts, C, isarc in configs():
        q, g = normalise(pts)
        M = len(q)
        N = gnorm(q[0])
        assert all(gnorm(z) == N for z in q)
        fac = factorint(N)
        assert all(p % 4 == 1 for p in fac), (name, fac)
        ok = True
        Q4 = 0  # 4*sum Q_p log p as exponent dict
        Qp = {}
        prof = {}
        for p, e in fac.items():
            pi = gprime_above(p); pib = gconj(pi)
            a = [gval(z, pi) for z in q]
            ab = [gval(z, pib) for z in q]
            assert all(a[j] + ab[j] == e for j in range(M))
            assert min(a) == 0 and max(a) == e  # normalisation (gcd 1)
            prof[p] = a
            # A(1): component check via chart coordinates
            for j in range(M):
                aj = a[j]
                if aj < e:
                    alpha = gdiv_exact(q[j], pi_pow(pi, aj))
                    assert gval(alpha, pi) == 0
                    beta_val = gval(gconj(q[j]), pi) - (e - aj - 1)
                    assert beta_val >= 1  # beta = v/pi^(e-a-1) vanishes mod pi -> on C_a
                else:
                    assert gval(gconj(q[j]), pi) == 0  # beta = v unit, alpha = u/pi^(e-1) vanishes -> on C_e
            # A(2): local intersection via chart coordinate vs v_pi(c_ij)
            for i in range(M):
                for j in range(i+1, M):
                    d = gsub(q[i], q[j]); gij = ggcd(q[i], q[j]); cij = gdiv_exact(d, gij)
                    for P, lev, other in ((pi, a, None), (pib, ab, None)):
                        if lev[i] != lev[j]:
                            ichart = 0
                        else:
                            L = lev[i]
                            if L < e:
                                ichart = gval(gsub(gdiv_exact(q[i], pi_pow(P, L)), gdiv_exact(q[j], pi_pow(P, L))), P)
                            else:
                                ichart = gval(gsub(gconj(q[i]), gconj(q[j])), P)
                        if ichart != gval(cij, P):
                            ok = False
            # A(3): vertical correction for random zero-sum charges and for 2*D0
            for trial in range(5):
                c = [rng.randint(-3, 3) for _ in range(M)]
                s = sum(c); c[0] -= s
                ch = [0]*(e+1)
                for j in range(M):
                    ch[a[j]] += c[j]
                lhs = chain_minus_phi2(e, ch)
                rhs = sum(sum(c[j] for j in range(M) if a[j] < t)**2 for t in range(1, e+1))
                if lhs != rhs:
                    ok = False
            ch = [0]*(e+1)
            for j in range(M):
                ch[a[j]] += 2
            ch[0] -= M; ch[e] -= M
            lhs = chain_minus_phi2(e, ch)
            S = [sum(1 for j in range(M) if a[j] < t) for t in range(1, e+1)]
            Q4p = sum((2*s - M)**2 for s in S)  # = 4 Q_p
            if lhs != Q4p:
                ok = False
            Qp[p] = Q4p
        # A(4): exact identity prod Norm(z_i-z_j)^4 = N^{M(M-2)} prod p^{4Q_p} prod Norm(c_ij)^4
        L = 1; Rr = N**(M*(M-2))
        for p, v in Qp.items():
            Rr *= p**v
        for i in range(M):
            for j in range(i+1, M):
                d = gsub(q[i], q[j]); gij = ggcd(q[i], q[j]); cij = gdiv_exact(d, gij)
                L *= gnorm(d)**4
                Rr *= gnorm(cij)**4
                # g_ij supported on primes of N
                ng = gnorm(gij)
                for pp in factorint(ng):
                    if pp not in fac:
                        ok = False
        ok = ok and (L == Rr)
        allok &= ok
        lines.append('%-24s M=%2d N=%s ok=%s' % (name, M, fac, ok))
    lines.append('ALL OK: %s' % allok)
    return lines

def pi_pow(pi, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, pi)
    return r

if __name__ == '__main__':
    out = run()
    open('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round20/kill_arakelov_checks/kdict_out.txt', 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
