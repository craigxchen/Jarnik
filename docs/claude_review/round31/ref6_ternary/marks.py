"""Raw marks F_j(t, T(t)) in F_9[t] for V = i*F_q^k (coords on 1,t,..,t^{k-1}), via subspace polynomials.
Own arithmetic (f9lib). q = 3^r, r odd."""
import numpy as np, itertools
from f9lib import F9P, pw

def frobq(y, q):
    """y(t)^q for y in F_9[t], q = 3^r with r odd: coefficients conjugated, t -> t^q"""
    n = len(y.re)
    re = np.zeros((n-1)*q + 1 if n else 0, np.int64); im = re.copy()
    re[::q] = y.re; im[::q] = -y.im       # a^q = conj(a) on F_9 since r odd
    return F9P(re, im)

def kernel_basis(lam):
    k = len(lam)
    vecs = [x for x in itertools.product(range(3), repeat=k) if sum(a*b for a, b in zip(lam, x)) % 3 == 0 and any(x)]
    basis = []
    span = {tuple([0]*k)}
    for x in vecs:
        if x not in span:
            basis.append(x)
            span = {tuple((a + c*b) % 3 for a, b in zip(s, x)) for s in span for c in range(3)}
        if len(basis) == k-1: break
    return basis

def ipoly(x):
    return F9P([0]*len(x), list(x))      # i * sum x_b t^b

def raw_mark(lam, eps, q, T):
    """F = A (A^((q-1)/2) - c^((q-1)/2))^2 with A = P_H(T), c = P_H(v0), eps*lam(v0) = 1"""
    k = len(lam)
    piv = [b for b in range(k) if lam[b]][0]
    v0 = [0]*k; v0[piv] = (eps * lam[piv]) % 3            # lam(v0) = lam[piv]^2*eps = eps; eps*lam(v0)=1
    vals = [T, ipoly(v0)] + [ipoly(u) for u in kernel_basis(lam)]
    # build P_H by adjoining basis vectors u of H one at a time: P <- P^q - P(u)^(q-1) P
    nb = len(vals) - 2
    for s in range(nb):
        u = vals[2 + s]
        uq1 = pw(u, q-1) if q > 1 else None
        vals = [frobq(y, q) - uq1 * y for y in vals]
    A, c = vals[0], vals[1]
    hq = (q-1)//2
    S = pw(A, hq) - pw(c, hq)
    return A * S * S, A, c

def check_family(forms, signs, q, L, Textra=None, gcd=False):
    t = F9P([0, 1])
    T = F9P([0]*L + [1]) if Textra is None else F9P([0]*L + [1]) + Textra
    marks = [raw_mark(l, e, q, T)[0] for l, e in zip(forms, signs)]
    return marks
