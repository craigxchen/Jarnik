# Referee's independent exact library (no author code imported).
import itertools, flint
def cube_points(d):
    return [p for p in itertools.product((-1,0,1), repeat=d)]
def slabbox(d, cuts):
    """cuts: list of (S, L) with S subset of range(d); |sum_{x in S} k_x| <= L"""
    return [k for k in cube_points(d) if all(abs(sum(k[x] for x in S)) <= L for S, L in cuts)]
def monomials(d, maxdeg):
    mons = [e for e in itertools.product((0,1,2), repeat=d) if sum(e) <= maxdeg]
    mons.sort(key=sum); return mons
def evalmon(k, e):
    v = 1
    for a, b in zip(k, e):
        if b: v *= a**b
    return v
class Exact:
    """exact V(A,i) bases over Z for i = 1..imax; nu(u) for integer u."""
    def __init__(self, A, imax=None):
        self.A = A; d = len(A[0]); self.d = d
        mons = monomials(d, 2*d)
        self.H = []
        self.N = {}
        n = len(A)
        i = 1
        while True:
            cols = [e for e in mons if sum(e) < i]
            M = flint.fmpz_mat([[evalmon(k, e) for e in cols] for k in A])  # n x #cols
            r = M.rank()
            self.H.append(r)
            if r == n:
                self.rho = i - 1; break
            # left nullspace: c with c^T M = 0  <=> M^T c = 0
            X, nul = M.transpose().nullspace()
            self.N[i] = [[X[a, j] for a in range(n)] for j in range(nul)]
            i += 1
            if imax and i > imax: self.rho = None; break
    def nu(self, u):
        best = 0
        for i, basis in self.N.items():
            lv = [sum(a*b for a, b in zip(k, u))**i for k in self.A]
            if any(sum(c*l for c, l in zip(cv, lv)) != 0 for cv in basis):
                best = max(best, i)
        return best
