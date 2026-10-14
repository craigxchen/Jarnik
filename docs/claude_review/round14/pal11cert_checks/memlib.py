"""Membership of a degree-8 form (given by coefficients on reduced monomials) in T_8 = J(A)_8^perp for pal11, mod p,
in the Z_5-isotypic coordinates of isotan.py.  A form F = sum_a F_a u^a lies in span{Phi_c : c in V(A,8)} iff the
vector vec(a) = F_a / multinom(8;a) lies in T_8 = {mu_8(c)} (fplib docstring).  Mod p non-membership is a rigorous
rational non-membership because rank_p K_{<=8} = rank_Q K_{<=8} = |A| (isotan.py; referee_fatpoint §2.5 criterion)."""
import pickle
from math import factorial
from flint import nmod_mat
class T8:
    def __init__(self, fn):
        r = pickle.load(open(fn, 'rb'))
        self.p = r['p']; self.w = r['w']; self.perm = r['perm']
        self.Mre = r['Mre8']; self.Mfix = r['Mfix8']; self.T = r['T']
        p = self.p
        self.chi = [[pow(self.w, (j * t) % 5, p) for t in range(5)] for j in range(5)]
        self.Tm = [nmod_mat(len(Tj), len(Tj[0]), [x for row in Tj for x in row], p) for Tj in self.T]
        self.rk = [len(Tj) for Tj in self.T]
    def coords(self, vec):
        """vec: dict reduced-exponent-tuple -> int mod p (missing = 0). Returns per-character coordinate lists
        (all scaled by 5: free rep a -> sum_t chi_j(g^t) vec(g^t a); fixed a (j = 0 only) -> 5 vec(a))."""
        p = self.p
        out = []
        for j in range(5):
            c = []
            for orb in self.Mre:
                s = 0
                for t in range(5):
                    v = vec.get(orb[t], 0)
                    if v: s += self.chi[j][t] * v
                c.append(s % p)
            if j == 0:
                c += [(5 * vec.get(a, 0)) % p for a in self.Mfix]
            out.append(c)
        return out
    def member(self, vec):
        cs = self.coords(vec)
        for j in range(5):
            if not any(cs[j]): continue
            Mj = nmod_mat(self.rk[j] + 1, len(cs[j]), [x for row in self.T[j] for x in row] + cs[j], self.p)
            if Mj.rank() > self.rk[j]: return False
        return True
def multinom(m, a):
    v = factorial(m)
    for ax in a: v //= factorial(ax)
    return v
def product_coeffs(Z):
    d = len(Z[0]); poly = {tuple([0] * d): 1}
    for z in Z:
        new = {}
        nz = [(x, zx) for x, zx in enumerate(z) if zx]
        for e, c in poly.items():
            for x, zx in nz:
                e2 = list(e); e2[x] += 1; e2 = tuple(e2)
                new[e2] = new.get(e2, 0) + c * zx
        poly = {e: c for e, c in new.items() if c}
    return poly
def product_vec(Z, p, m=8):
    """returns None if some monomial has an exponent > 2 (then not a member: T_8 lives on reduced monomials)."""
    poly = product_coeffs(Z)
    vec = {}
    for e, c in poly.items():
        if max(e) > 2: return None
        vec[e] = (c * pow(multinom(m, e), -1, p)) % p
    return vec
