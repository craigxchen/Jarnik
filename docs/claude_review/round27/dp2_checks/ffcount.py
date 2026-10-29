"""Count points of S: w^2 = f(x,y,z) over F_q, q = p^n (p odd), with numpy log tables.
#S(F_q) = sum_{pts of P^2(F_q)} (1 + chi(f(pt))).  For a rational surface: #S(F_q) = q^2 + q*Tr(Frob^n|NS) + 1.
"""
import numpy as np, itertools, sys, time
import flint

def _exp_table(p, n, mod):
    q = p**n
    exp = np.zeros(q - 1, dtype=np.int64)
    cur = [1] + [0]*(n - 1)
    pw = p ** np.arange(n)
    for k in range(q - 1):
        exp[k] = int(np.dot(cur, pw))
        carry = cur[-1]
        cur = [0] + cur[:-1]
        for i in range(n):
            cur[i] = (cur[i] - carry*mod[i]) % p
    return exp

def find_primitive_poly(p, n):
    """monic poly of degree n over F_p whose root has multiplicative order p^n - 1 (hence irreducible, primitive)"""
    q = p**n
    for coeffs in itertools.product(range(p), repeat=n):
        if coeffs[0] == 0:
            continue
        mod = list(coeffs) + [1]
        exp = _exp_table(p, n, mod)
        if len(np.unique(exp)) == q - 1 and 0 not in exp:
            return mod, exp
    raise ValueError

class GF:
    def __init__(self, p, n):
        self.p, self.n, self.q = p, n, p**n
        q = self.q
        mod, exp = find_primitive_poly(p, n)
        pw = p ** np.arange(n)
        self.exp = exp
        log = np.full(q, -1, dtype=np.int64)
        log[exp] = np.arange(q - 1)
        assert (log[1:] >= 0).all()
        self.log = log
        # digit decomposition for addition
        digits = np.zeros((q, n), dtype=np.int64)
        for i in range(n):
            digits[:, i] = (np.arange(q) // p**i) % p
        self.digits = digits
        self.pw = pw

    def add(self, a, b):
        da, db = self.digits[a], self.digits[b]
        return ((da + db) % self.p) @ self.pw

    def mul(self, a, b):
        res = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        a = np.broadcast_to(a, res.shape); b = np.broadcast_to(b, res.shape)
        nz = (a != 0) & (b != 0)
        res[nz] = self.exp[(self.log[a[nz]] + self.log[b[nz]]) % (self.q - 1)]
        return res

    def scal(self, k):
        """embed integer k into F_q"""
        return int(k % self.p)

    def power(self, a, e):
        res = np.zeros_like(a)
        if e == 0:
            return np.ones_like(a)
        nz = a != 0
        res[nz] = self.exp[(self.log[a[nz]] * e) % (self.q - 1)]
        return res

    def chi(self, a):
        out = np.zeros_like(a)
        nz = a != 0
        out[nz] = np.where(self.log[a[nz]] % 2 == 0, 1, -1)
        return out

def eval_form(F, monos, X, Y, Z):
    """monos: list of ((i,j,k), coeff int)"""
    acc = np.zeros(X.shape, dtype=np.int64)
    for (i, j, k), co in monos:
        cs = F.scal(co)
        if cs == 0:
            continue
        term = F.mul(F.mul(F.power(X, i), F.power(Y, j)), F.power(Z, k))
        term = F.mul(term, np.full_like(term, cs))
        acc = F.add(acc, term)
    return acc

def count_S(monos, p, n, chunk=150):
    F = GF(p, n)
    q = F.q
    total = 0
    allx = np.arange(q, dtype=np.int64)
    # affine chart z=1: points (x,y,1)
    for x0 in range(0, q, chunk):
        xs = allx[x0:x0+chunk]
        X = np.repeat(xs, q); Y = np.tile(allx, len(xs)); Z = np.ones_like(X)
        val = eval_form(F, monos, X, Y, Z)
        total += len(X) + int(F.chi(val).sum())
    # line z=0: (x,1,0) and (1,0,0)
    X = allx.copy(); Y = np.ones_like(X); Z = np.zeros_like(X)
    val = eval_form(F, monos, X, Y, Z)
    total += len(X) + int(F.chi(val).sum())
    val = eval_form(F, monos, np.array([1]), np.array([0]), np.array([0]))
    total += 1 + int(F.chi(val).sum())
    return total

if __name__ == '__main__':
    import sympy as sp
    from geom import default_f
    x, y, z = sp.symbols('x y z')
    f = default_f()
    monos = [((i, j, k), int(co)) for (i, j, k), co in sp.Poly(f*1, x, y, z).terms()]
    p = int(sys.argv[1]); nmax = int(sys.argv[2])
    for n in range(1, nmax + 1):
        t0 = time.time()
        N = count_S(monos, p, n)
        q = p**n
        tr = (N - q*q - 1) / q
        print(n, q, N, 'trace', tr, '%.1fs' % (time.time() - t0), flush=True)
