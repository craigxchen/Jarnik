"""Exact first-order (dual number) Jacobian of t=(x1,x2,x3) [x4=1] -> contact point B_C(t) in M_{0,5},
C = {x_j, a}: B_C = (x1,x2,x3,b'_j,c'_j) where b'_j, c'_j are the other preimages of x_j^2.
Coordinates on M_{0,5}: images of b', c' under the Moebius map x1->inf, x2->0, x3->1."""
from fractions import Fraction as Fr

class Du:
    def __init__(self, v, g):
        self.v = Fr(v); self.g = [Fr(x) for x in g]
    @staticmethod
    def c(v, n):
        return Du(v, [0] * n)
    def _w(self, o):
        return o if isinstance(o, Du) else Du(o, [0] * len(self.g))
    def __add__(self, o):
        o = self._w(o); return Du(self.v + o.v, [a + b for a, b in zip(self.g, o.g)])
    __radd__ = __add__
    def __neg__(self):
        return Du(-self.v, [-a for a in self.g])
    def __sub__(self, o):
        return self + (-self._w(o))
    def __rsub__(self, o):
        return self._w(o) - self
    def __mul__(self, o):
        o = self._w(o); return Du(self.v * o.v, [self.v * b + o.v * a for a, b in zip(self.g, o.g)])
    __rmul__ = __mul__
    def __truediv__(self, o):
        o = self._w(o); return Du(self.v / o.v, [(a * o.v - self.v * b) / (o.v * o.v) for a, b in zip(self.g, o.g)])
    def __rtruediv__(self, o):
        return self._w(o) / self
    def __pow__(self, n):
        r = Du(1, [0] * len(self.g))
        for _ in range(n): r = r * self
        return r

def contact_coords(x, j):
    x1, x2, x3, x4 = x
    e1 = x1 + x2 + x3 + x4
    e2 = x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4
    e3 = x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4
    e4 = x1 * x2 * x3 * x4
    def PQ(lam):
        return [lam * e4, -lam * e3, 1 + lam * e2], [1, lam * e1, -lam]
    lb = 4 * e4 / (e3 * e3 - 4 * e2 * e4); lc = -4 / (e1 * e1)
    xj = x[j]; T = xj * xj
    def other(lam):
        P, Q = PQ(lam)
        A1 = P[1] - T * Q[1]; A2 = P[2] - T * Q[2]
        return -A1 / A2 - xj
    bp, cp = other(lb), other(lc)
    M = lambda z: ((z - x2) * (x3 - x1)) / ((z - x1) * (x3 - x2))
    return M(bp), M(cp)

def rank(rows):
    rows = [list(r) for r in rows]; rk = 0; ncol = len(rows[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

if __name__ == '__main__':
    for base in ([2, 3, -5], [Fr(1, 2), 7, -3], [4, -7, 9]):
        x = [Du(base[0], [1, 0, 0]), Du(base[1], [0, 1, 0]), Du(base[2], [0, 0, 1]), Du(1, [0, 0, 0])]
        for j in range(3):
            u, v = contact_coords(x, j)
            J = [u.g, v.g]
            print('base', [str(b) for b in base], 'contact C={x_%d,a}:' % (j + 1), 'B_C =', (str(u.v), str(v.v)), ' Jacobian rank =', rank(J))
