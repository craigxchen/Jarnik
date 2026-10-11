"""milp2 with explicit support indicators y_S (stronger LP relaxation):
   om_S >= y_S, om_S <= Ubar y_S, sum_S y_S >= M-1 (rank), covering rows and cut deactivation written
   with y_S, lift rows of cuts use Bm * y_S for patterns outside Good(C).  Same validity argument as
   milp2 (every constraint holds at every feasible point with Omega <= Ubar, y_S = [om_S > 0]).
usage: python3 milp3.py M rho0 rho1 Ubar [time_limit] [strict] [img_cap]
"""
import sys, os, json
import numpy as np
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
import milp2


class Outer3(milp2.Outer):
    def solve(self, time_limit):
        n = self.n
        cut_y = []
        ny = 0
        for k in self.cuts:
            cut_y.append(list(range(ny, ny + len(k))))
            ny += len(k)
        nv = 2 * n + ny
        Y = n                      # offset of y_S
        Z = 2 * n                  # offset of cut binaries
        cons = []
        A = np.array([r for r, _ in self.int_rows]); b = np.array([rhs for _, rhs in self.int_rows]) + self.eps
        A1 = np.zeros((A.shape[0], nv)); A1[:, :n] = A
        cons.append(LinearConstraint(A1, lb=b, ub=np.inf))
        Cv = np.array(self.cover_rows)
        A2 = np.zeros((Cv.shape[0], nv)); A2[:, Y:Y + n] = Cv
        cons.append(LinearConstraint(A2, lb=1.0, ub=np.inf))
        L = np.zeros((2 * n + 2, nv))
        for S in range(n):
            L[S, S] = 1; L[S, Y + S] = -1                     # om - y >= 0
            L[n + S, S] = -1; L[n + S, Y + S] = self.Ubar      # Ubar y - om >= 0
        L[2 * n, :n] = 1                                       # Omega <= Ubar
        L[2 * n + 1, Y:Y + n] = 1                              # sum y >= M-1
        lb = np.concatenate([np.zeros(2 * n), [-np.inf, self.M - 1]])
        ub = np.concatenate([np.full(2 * n, np.inf), [self.Ubar, np.inf]])
        cons.append(LinearConstraint(L, lb=lb, ub=ub))
        if self.cuts:
            nrows = len(self.cuts) + sum(len(self.lifts[c]) for k in self.cuts for c in k)
            Ak = lil_matrix((nrows, nv)); bk = np.zeros(nrows)
            i = 0
            for k, ys in zip(self.cuts, cut_y):
                gmask = np.ones(n, dtype=bool)
                for c in k:
                    gmask &= self.good[c]
                for j in np.nonzero(~gmask)[0]:
                    Ak[i, Y + j] = 1.0
                for y in ys:
                    Ak[i, Z + y] = 1.0
                bk[i] = 1.0
                i += 1
                for c, y in zip(k, ys):
                    gm = self.good[c]
                    for row, r, Bm in self.lifts[c]:
                        for j in np.nonzero(row)[0]:
                            if gm[j]:
                                Ak[i, j] = row[j]
                            else:
                                Ak[i, Y + j] = Bm
                        Ak[i, Z + y] = -Bm
                        bk[i] = r - Bm + self.eps
                        i += 1
            cons.append(LinearConstraint(csr_matrix(Ak), lb=bk, ub=np.inf))
        lo = np.zeros(nv); hi = np.ones(nv); integ = np.ones(nv)
        hi[:n] = self.Ubar; integ[:n] = 0
        cobj = np.zeros(nv); cobj[:n] = 1
        res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lo, hi),
                   options={'time_limit': time_limit, 'mip_rel_gap': 1e-9, 'disp': False})
        return res


milp2.Outer = Outer3

if __name__ == '__main__':
    M = int(sys.argv[1]); rho0 = Fraction(sys.argv[2]); rho1 = Fraction(sys.argv[3]); Ubar = float(sys.argv[4])
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 3600
    strict = len(sys.argv) > 6 and sys.argv[6] == 'strict'
    cap = int(sys.argv[7]) if len(sys.argv) > 7 else 200
    os.makedirs('results', exist_ok=True)
    tag = f'M{M}_r{sys.argv[2].replace("/", "_")}_{sys.argv[3].replace("/", "_")}{"_s" if strict else ""}'
    out = milp2.run(M, rho0, rho1, Ubar, tl, strict, img_cap=cap, log=f'results/c_{tag}.json')
    print('RESULT', json.dumps(out, default=str))
