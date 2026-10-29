"""K7: exact verification of the order-7 lattice zonotopes found by k6 in syl8 and sdc7 (L = 2):
(i) all 2^7 points k0 + sum t_i v_i (t in {0,1}^7) lie in A; (ii) c = prod_i (T_{v_i} - 1) delta_{k0} lies in
V(A,7) exactly (all moments of order <= 6 vanish) and not in V(A,8); (iii) P_c(Y) = Y^{k0} prod (Y^{v_i} - 1);
(iv) height |c|_inf and support size.  Exact integer arithmetic."""
import sys, itertools
sys.path.insert(0, '.')
from klib import build_box, monos_lt, mval
exec(open('k6_zono.py').read().split("model = sys.argv[1]")[0])   # reuse profile builders only (own file)
cases = {
 'syl8': ((0, 0, -1, -1, 1, -1, 0), [(0,0,0,0,0,0,1), (0,0,0,0,0,1,-1), (0,0,0,0,0,1,0), (0,0,0,1,-1,0,0),
                                     (0,0,0,1,0,0,0), (0,0,1,0,-1,0,0), (0,0,1,0,0,0,0)]),
 'sdc7': ((0, -1, 0, 0, 0, 0, 0), [(0,0,0,0,0,0,1), (0,0,0,0,0,1,-1), (0,0,0,0,1,-1,0), (0,0,0,1,-1,0,0),
                                   (0,0,1,-1,0,0,0), (0,1,-1,0,0,0,0), (0,1,0,0,0,0,0)]),
}
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for name, (k0, gens) in cases.items():
    slabs = [(S, 2) for S in (syl8() if name == 'syl8' else sdc7())]
    A = set(build_box(8, slabs))
    pts = [tuple(k0[x] + sum(t[i] * gens[i][x] for i in range(7)) for x in range(7)) for t in itertools.product((0, 1), repeat=7)]
    allin = all(p in A for p in pts)
    c = {}
    for t in itertools.product((0, 1), repeat=7):
        p = tuple(k0[x] + sum(t[i] * gens[i][x] for i in range(7)) for x in range(7))
        c[p] = c.get(p, 0) + (-1) ** (7 - sum(t))
    c = {k: v for k, v in c.items() if v}
    mom = {}
    ok7 = all(sum(v * mval(k, a) for k, v in c.items()) == 0 for a in monos_lt(7, 7))
    ok8 = all(sum(v * mval(k, a) for k, v in c.items()) == 0 for a in monos_lt(7, 8))
    # rank of generators (are they independent?) and which are in L_pair
    pr("%s: k0=%s gens=%s" % (name, k0, gens))
    pr("   all 128 zonotope points in A: %s; supp c size %d; |c|_inf = %d; c in V(A,7): %s; c in V(A,8): %s" % (
        allin, len(c), max(abs(v) for v in c.values()), ok7, ok8))
    import sympy as sp
    Y = sp.symbols('y0:7')
    Pc = sum(v * sp.prod([Y[x] ** k[x] for x in range(7)]) for k, v in c.items())
    prod = sp.prod([Y[x] ** k0[x] for x in range(7)]) * sp.prod([sp.prod([Y[x] ** g[x] for x in range(7)]) - 1 for g in gens])
    pr("   P_c == Y^k0 prod (Y^v_i - 1): %s" % (sp.simplify(sp.expand(Pc - prod)) == 0))
open("k7_zono_verify_out.txt", "w").write("\n".join(out) + "\n")
