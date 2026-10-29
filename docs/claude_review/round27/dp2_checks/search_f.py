"""Search a small-coefficient quartic f with f(0,0,1)=1 such that:
 smooth mod 3 (Groebner over GF(3): leading monomials contain pure powers of x,y,z) => smooth over Q, good red. at 3;
 F3 (tangent cone of D at P) irreducible over Q; Res(F3,F4') != 0 (D irreducible, P on no line);
 then report Frobenius-3 traces on NS for n=1..8 (computed separately)."""
import random, sys
import sympy as sp
from sympy import symbols, Poly, groebner, diff, factor_list, resultant, expand
from geom import analyse
x, y, z, t = symbols('x y z t')

monos = [(i, j, 4 - i - j) for i in range(5) for j in range(5 - i) if not (i == 0 and j == 0)]

def smooth_mod(f, p):
    G = groebner([diff(f, v) for v in (x, y, z)], x, y, z, order='grevlex', modulus=p)
    pure = set()
    for g in G.exprs:
        mon = Poly(g, x, y, z).monoms(order='grevlex')[0]
        nz = [i for i, e in enumerate(mon) if e > 0]
        if len(nz) == 1:
            pure.add(nz[0])
    return len(pure) == 3

def candidate(seed):
    rnd = random.Random(seed)
    f = z**4 + sum(rnd.choice([-2, -1, 0, 0, 1, 1, 2])*x**i*y**j*z**k for (i, j, k) in monos)
    return expand(f)

if __name__ == '__main__':
    found = 0
    for seed in range(int(sys.argv[1]), int(sys.argv[2])):
        f = candidate(seed)
        if not smooth_mod(f, 3):
            continue
        out, w, T2hat, F4 = analyse(f)
        fl = out['F3_factor'][1]
        if len(fl) != 1 or fl[0][1] != 1 or Poly(fl[0][0], x, y).total_degree() != 3:
            continue
        if out['Res_F3_F4p'] == 0 or 0 in out['inf_check']:
            continue
        print(seed, 'f =', f, '| discF3 =', out['disc_F3'], 'real roots', out['F3_real_roots'], flush=True)
        found += 1
        if found >= 6:
            break
