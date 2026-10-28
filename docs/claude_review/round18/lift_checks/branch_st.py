"""Check: no pencil member f_lam has branch pair {s, t'} (Res_lam(Delta_lam(s), Delta_lam(t')) != 0), and list the
lam with s or t' as a branch value (genus-5 members) at t0 and item 109."""
import sympy as sp
from fractions import Fraction as Fr
from fam import family, disc_T, check
L = sp.symbols('lam'); Tt = sp.symbols('T')
import json
TS = [r['t'] for r in json.load(open('batch_beta7.json')) if r['ok']]
for xs in [[Fr(v) for v in t] for t in TS]:
    F, D, ok = check([Fr(x) for x in xs], verbose=False)
    e1, e2, e3, e4 = [sp.Rational(str(v)) for v in F['e']]
    s = -sp.Rational(str(D['b'][1])) / sp.Rational(str(D['b'][2])); tp = -sp.Rational(str(D['c'][0])) / sp.Rational(str(D['c'][1]))
    A2 = (1 + L * e2) + L * Tt; A1 = -L * e3 - L * e1 * Tt; A0 = L * e4 - Tt
    Dl = sp.expand(A1 ** 2 - 4 * A2 * A0)
    ps, pt = sp.Poly(Dl.subs(Tt, s), L), sp.Poly(Dl.subs(Tt, tp), L)
    lb, lc = sp.Rational(str(F['lam']['b'])), sp.Rational(str(F['lam']['c']))
    print('t =', [str(x) for x in xs], ' Res(Delta(s),Delta(t\')) != 0:', sp.resultant(ps, pt) != 0,
          ' roots Delta(s)=0:', sp.roots(ps), ' (lam_b =', lb, ') roots Delta(t\')=0:', sp.roots(pt), ' (lam_c =', lc, ')')
