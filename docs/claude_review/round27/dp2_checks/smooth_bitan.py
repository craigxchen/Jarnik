"""(i) smoothness of the quartic f over Q and mod small primes (Groebner: leading monomials contain pure powers);
(ii) no line through p=(0:0:1) is a bitangent (so P lies on none of the 56 lines);
(iii) the 12 lines through p tangent to B={f=0}: field of definition.
"""
import sympy as sp
from sympy import symbols, groebner, Poly, expand, resultant, gcd, factor_list, diff
from geom import default_f, parts
x, y, z, t, m, c = symbols('x y z t m c')
f = default_f()

def zero_dim_only_origin(polys, gens, modulus=None):
    G = groebner(polys, *gens, order='grevlex', modulus=modulus) if modulus else groebner(polys, *gens, order='grevlex')
    lms = [Poly(g, *gens).monoms(order='grevlex')[0] for g in G.exprs]
    pure = {}
    for mon in lms:
        nz = [i for i, e in enumerate(mon) if e > 0]
        if len(nz) == 1:
            pure[nz[0]] = min(pure.get(nz[0], 99), mon[nz[0]])
    return len(pure) == len(gens), pure

grads = [diff(f, v) for v in (x, y, z)]
for mod in [None, 3, 5, 7]:
    ok, pure = zero_dim_only_origin(grads, (x, y, z), mod)
    print('smooth' if ok else 'NOT certified smooth', 'over', 'Q' if mod is None else 'F_%d' % mod, pure)

# (ii) bitangents through p. Line l_t: (x,y) = m*(1,t) (plus vertical line x=0: (x,y)=m*(0,1)).
def binary_quartic_coeffs(a1, b1):
    g = expand(f.subs({x: m*a1, y: m*b1, z: c}))
    P = Poly(g, m, c)
    return [P.coeff_monomial(m**j * c**(4-j)) for j in range(5)]

def bitangent_conditions(fj):
    f0, f1, f2, f3, f4 = fj  # f0 = 1
    beta = f1/2
    gam = (f2 - f1**2/4)/2
    return expand(f3 - 2*beta*gam), expand(f4 - gam**2)

fj = binary_quartic_coeffs(1, t)
assert fj[0] == 1
c1, c2 = bitangent_conditions(fj)
R = resultant(Poly(c1, t), Poly(c2, t))
print('bitangent-through-p resultant (nonzero => no bitangent through p with finite slope):', R)
fj0 = binary_quartic_coeffs(0, 1)
c10, c20 = bitangent_conditions(fj0)
print('vertical line conditions (must not both vanish):', c10, c20)

# (iii) tangent lines from p: binary quartic g_t(m,c) has a repeated root <=> disc_{m:c}(g_t) = 0
g_t = sum(fj[j]*m**j for j in range(5))  # dehomogenised c=1
disc = sp.discriminant(Poly(g_t, m))
disc = sp.factor(disc)
print('discriminant in slope t (tangent lines from p):', disc)
fl = factor_list(sp.numer(sp.together(disc)))
print('factor degrees:', [(Poly(q, t).degree(), e) for q, e in fl[1]])
