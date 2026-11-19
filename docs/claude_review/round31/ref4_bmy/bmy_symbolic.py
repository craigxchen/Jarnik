"""Independent symbolic re-derivation of the ledger behind (IV.2).

Everything is derived from intersection numbers, not from Astra's formulas.
"""
import sympy as sp

a, b0, b, k, m, S, I, s, delta, Ebar, L, e = sp.symbols('a b0 b k m S I s delta E L e')

# ---------- 1. initial surface F_e, sections C_i = C0 + b_i F --------------
# K = -2C0-(2+e)F ; C_i^2 = 2b_i - e ; K.C_i = -2 - C_i^2 ; K^2 = 8
# (K + a sum C)^2 = 8 + 2a sum K.C_i + a^2 (S + 2I)
KB2 = 8 + 2*a*(-2*m - S) + a**2*(S + 2*I)
print('initial (K+aC)^2 :', sp.expand(KB2))
# virtual orbifold Euler number  e(X) - sum w e(D) + sum w w' D.D'
eorb = 4 - a*2*m + a**2*I
slack0 = sp.expand(3*eorb - KB2)
print('initial 3e-(K+B)^2 :', slack0)

# ---------- 2. one blow-up, derived from classes ---------------------------
# center on k boundary sections (weight a) and one vertical curve of weight b0,
# new exceptional weight b.  Pullback: C' = C - E, V' = V - E, K' = K + E.
# (K'+B')^2 = (K+B)^2 + (1 + b - b0 - a k)^2 * E^2  with E^2=-1
dKB2 = -(1 + b - b0 - a*k)**2
# Delta e_orb (virtual): e(Y)+1, new curve -2b, crossings:
#   each section moves its incidence from V (a*b0) to E (a*b); V'.E = 1 -> b*b0;
#   pair intersections among the k sections drop by binom(k,2) -> -a^2 binom(k,2)
deorb = 1 - 2*b + b*b0 + a*k*(b - b0) - a**2*k*(k-1)/2
dslack = sp.expand(3*deorb - dKB2)
# Astra's F re-adds a^2 per consumed pair (she charges a^2 per final transverse pair)
F = sp.expand(dslack + a**2*k*(k-1)/2)
Fastra = sp.expand(3*(1-2*b+b*b0)+3*a*k*(b-b0)+(1+b-b0-a*k)**2-2*a**2*k*(k-1)/2)
print('F matches Astra:', sp.simplify(F - Fastra) == 0)
A = sp.Rational(2, 3)
print('F(k,0,0)   =', sp.factor(F.subs({a: A, b0: 0, b: 0})))
print('F(2,0,a)   =', F.subs({a: A, b0: 0, b: A, k: 2}))
print('F(2,a,a)   =', F.subs({a: A, b0: A, b: A, k: 2}))
# chain of length e: total of F, minus 4/9 per consumed unit
chain = sp.Rational(8, 9) + (e-1)*sp.Rational(-4, 9)
print('chain total F =', sp.expand(chain), ' ; minus (4/9)e =', sp.expand(chain - sp.Rational(4, 9)*e))
big = sp.expand(F.subs({a: A, b0: 0, b: 0}) - sp.Rational(4, 9)*k*(k-1)/2)
print('k>=3 center, F-(4/9)C(k,2) =', sp.factor(big), ' = -(4/9)*', sp.factor(-big*sp.Rational(9, 4)))

# ---------- 3. four generic half fibres ------------------------------------
# fibre f, weight 1/2: e_orb: -1/2*2 + m*(a/2) ; (K+B)^2 += 2*(1/2)*f.(K+B) = -2 + a m
d_fib = sp.expand(3*(-1 + a*m/2) - (-2 + a*m))
print('per half fibre :', d_fib)
base = sp.expand((slack0 - a**2*I + 4*d_fib).subs(a, A))
print('constant + S part after 4 fibres :', base)

# ---------- 4. phi superadditivity -----------------------------------------
u, v = sp.symbols('u v')
psi = lambda x: x*(x-1)/2 + 2*x - 9
print('psi(u+v)-psi(u)-psi(v) =', sp.expand(psi(u+v) - psi(u) - psi(v)))
phi = lambda x: max(0, x*(x-1)//2 + 2*x - 9)
bad = [(x, y) for x in range(1, 300) for y in range(1, 300) if phi(x+y)-phi(x)-phi(y) > 2*x*y]
print('phi violations (u,v<300):', bad[:5])

# ---------- 5. final substitution ------------------------------------------
# 0 <= 8/9 S + 4/9 (I - sum C(s_h,2)) + sum(4 - 8 s_h/9) + 8/9 sumE, I=(m-1)S/2,
# S <= sum s_h - mL/4 + m delta (coefficient of S positive)
coefS = sp.expand(sp.Rational(8, 9) + sp.Rational(4, 9)*(m-1)/2)
print('coefficient of S :', sp.factor(coefS))
perlayer = sp.expand(coefS*(s - m/4) - sp.Rational(4, 9)*s*(s-1)/2 + 4 - sp.Rational(8, 9)*s)
target = 4 - m/6 - sp.Rational(2, 9)*(s - m/2)**2
print('per-layer identity holds:', sp.expand(perlayer - target) == 0)
print('delta coefficient:', sp.factor(coefS*m))
# multiply by 9/2 and rearrange
print('constant per layer *9/2 :', sp.expand(sp.Rational(9, 2)*(4 - m/6)), '= -(3/4)(m-24)')
mm = sp.symbols('mm', positive=True)
r = 3*(mm-24)/(4*mm*(mm+3))
crit = sp.solve(sp.diff(r, mm), mm)
print('argmax of 3(m-24)/(4m(m+3)):', crit, [sp.nsimplify(r.subs(mm, c)) for c in crit],
      'r(48)=', r.subs(mm, 48), 'r(47)=', r.subs(mm, 47), 'r(49)=', r.subs(mm, 49))
