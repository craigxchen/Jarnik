# Exact inputs of Prop 6.2 recomputed from scratch: c spanning V(A,8), parity, Phi = sum c_k <k,u>^8,
# irreducibility over Q, sign witnesses (exact), and a rigorous zeta <= 7 bound.
import json, itertools, math, random
import flint
from rmodels import slab_box, moment_matrix
cuts = [frozenset(c) for c in json.load(open('r_sdc7_cuts.json'))]
n = 8; base = 7
A, others = slab_box(n, cuts, 1, 2, base=base)
rows, _ = moment_matrix(A, others, 8, 1)
M = flint.fmpz_mat(rows)
print("|A| =", len(A), " exact rank Poly_<8|_A =", M.rank(), " => dim V(A,8) =", len(A) - M.rank())
rows9, _ = moment_matrix(A, others, 9, 1)
print("exact rank Poly_<9|_A =", flint.fmpz_mat(rows9).rank(), "(= |A| means rho = 8)")
N, nul = M.nullspace()
c = [int(N[i, 0]) for i in range(len(A))]
g = 0
for x in c: g = math.gcd(g, x)
c = [x // g for x in c]
idx = {k: i for i, k in enumerate(A)}
even = all(c[idx[k]] == c[idx[tuple(-t for t in k)]] for k in A)
print("primitive c: sup-norm", max(map(abs, c)), " support", sum(1 for x in c if x), " even:", even)
# Phi(u) = sum_k c_k <k',u>^8  with u in coordinates u_x = delta_x - delta_base (x in others)
U = flint.fmpz_mpoly_ctx.get(tuple(f"u{i}" for i in range(7)), 'degrevlex')
us = U.gens()
Phi = U.from_dict({})
for ck, k in zip(c, A):
    if ck:
        lin = sum(k[x] * us[i] for i, x in enumerate(others))
        Phi += ck * lin ** 8
cont = Phi.content()
print("Phi != 0:", not Phi.is_zero(), " content(Phi) =", cont, " 8! =", math.factorial(8))
gpoly = Phi / cont if cont != 0 else Phi
print("g = Phi/content: terms", len(list(gpoly.terms())), " coefficient set", sorted(set(int(v) for v in gpoly.coeffs())))
const, facs = gpoly.factor()
print("factorisation of g over Q: const", const, " factors (degree, mult):", [(f.total_degree(), e) for f, e in facs])
# sign witnesses: search small integer points for positive and negative values
def ev(p, pt):
    return int(p(*[flint.fmpz(v) for v in pt]))
random.seed(1)
pos = neg = None
vals = []
for pt in itertools.product(range(-2, 3), repeat=7):
    if not any(pt): continue
    v = ev(gpoly, pt)
    if v > 0 and pos is None: pos = (pt, v)
    if v < 0 and neg is None: neg = (pt, v)
    vals.append(v)
print("on [-2,2]^7 \\ 0: #positive", sum(1 for v in vals if v > 0), " #negative", sum(1 for v in vals if v < 0), " #zero", sum(1 for v in vals if v == 0))
print("first negative witness", neg)
if pos is None:
    # random search on larger integer points
    for trial in range(200000):
        pt = tuple(random.randint(-30, 30) for _ in range(7))
        v = ev(gpoly, pt)
        if v > 0:
            pos = (pt, v); break
print("positive witness", pos)
# rigorous zeta bound: every zonotope generator v (in the box: |v_x|<=2 for x != base) has
# sum_p |v(S_p)| >= s_min; zonotope in A => sum_i L_i sum_p |v_i(S_p)| <= 14 * range(k(S_p)) = 14*4.
smin = None
for vp in itertools.product(range(-2, 3), repeat=7):
    if not any(vp): continue
    v = [0] * n
    for x, t in zip(others, vp): v[x] = t
    v[base] = -sum(vp)
    s = sum(abs(sum(v[x] for x in S)) for S in cuts)
    smin = s if smin is None else min(smin, s)
rng = [max(sum(k[x] for x in S) for k in A) - min(sum(k[x] for x in S) for k in A) for S in cuts]
print("min_v sum_p |v(S_p)| over generators =", smin, "; ranges of k(S_p) on A:", sorted(set(rng)),
      "; => zeta <=", sum(rng) // smin, "(floor of", sum(rng), "/", smin, ")")
json.dump({'c': c, 'A': A, 'others': others}, open('r_sdc7_c.json', 'w'))
