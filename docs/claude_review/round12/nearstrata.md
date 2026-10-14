# Round 12: the quantitative half of the T = 1 programme (near-strata Lojasiewicz bound and arithmetic separation)

## 0. Summary

**No growth improvement and no uniform bound is proved.** `M <= (2/3 + eps) log R/loglog R` stands. The task was the
quantitative half of the `T = 1` programme of R11 §8.2: (1) a Lojasiewicz-type lower bound for
`kappa(u) = sup_c |Phi_c(u)|/(|c|_1 d^m)` away from the strata, (2) the arithmetic separation lemma, (3) the assembly. The
answers, with their status:

1. **What kappa is (Prop 2.1, proved).** `d^m kappa(u)` is the best uniform approximation error of `<k', u>^m` on the point set
   `A` by polynomials of degree `< m` (exactly for `V(A,m)`; with `V_H^perp` in place of `Poly_{<m}|_A` for `V_H`).
2. **Effective Lojasiewicz inequalities in two models (Theorems 2.4-2.5, proved; exact check of the construction).** On the
   uncut ternary cube and on every one-cut box, `kappa_H(u) >= (m!/(2d)^m) (dist_inf(u, TC_H)/(n-1))^m`, with the tangent cone an
   explicit finite union of rational subspaces. Exponent `m`, constant `e^{-O(m ln n)}`. These boxes are zonotopal, so this
   mechanism (products of differences) never applies to an actual cluster in the regime (Remark 2.10).
3. **The right qualitative input is ZE-TC^braid, and degree bounds cannot quantify it (Props 2.6-2.7, proved).**
   `Lambda_A(g) := min{kappa_H(u) : |u|_inf = 1, Gamma(u) <= g}` (`Gamma` = normalised pair defect) is positive for all `g` iff the real
   tangent cone lies in the braid arrangement; this gives a collision Lojasiewicz inequality `kappa_H >= gamma G^theta`
   ineffectively. An explicit pair of degree-`m` integer forms with zero set in the braid arrangement has
   `Lambda <= e^{-(m-2) N_2 g/2 + O(m N_2 ln n)}`, so **no effective Lojasiewicz theorem depending only on dimension, degree,
   heights and the qualitative zero set (Solerno/Kollar type) can supply what is needed** (exponent `d/2` in `g`).
4. **Arithmetic separation (Lemmas 3.1-3.3, Theorem 3.4, Corollaries 3.5, 3.7; proved, any cluster, any norm).** The actual
   direction satisfies `Gamma(hat u) <= abar = W/(2(n-1)) + 2 ln C'` (a discriminant bound from the pair gaps and the Plotkin
   average); hence it is at distance `>= (1/2) e^{-N_2 abar/(2P(L))}` from every subspace with `P(L)` collision pairs (Markov's
   inequality, replacing the judge's heuristic). Every cheap span at `m' = ceil(4w)` has `P >= ~n^2/32`, so the separation is
   `n^{-O(1)}` for near-extremal clusters and `n^{-(8/c)(1+o(1))}` for any cluster with `n >= c log R/loglog R`; and it beats the
   smallest radius any near-strata statement can have near cheap spans for `n >= 5185` (`r = n`, scan to `2 * 10^4`).
5. **Assembly (Theorem 4.1, Corollary 4.2, Theorem 4.3; proved implications).** Step 0 + the discriminant bound + the single
   inequality **(MIN)** `Lambda_A(abar) > C d e^{-W/4}/(m'+1)` contradict [D] Thm 4.1 at `T = 1`. If (MIN) holds for all
   regime profiles on `n >= n_1` points, then `n = O(log W) = O(loglog R)` **on all norms** (no rank lemma, no nondegenerate
   subcluster). (MIN) follows from a collision Lojasiewicz inequality with `theta < d` and `ln(1/gamma) <= (3/8 - o(1)) n ln n`,
   or from R11's ZE-TC_eps in a corrected form (all rational strata, each with `>= N_2/Q_0` collision pairs).
6. **The minimal open statement is MIN\*(A)** (§4.3): `Lambda_A(g) > (d/(m'+1)) e^{-dg/2}` for `g >= 1.5 ln n - 7` (`C <= 1`), for the
   `T = 1` slab box at `m' = ceil(4w)`. It is minimal relative to the arithmetic used: every `u` with `Gamma(u) <= abar` is
   compatible with all averaged pair data.
7. **Necessary conditions and obstructions (Props 5.1-5.3, proved; thresholds by scan).** MIN\* forces every linear subspace `L`
   of the tangent cone with transversal order `o` to have `P(L) >= o n/2` collision pairs; a collision-free rational stratum
   kills it. At orders `m'' in (2jk, m_max]` on Sylvester profiles (`n = 2^k`, nonempty from `k = 16`) **(MIN) is false**: the
   `k`-column span is collision-free and contains a point with `Gamma <= 3`. At `m' = ceil(4w)` every cheap span satisfies the
   necessary condition whatever its true transversal order for `n >= 473,716` (`r = n`), and with the proved orders for
   `n >= 5,185`.

| task item | outcome | where |
| --- | --- | --- |
| (1) explicit `f(eta, n, m)` | proved in the uncut and one-cut models (`f = (m!/(2d)^m)(eta/(n-1))^m`); in general only ineffective, and provably not obtainable from degree/height bounds; the needed qualitative input is ZE-TC^braid plus the size bound MIN\* | Thms 2.4-2.5, Props 2.6-2.7 |
| (2) separation lemma | **proved** (any cluster, any norm, growth regime included), for every subspace with many collision pairs, in particular every cheap span | Lemmas 3.1-3.3, Thm 3.4, Cors 3.5, 3.7 |
| (3) assembly | **proved implication**: MIN for regime profiles `=>` `n = O(loglog R)` on all norms | Thm 4.1, Cor 4.2, Thm 4.3 |
| minimal open statement | MIN\*(A) at `m' = ceil(4w)`: ZE-TC^braid plus `ln(1/Lambda_A(g)) < dg/2 - ln(d/(m'+1))` | §4.3 |
| growth / uniform bound | **not proved** | — |

Fake test. The inputs used about positions are the pair gaps and their Plotkin average, which every fake satisfies; the
input that fakes violate is [D] Thm 4.1 (rigidity), entering through Step 0. MIN\* is a statement about the profile's slab
box alone; a gap-admissible fake whose `hat u` violates Step 0's conclusion is no obstruction, while one whose `hat u` satisfies
`kappa_H(hat u) <= eps_0` with `Gamma(hat u) <= abar` would refute MIN\* (R11 §8.3 item 4 in this language).

## 1. Setting, conventions, imported facts

**Conventions.** *Theorem / Proposition / Lemma*: proved here with all quantifiers, or quoted with an exact
source. *Exact*: integer or exact-rational computation. *Evidence*: floating point or a finite experiment that no
proof uses. *Hypothesis*: unproved. Scripts are in `nearstrata_checks/` (each `x.py` writes `x_out.txt`).
[D] = `round10/determinant.md`; R11 = `round11/README.md`; `theory.md`, `compute.md` are the round-11 write-ups
(refereed). Two companion round-12 write-ups in this folder, `t1rigour.md` and `zetc.md`, are **unrefereed**; I
cite them only for comparison, and every statement of theirs that I use is re-proved here.

**Clusters.** `N >= 2`, `W = log N`; `n` distinct `z_x in Z[i]` with `|z_x|^2 = N` on an arc of length `<= C sqrt R`,
`R = sqrt N`; `C' = max(C, 1)`. `delta_x` = angle of `z_x` minus the arc midpoint; `Theta <= C e^{-W/4}` is the angular
width. For a split `p = pi conj(pi)` with `p^{e_p} || N`: `a_{p,x} = v_pi(z_x)`, `beta_p = (2 a_p - e_p 1)/e_p in [-1,1]^n`;
`r` = number of split primes; `d_xy = sum_p log p |a_{p,x} - a_{p,y}|` (the *pair distance*). Nothing below assumes
`N` squarefree, normalised or odd, and nothing assumes near-extremality unless stated.

**The T = 1 slab box** (R11 §5.3, J3). `d = n - 1`, base point `n`, `k' = (k_1, ..., k_d)`, `k_n = -sum k'`. For `r >= 1`:

```text
w = sqrt(8 d ln 4r),   j = floor(w/2),   m' = ceil(4w),   N_2 := binom(n, 2),
A = { k in Z^n_0 : k' in {-1,0,1}^d,  |<k, beta_p>| <= w for every split p },     U = |A|,  H = U^2,
m_max = largest m with 8 binom(m + d - 1, d) <= 3^d.
```

`V(A, m')` = the `c in Q^A` annihilating every polynomial of degree `< m'` in `k'`; `V_H` = the `Q`-span of the integer
`c in V(A, m')` with `|c|_inf <= H`. The *regime conditions* are

```text
(R1)  m' <= m_max,            (R2)  W > W_1 := 24 d ln 3 / m' + 8 ln(e d C' / m').
```

(R1) holds exactly for `n >= 2900` (`r <= n`), `n >= 3623` (`r <= n ln n`), `n >= 5809` (`r <= n^2`) in the scanned range
`n <= 12000`, and for all `n` with `n - 1 >= 316 ln 4r` (J3; `t1rigour.md` Thm 2.1 gives the closed form; both
unrefereed; I use (R1) as a hypothesis). (R2) is `O(sqrt(n/ln r))` and is far below `W >= 3 n ln n - O(n)`.

**Positions and forms.** `u_x = delta_x - delta_n` (`x <= d`), and **always `u_n := 0`**, so that `<k, delta> = <k', u>` and
pair differences are `u_x - u_y = delta_x - delta_y` for all `x, y in [n]`. For `u in R^d`:

```text
Phi_c(u) = sum_{k in A} c_k <k', u>^{m'},     kappa_c(u) = |Phi_c(u)| / (|c|_1 d^{m'}),     kappa_H(u) = sup_{c in V_H \ 0} kappa_c(u),
TC_H = { u in R^d : Phi_c(u) = 0 for all c in V_H }.
```

(`d = D/2` at `T = 1`, since `|k'|_1 <= d` on `A`.) `|.|_inf` is the sup norm on `R^d`, `hat u = u/|u|_inf`, and
`dist(u, L) = min_{y in L} |u - y|_inf`. For a linear subspace `L` of `R^d` its **collision pairs** are the pairs
`{x, y}` of `[n]` with `y_x = y_y` for all `y in L` (with `y_n = 0`); `P(L)` is their number.

**The pair-defect functional.** For `u in R^d` (with `u_n = 0`) put

```text
Gamma(u) = (2/N_2) sum_{x<y in [n]} ln^+ (1/|u_x - u_y|)   in [0, +infinity],        G(u) = ( prod_{x<y} |u_x - u_y| )^{1/N_2}.
```

`Gamma(u) = +infinity` iff `u` lies in the **braid arrangement** `Br = {u : u_x = u_y for some x != y in [n]}`. Since
`ln^+ >= ln`, `G(u) >= e^{-Gamma(u)/2}` always.

**Imported, refereed.**
* (I1) [D] Thm 4.1, proof: for integer `c` supported in `A` with `F_c(z) != 0`, `Norm F_c(z) >= N^{2d} prod_p p^{-h_{supp c}(a_p)}`, and
  `<k, a_p> = (e_p/2) <k, beta_p>` on `Z^n_0`.
* (I2) `theory.md` Lemma 2.4 (Taylor remainder for `c in V(S,m)`), Prop 2.1 (tangent cone), Prop 2.2 (spans of cheap
  families lie in the base locus and in `TC`).
* (I3) `compute.md` Prop 2.1 (zonotope certificate): if `k_0 + sum_i t_i v_i in A` for all integers `0 <= t_i <= L_i`
  (any `v_i in Z^d`, dependence allowed) and `0 <= alpha_i <= L_i`, `sum alpha_i = m`, then the coefficient vector `c` of
  `w^{k_0} prod_i (w^{v_i} - 1)^{alpha_i}` is an integer vector in `V(A, m)` with `|c|_1 <= 2^m`.

**Lemma 1.1 (T = 1 rigidity and Step 0).** Let a cluster satisfy (R1)-(R2). Then `V_H != 0`, and for every `c in V_H`:
`sum_k c_k e^{i<k,delta>} = 0` and `kappa_c(hat u) <= eps_0 := C d e^{-W/4} / (m'+1)`. Also `e^{-W/2} <= |u|_inf <= C e^{-W/4}`.

*Proof.* (This is J3 and `t1rigour.md` Thms 3.1, 4.2; the argument is [D] Thms 2.1, 4.1 with `T = 1`.) *Size of `A`:* for
`k'` uniform on `{-1,0,1}^d`, `<k, beta_p> = sum_{x<n} k_x (beta_{p,x} - beta_{p,n})` is a sum of independent centred terms
with ranges `<= 4`, so Hoeffding gives `P(|<k,beta_p>| > w) <= 2 e^{-w^2/(8d)} = 1/(2r)`; a union bound gives
`U >= 3^d/2`. *`V_H != 0`:* by (R1), `E := binom(m'+d-1, d) <= 3^d/8 <= U/4`. The moment matrix (rows `k'^alpha`,
`|alpha| < m'`) has entries in `{-1,0,1}`, so the Siegel count of [D] Thm 2.1(c) with `B = 1` gives
`dim V_H >= U(1.25 ln U - 0.275)/(2 ln U + 1.1) >= U/2` once `ln U >= 3.3`. *Rigidity:* let `c in V_H` be an integer
generator with `F_c(z) := N^d sum_k c_k e^{i<k,delta>} != 0` (a Gaussian integer, since `|k|_1 <= 2d` on `A`). Taylor:
`|F_c(z)| <= N^d |c|_1 (d Theta)^{m'}/m'!` (moments of order `< m'` vanish and `|<k,delta>| <= d Theta`). (I1):
`h_{supp c}(a_p) <= e_p w`, so `Norm F_c(z) >= N^{2d} e^{-wW}`. Comparing, `(m'/2 - w) W <= 2 ln|c|_1 + 2m' ln(dC') - 2 ln m'!`;
with `m' >= 4w`, `|c|_1 <= U H <= 3^{3d}` and `ln m'! >= m' ln(m'/e)` this gives `W <= W_1`, contradicting (R2). So
`F_c(z) = 0` on generators, hence on `V_H`. *Step 0:* `theory.md` Lemma 2.4 with `t = |u|_inf`, `D_{hat u} <= d` gives
`kappa_c(hat u) <= d |u|_inf/(m'+1)`, and `|u|_inf <= Theta <= C e^{-W/4}`. The lower bound on `|u|_inf` is Lemma 3.1 below
with `d_{xn} <= W`. QED

The **Chebyshev ceiling** (J2): `kappa_c(u) <= 2^{1-m'}` for every `|u|_inf <= 1` and every `c in V(A, m')`. So Step 0 is
informative because `eps_0 = e^{-W/4 + O(log n)}` lies far below `2^{1-m'} = e^{-O(sqrt(n log n))}`.

## 2. The Lojasiewicz side: what kappa is, what can be proved, what must be assumed

Throughout §2, `A` is a finite subset of `{k in Z^n_0 : k' in {-1,0,1}^d}` (so `|<k',u>| <= d |u|_inf`), `m >= 1`, and `V_H`,
`Phi_c`, `kappa_c`, `kappa_H`, `TC_H` are defined from `(A, m)` as in §1 (with `m` in place of `m'`).

### 2.1 What kappa measures

**Proposition 2.1 (duality).** For every subspace `V` of `R^A` and every `g in R^A`,

```text
sup_{c in V \ 0} |<c, g>| / |c|_1  =  dist_inf(g, V^perp),        V^perp = { f in R^A : <c, f> = 0 for all c in V }.
```

Consequently, with `l_u(k) = <k', u>`,

```text
d^m kappa_H(u) = dist_inf( l_u^m |_A ,  V_H^perp ),        d^m kappa_V(u) = dist_inf( l_u^m |_A ,  Poly_{<m}|_A )   for V = V(A, m),
```

and `kappa_H <= kappa_V`, with equality when `V_H = V(A,m)`.

*Proof.* If `f in V^perp` then `|<c,g>| = |<c, g - f>| <= |c|_1 |g - f|_inf`; this gives `<=`. Conversely let
`r = dist_inf(g, V^perp) > 0` (if `r = 0` there is nothing to prove). The open convex set
`K = V^perp + {h : |h|_inf < r}` does not contain `g`, so some nonzero linear functional `c` has `<c, g> >= sup_K <c, .>`.
The supremum is finite only if `c` vanishes on the subspace `V^perp`, i.e. `c in (V^perp)^perp = V`, and then it equals
`r |c|_1`. So `<c, g> >= r |c|_1`. For the displayed consequences: `V(A,m)` is (over `R`) the annihilator of
`Poly_{<m}|_A`, so `V(A,m)^perp = Poly_{<m}|_A`; the supremum over the `Q`-span `V_H` equals the supremum over its real
span by continuity of `c -> |<c,g>|/|c|_1` on the real span minus `0`. QED

So `kappa` small means: **the `m`-th power of the linear form `<k', u>` nearly coincides on `A` with a polynomial of lower
degree** (in sup norm, relative to the trivial size `d^m |u|_inf^m`). The tangent cone `TC_V = {kappa_V = 0}` is the case of
exact coincidence (`compute.md` R2), and ZE-TC_eps is a statement about *near* coincidence.

**Lemma 2.2 (regularity).** `kappa_H` is the maximum of finitely many functions `|Phi_c|/d^m` (`c` ranging over the vertices
of the polytope `{c in V_H tensor R : |c|_1 <= 1}`). Hence it is continuous, semialgebraic, homogeneous of degree `m`, and
`TC_H` is a real algebraic cone.

*Proof.* `c -> |<c, l_u^m>|` is convex, so its maximum over a polytope is attained at a vertex; the polytope is the
intersection of a subspace with the cross-polytope and has finitely many vertices, independent of `u`. QED

### 2.2 Explicit certificates: products of differences

**Lemma 2.3 (zonotope lower bound).** Let `k_0 + sum_{i=1}^s t_i v_i in A` for all integers `0 <= t_i <= L_i`
(`v_i in Z^d \ 0`, dependence and repetition allowed), let `0 <= alpha_i <= L_i` with `sum_i alpha_i = m`, and suppose
`H >= 2^m`. Then for every `u in R^d`

```text
kappa_H(u)  >=  (m! / (2d)^m)  prod_i |<v_i, u>|^{alpha_i}.
```

*Proof.* By (I3) the coefficient vector `c` of `w^{k_0} prod_i (w^{v_i} - 1)^{alpha_i}` is an integer vector of `V(A,m)` with
`|c|_inf <= |c|_1 <= 2^m <= H`, so `c in V_H`. As formal power series in `z`, `sum_k c_k e^{<k',z>} = e^{<k_0,z>} prod_i (e^{<v_i,z>} - 1)^{alpha_i}`.
The right side starts in degree `m` with `prod_i <v_i,z>^{alpha_i}`; the degree-`m` part of the left side is `Phi_c(z)/m!`.
So `Phi_c(u) = m! prod_i <v_i,u>^{alpha_i}` and `kappa_c(u) >= m! prod |<v_i,u>|^{alpha_i} / (2^m d^m)`. QED

For the slab boxes of §1, `H = U^2 >= 2^{m'}` (as `U >= 3^d/2` and `m' <= m_max < d`), and the same holds in the models
below. Lemma 2.3 is the "explicit construction of `c` from the leading-form ideal" asked for in the task, in the only
form I can make explicit: products of first differences. Theorems 2.4-2.5 show that where such certificates suffice,
the Lojasiewicz inequality holds with **exponent `m` and constant `m!/(2d)^m`**, far better than needed (§4.3). Remark
2.10 shows why this can never be the mechanism for an actual cluster in the regime.

### 2.3 Two models with a complete effective Lojasiewicz inequality

**Theorem 2.4 (uncut ternary cube).** Let `A = B_1 = {k in Z^n_0 : k' in {-1,0,1}^d}`, `1 <= m <= 2d`, `H >= 2^m`, and
`Sigma_m = {u in R^d : |supp u| < m/2}` (a finite union of coordinate subspaces). Then `TC_H(B_1, m) = Sigma_m` and,
for every `u in R^d`,

```text
kappa_H(u)  >=  (m!/(2d)^m) dist_inf(u, Sigma_m)^m .
```

*Proof.* `Sigma_m` is contained in `TC_H`: for `|X| < m/2` the coordinate functionals `k -> k_x` (`x in X`) have width
`<= 2` on `B_1` and total width `2|X| < m`, so `R^X` is the span of a cheap family ((I2), Prop 2.2). For the inequality put
`rho = dist_inf(u, Sigma_m)` and assume `rho > 0`. Let `X = {x <= d : |u_x| >= rho}`. If `2|X| < m`, the vector `u 1_X` lies
in `Sigma_m` and `|u - u 1_X|_inf = max_{x notin X} |u_x| < rho`, a contradiction. So `|X| >= m/2`. Choose
`X' in X` with `|X'| = ceil(m/2)` and the zonotope `-sum_{x in X'} e_x + sum_{x in X'} [0,2] e_x` (it lies in `B_1`; `k_n`
is unconstrained at `T = 1`). Its lengths sum to `2|X'| >= m`; take `alpha_x = 2`, except `alpha_x = 1` for one `x` if `m` is
odd. Lemma 2.3 gives `kappa_H(u) >= (m!/(2d)^m) prod_{x in X'} |u_x|^{alpha_x} >= (m!/(2d)^m) rho^m`. In particular
`kappa_H > 0` off `Sigma_m`, so `TC_H = Sigma_m`. QED

(The equality `TC = Sigma_m` is R11 §8.1 / referee Prop R1 of round 10; the inequality is new.)

**One cut.** Fix `S in [d]` with `s = |S|`, an integer `1 <= L <= s - 1`, `R = [d] \ S`, and

```text
A_{S,L} = { k in Z^n_0 : k' in {-1,0,1}^d,  |k(S)| <= L },        k(S) = sum_{x in S} k_x .
```

This is the `T = 1` slab box of a profile with exactly one cutting column, a sign column whose side not containing the
base point is `S` (then `<k, beta> = 2k(S)` on `Z^n_0` and `L = j`). Define the rational subspaces

```text
(I)   Sigma^I_X   = { y : supp y in X },                                X in [d],  2|X| < m,
(II)  Sigma^II_X  = { y : supp(y - a 1_S) in X for some a in R },       X in [d],  2L + 2|X| < m,
(III) Sigma^III_Y = { y : supp(y|_R) in Y }   (y|_S arbitrary),         Y in R,    s + L + 2|Y| < m,
```

and let `Sigma(S,L,m)` be their union. For `u in R^d` (with `u_n = 0`) and `eta > 0`, the **`eta`-clusters** of `u` are the
classes of the equivalence relation on `[n]` generated by `|u_x - u_y| < eta`. Points in different clusters satisfy
`|u_x - u_y| >= eta`, and a cluster with `q` points has diameter `<= (q - 1) eta`, hence `< (n - 1) eta`. Let `C_0` be the cluster of `n` and

```text
sigma = |S \ C_0|,    t = s - max_C |S n C|  (max over all clusters, C_0 included),    a = |R \ C_0|,
nu_eta(u) = 2a + min(2 sigma, 2L + 2t, s + L).
```

**Lemma 2.5 (the strata lie in the tangent cone).** `Sigma(S, L, m)` is contained in `TC_H(A_{S,L}, m)` (and in the base locus).

*Proof.* (I) and (II) are spans of cheap families: `k -> k_x` has width `<= 2` and `k -> k(S)` has width `<= 2L` on `A_{S,L}`;
apply (I2) Prop 2.2. (III): let `Y' = S u Y` and `pi(k) = k_{Y'}`. Every function on `pi(A_{S,L}) = A2 x {-1,0,1}^Y`, with
`A2 = {k in {-1,0,1}^S : |k(S)| <= L}`, is a polynomial of degree `<= (s + L) + 2|Y|`: a function on a product is a sum of
products of functions of the factors; on `{-1,0,1}^Y` every function has degree `<= 2|Y|` (Lagrange in each coordinate);
on `A2` every function has degree `<= s + L` by Newton interpolation across the levels `Lambda_l = {k(S) = l}`
(`|l| <= L`). [Every function on `Lambda_l` has degree `<= s - |l|`: for `l >= 0` and `kappa in Lambda_l` with zero set `Z`
and minus set `M`, the polynomial `prod_{x in Z} (1 - k_x^2) prod_{x in M} (k_x^2 - k_x)/2` has degree `|Z| + 2|M| = s - l`,
equals `1` at `kappa`, and vanishes at the other points of `Lambda_l` (it forces `k = 0` on `Z`, `k = -1` on `M`, and then the
remaining `#{+1}` coordinates must sum to `l + |M|`, so all are `+1`); `l < 0` by symmetry. Order the levels
`0, 1, -1, ..., L, -L` and write `f = sum_i prod_{i' < i} (k(S) - l_{i'}) q_i`; the `i`-th term has degree
`<= i + s - ceil(i/2) <= s + L`.] For `u in R^{Y'}`, `l_u^m` is a function of `k_{Y'}`, hence agrees on `A_{S,L}` with a
polynomial of degree `<= s + L + 2|Y| < m`, which every `c in V(A,m)` annihilates; the same holds for `e^{i l_u}`. QED

**Theorem 2.5 (one cut: effective ZE-TC_eps).** Let `m >= 1`, `H >= 2^m`, `A = A_{S,L}`. For every `u in R^d` and `eta > 0`:
* (a) if `nu_eta(u) >= m`, then `kappa_H(u) >= (m!/(2d)^m) eta^m`;
* (b) if `nu_eta(u) < m`, then `dist_inf(u, Sigma(S,L,m)) < (n-1) eta`.

Consequently `TC_H(A_{S,L}, m) = Sigma(S, L, m)` and, for every `u in R^d`,

```text
kappa_H(u)  >=  (m!/(2d)^m) ( dist_inf(u, TC_H(A_{S,L}, m)) / (n-1) )^m .                                            (2.1)
```

*Proof.* (a) We build a lattice zonotope inside `A` from the following **gadgets**, placed on pairwise disjoint sets of
coordinates; each gadget lists generators, lengths, base point, its *degree* (sum of lengths) and its *cost* (the
largest `|k(S)|` over its points), and every generator `v` satisfies `|<v,u>| >= eta`:

| gadget | coordinates | generators (length) | base | degree | cost | factor `prod <v,u>` (all `>= eta` in modulus) |
| --- | --- | --- | --- | --- | --- | --- |
| R-single | `x in R \ C_0` | `e_x` (2) | `-e_x` | 2 | 0 | `u_x^2` |
| S-single | `x in S \ C_0` | `e_x` (2) | `-e_x` | 2 | 1 | `u_x^2` |
| pair | `x, y in S`, different clusters | `e_x - e_y` (2) | `e_y - e_x` | 2 | 0 | `(u_x - u_y)^2` |
| hexpair | `x, y in S \ C_0`, different clusters | `e_x, e_y, e_x - e_y` (1 each) | `-e_x` | 3 | 1 | `u_x u_y (u_x - u_y)` |

(The hexpair's eight points are `(-1,0), (0,0), (-1,1), (0,-1), (0,1), (1,-1), (0,0), (1,0)` in `(k_x, k_y)`; all have
`|k_x + k_y| <= 1`.) `x notin C_0` gives `|u_x| = |u_x - u_n| >= eta`. A sum of gadgets on disjoint coordinates has all
coordinates in `{-1,0,1}` and `|k(S)| <= ` total cost, so it lies in `A` when the total cost is `<= L`. By Lemma 2.3 it is
enough to find gadgets of total cost `<= L` and total degree `>= nu_eta(u)`. The `a` R-singles give `2a` at cost `0`. On `S`
the classes are the sets `S n C`; the *zero class* is `S n C_0`, of size `z = s - sigma`; `M = s - t` is the largest class
size. Let `tau = min(2 sigma, 2L + 2t, s + L)`.
* *`tau = 2 sigma`.* From `2 sigma <= s + L`, `sigma <= z + L`. Pair `min(sigma, z)` coordinates of `S \ C_0` with distinct
  coordinates of the zero class; use the remaining `sigma - min(sigma, z) <= L` coordinates of `S \ C_0` as S-singles.
  Degree `2 sigma`, cost `<= L`.
* *`tau = 2L + 2t < 2 sigma`.* A largest class `S n C*` is not the zero class (otherwise `t = sigma`). From
  `2L + 2t <= s + L`, `M >= t + L`. Pair the `t` coordinates outside `C*` with distinct coordinates of `S n C*`, and use `L`
  further coordinates of `S n C*` (not in `C_0`) as S-singles. Degree `2t + 2L`, cost `L`.
* *`tau = s + L < min(2 sigma, 2L + 2t)`.* Then `z < q/2` with `q = s - L >= 1`, and every class has size `< (s + L)/2`;
  `sigma > (s + L)/2 >= L`; and there are at least two nonzero classes (one nonzero class would give
  `sigma <= M < (s+L)/2`). We use the fact that a complete multipartite graph on `2h` vertices with all parts `<= h` has a
  perfect matching (list the vertices part by part and match the `i`-th with the `(i+h)`-th).
  * `q` even, `h = q/2`: the total excess `E = sum_{nonzero classes} (|class| - h)^+` is `<= L` (one exceeding class has
    excess `< (s+L)/2 - h = L`; two or more have total excess `<= s - 2h = L`). Use the excess coordinates and then
    `L - E` further coordinates of `S \ C_0` as S-singles; the remaining `q` coordinates have all classes `<= h` (the zero
    class has `z < h`), so they carry `h` pairs. Degree `2L + q = s + L`, cost `L`.
  * `q` odd, `h_0 = (q-1)/2`: now one exceeding class has size `<= (s+L-1)/2`, hence excess over `h_0` at most `L`, and two
    or more have total excess `<= s - 2h_0 = L + 1`; the zero class has `z <= h_0`. Remove `L + 1` coordinates of `S \ C_0`
    (possible: `sigma >= L + 1`): all excess coordinates, and in any case two coordinates `x, y` from different nonzero
    classes (from two exceeding classes if there are two; otherwise `x` from the exceeding class, if any, and `y` from
    another nonzero class). Make `{x, y}` a hexpair and the other `L - 1` removed coordinates S-singles; the remaining
    `q - 1 = 2h_0` coordinates have all classes `<= h_0` and carry `h_0` pairs. Degree `3 + 2(L-1) + (q-1) = s + L`, cost `L`.

(b) Suppose `nu_eta(u) < m`; every cluster has diameter `< (n-1) eta`. If `tau = 2 sigma`: `X = (S u R) \ C_0` has
`2|X| = 2a + 2 sigma < m`, and `y = u 1_X in Sigma^I_X` has `|u - y|_inf <= max_{x in C_0} |u_x - u_n| < (n-1) eta`. If
`tau = 2L + 2t`: with `C*` as above and `a* = u_{x*}` for some `x* in S n C*`, let `y = a*` on `S n C*`, `y = 0` on `R n C_0`,
`y = u` elsewhere; `y - a* 1_S` is supported on `X = (S \ C*) u (R \ C_0)`, `2L + 2|X| = 2L + 2t + 2a < m`, so
`y in Sigma^II_X`, and `|u - y|_inf < (n-1) eta`. If `tau = s + L`: `Y = R \ C_0`, `y = u` off `R n C_0` and `0` on it;
`s + L + 2|Y| < m`, `y in Sigma^III_Y`, `|u - y|_inf < (n-1) eta`.

The consequences: by (a)-(b), `kappa_H(u) > 0` whenever `dist(u, Sigma) > 0` (take `eta = dist/(n-1)`), so
`TC_H in Sigma`, and Lemma 2.5 gives equality; (2.1) is (a) with `eta = dist_inf(u, Sigma)/(n-1)`, since then (b) fails. QED

*Remarks.* (1) The equality `TC = Sigma(S,L,m)` is `zetc.md` Cor 2.2(1) (unrefereed there; proved independently here,
with the slightly different gadget set {single, pair, hexpair}). (2) The exponent `m` in (2.1) is the natural one: along a
line `u_0 + eta v` not contained in `TC_H`, some `Phi_c` restricts to a nonzero polynomial of degree `<= m` in `eta`, so the
vanishing order of `kappa_H` along any line is `<= m`; Theorems 2.4-2.5 show that in these models the order `m` is also
achieved *uniformly*, with an explicit constant. (3) The loss `(n-1)^{-m}` comes only from the
chaining in the `eta`-clusters and is irrelevant at the scales of §4 (`e^{m ln n}` against `e^{(3/4) n ln n}`).

*Exact check* (`one_cut_gadgets.py`): for 3000 random instances (`d <= 7`, all `1 <= L < s <= d`, rational `u` with planted
clusters, `eta = 1/10`; all four cases A, B, C-even, C-odd occur), the gadget zonotope of the proof lies in `A_{S,L}`, has
degree `nu_eta(u)`, all its generators have `|<v,u>| >= eta`, its product vector satisfies
`sum c_k <k',u>^m = m! prod <v,u>^{alpha}` with `|c|_1 <= 2^m` and annihilates the powers `<lambda,k'>^j`, `j < m`, for random
integer `lambda`; and when `nu_eta(u) < m` the vector `y` of part (b) lies in the claimed stratum within `(n-1) eta`.
(A check of the construction, not a substitute for the proof.)

### 2.4 The general case: what the qualitative hypothesis gives, and what it cannot give

For a general `T = 1` slab box no explicit certificate is available (Remark 2.10), and the lower bound must be assumed.
The right object is the following one-variable function. For `g in [0, +infinity)` put

```text
Lambda_A(g)  =  min { kappa_H(u) : u in R^d, |u|_inf = 1, Gamma(u) <= g } .                                        (2.2)
```

The minimum exists: `Gamma` is a finite sum of continuous maps into `[0, +infinity]`, so `{|u|_inf = 1, Gamma <= g}` is
compact, and `kappa_H` is continuous (Lemma 2.2). `Lambda_A` is non-increasing in `g`. The arithmetic of §3 says exactly that
the actual direction lies in `{Gamma <= abar}`; so `Lambda_A` plays the role of `f(eta, n, m)` in the task, with the
"distance from the strata" measured by the pair defect `Gamma`. Two hypotheses are natural:

> **ZE-TC^braid(A, m).** `TC_H(A, m)` is contained in the braid arrangement `Br = {u : u_x = u_y for some x != y in [n]}`
> (with `u_n = 0`).
>
> **CL(A, m; theta, gamma)** (collision Lojasiewicz inequality). `kappa_H(u) >= gamma G(u)^theta` for every `u` with
> `|u|_inf = 1`.

**Proposition 2.6 (qualitative versus quantitative).**
* (i) `Lambda_A(g) > 0` for every `g < infinity` **iff** ZE-TC^braid(A, m).
* (ii) ZE-TC^braid(A, m) implies CL(A, m; theta, gamma) for some `theta, gamma > 0` (ineffective). Conversely CL with any
  `gamma > 0` implies ZE-TC^braid.
* (iii) CL(A, m; theta, gamma) implies `Lambda_A(g) >= gamma e^{-theta g/2}` for all `g`.
* (iv) ZE-TC^braid is implied by ZE-TC^rat (R11 §8.2 corrected in `t1rigour.md` §6.5: every real point of `TC_H` lies in a
  rational subspace contained in `TC_H`) together with: every rational subspace contained in `TC_H` has a collision pair.

*Proof.* (i) If `u in TC_H \ Br`, `u != 0`, then `hat u` has `Gamma(hat u) < infinity` and `kappa_H(hat u) = 0`, so
`Lambda_A(Gamma(hat u)) = 0`. Conversely a minimiser `u` with `kappa_H(u) = 0` and `Gamma(u) <= g < infinity` lies in
`TC_H \ Br`. (ii) On the compact set `K = {|u|_inf = 1}`, `f = kappa_H` and `h = prod_{x<y} |u_x - u_y| = G^{N_2}` are
continuous semialgebraic, and ZE-TC^braid says `f^{-1}(0) in h^{-1}(0)`. The Lojasiewicz inequality for continuous
semialgebraic functions on a compact set (Bochnak-Coste-Roy, *Real Algebraic Geometry*, Cor. 2.6.7) gives `N', c > 0` with
`h^{N'} <= c f` on `K`; take `theta = N_2 N'`, `gamma = 1/c`. The converse is clear since `G = 0` exactly on `Br`.
(iii) `G(u) >= e^{-Gamma(u)/2}`. (iv) Every point of `TC_H` lies in a rational subspace inside `TC_H`, which has a collision
pair `{x,y}`, so the point has `u_x = u_y`. QED

So the **qualitative input that is needed is ZE-TC^braid**, not ZE-TC^rat: a rational subspace of the tangent cone without
collision pairs (and an actual cluster has no collisions) is as harmful as a non-rational component (Prop 5.1(i)
below makes this quantitative). The quantitative input is a lower bound for `Lambda_A` at one level `g = abar`; §4 shows
which bound suffices, and §5 which bounds are excluded.

**Proposition 2.7 (degree, height and the qualitative hypothesis do not suffice).** Let `n >= 4` and `m >= 3`. The two forms
`F_1 = (u_1 - 2u_2)^2 u_3^{m-2}` and `F_2 = (u_1 - u_2)^2 u_3^{m-2}` in `u_1, ..., u_d` have integer coefficients of absolute
value `<= 4`, degree `m`, and their common real zeros lie in `Br`. Yet, with `g_0(n) := 2 ln(2n)`, for every `g >= g_0(n) + 1`,

```text
min { max(|F_1(u)|, |F_2(u)|) : |u|_inf = 1, Gamma(u) <= g }  <=  exp( -(m-2) (N_2/2) (g - g_0(n)) ).
```

*Proof.* Common zeros: `u_3 = 0 = u_n`, or `u_1 = 2u_2` and `u_1 = u_2`, i.e. `u_1 = u_2 = 0 = u_n`; both lie in `Br`. Test
point: `u_1 = 1`, `u_2 = 1/2`, `u_3 = tau` with `0 < tau <= 1/8`, and `u_4, ..., u_d` in `[-1, -1/4]` with mutual gaps
`>= 1/(2n)` (they fit, as `(d-4)/(2n) < 3/4`). Then `|u|_inf = 1`, every pair other than `{3, n}` is at distance
`>= min(1/4, 1/(2n)) = 1/(2n)` and contributes `<= ln(2n)` to `(N_2/2) Gamma`, while `{3, n}` contributes `ln(1/tau)`. So
`Gamma(u) <= 2 ln(2n) + (2/N_2) ln(1/tau) = g` for `tau = exp(-(N_2/2)(g - g_0))`, and `tau <= 1/8` because
`(N_2/2)(g - g_0) >= N_2/2 >= 3 > ln 8`. At this point `F_1 = 0` and `|F_2| = tau^{m-2}/4`. QED

*Reading.* Any "effective Lojasiewicz inequality" whose constants depend only on the number of variables, the degree, the
heights of the coefficients and the qualitative position of the zero set (Solerno, *Effective Lojasiewicz inequalities in
semialgebraic geometry*, AAECC 1991; Kollar, *An effective Lojasiewicz inequality for real polynomials*, Period. Math. Hungar.
1999; quoted from memory) applies to `{F_1, F_2}` as well, so the lower bound it gives for `Lambda` can never exceed
`exp(-(m-2) N_2 (g - g_0)/2)`.
Section 4 needs `Lambda_A(abar) >= e^{-(d/2) abar + O(ln n)}`, i.e. an exponent `d/2`, against `(m-2) N_2/2` here: a factor
`~ m n/2` too weak. So **the compactness-plus-degree-bound route cannot reach the target, whatever the constants in those
theorems**; independently, their exponents are of size `m^{O(d)}`. The obstruction in Prop 2.7 is visible: the forms vanish
to order `m - 2` transversally to the hyperplane `{u_3 = 0}`, which carries only **one** collision pair. Proposition 5.1 shows
that this is the general mechanism: a lower bound of the needed strength forces, along every rational stratum `L` of the
tangent cone, `(transversal order) * N_2 / P(L) < d (1 + o(1))`. Any proof must therefore use the structure of the slab
system, through the collision counts of its strata.

**Remark 2.8 (what an effective proof would have to supply).** By Prop 2.6(iii), Theorem 4.1 below needs
`CL(A, m'; theta, gamma)` with

```text
theta < d     and     ln(1/gamma) + theta ln C'  <  (1 - theta/d) W/4 - ln(C d/(m'+1)) ,
```

e.g. `theta <= d/2` and `ln(1/gamma) <= (3/8 - o(1)) n ln n` at `W ~ 3 n ln n`. Against the Chebyshev ceiling `2^{1-m'}`
this allows a loss of `e^{(3/8) n ln n}`: the inequality is enormously weaker than the model bounds of Theorems 2.4-2.5
(which have `ln(1/gamma) = O(m ln n)`), but its exponent is constrained: `theta < d` is far below anything compactness gives.

**Remark 2.9 (exponent versus constant).** In Theorems 2.4-2.5 the exponent with respect to `dist_inf(u, TC_H)` is `m`. In
`G`-form, a stratum `L` with `P(L)` collision pairs contributes `G <~ eta^{P(L)/N_2}` at distance `eta`, so the `G`-exponent
needed near `L` is (transversal order) `* N_2/P(L)`. For the strata of Theorem 2.5, `P(L) >= (n - m)(n - m - 2)/4` (Corollary 3.6), so
the `G`-exponent is `<= 2m (1 + O(m/n))`, which is `< d` once `n > 2m (1 + o(1))`; so CL holds in the one-cut models with
`theta ~ 2m` and `ln(1/gamma) = O(m ln n)` (Corollary 3.6).

**Remark 2.10 (explicit product certificates never apply to an actual cluster).** Let a cluster satisfy (R1)-(R2) and let
`Z = k_0 + sum_i [0, L_i] v_i` be a lattice zonotope in its own slab box with `sum L_i >= m'`. Then some generator has
`<v_i, delta> = 0`. Indeed `h_Z(beta_p) = sum_i L_i |<v_i, beta_p>| <= h_A(beta_p) <= 2w`, and with
`X_v = sum_p log p |<v, a_p>|`, `sum_i L_i X_{v_i} = sum_p (e_p log p/2) h_Z(beta_p) <= w W <= m' W/4`, so some `X_{v_i} <= W/4`;
if `<v_i, delta> != 0`, the character inequality (`determinant.md` Prop 5.3(iii), constant as corrected by its referee)
gives `X_{v_i} >= W/2 - 2 ln(C |v_i|_1 / (2 sqrt 2))`, and `|v_i|_1 <= 4d` at `T = 1` (both `k_0` and `k_0 + L_i v_i` lie in the
box), so `W <= 8 ln(sqrt 2 d C')`, contradicting `W >= 3n ln n - (K_0 + 4 log^+ C) n` (`two_thirds` (T3)) once `n >= n(C)`.
(This is `compute.md` Prop 2.2 at `T = 1`.) Hence at the
actual `hat u` Lemma 2.3 always gives `0`, for nondegenerate `delta` there is no such `Z` at all, and the boxes of
Theorems 2.4-2.5 (which are zonotopal: `zeta = rho`) are never the slab boxes of actual clusters. They are models, in which
the quantitative statement can be checked; they are not evidence about regime boxes.

## 3. The arithmetic side: separation of the actual point (any cluster, any norm)

### 3.1 The discriminant bound

**Lemma 3.1 (pair gap).** For an actual cluster and `x != y` in `[n]`: `e^{-d_xy/2} <= |delta_x - delta_y| <= Theta <= C e^{-W/4}`.
Hence `eps_xy := d_xy - W/2 + 2 ln C' >= 0` and `|hat u_x - hat u_y| >= e^{-eps_xy/2}`.

*Proof.* `g = gcd(z_x, z_y)` has `v_pi(g) = min(a_{p,x}, a_{p,y})` and `v_{conj pi}(g) = e_p - max(a_{p,x}, a_{p,y})` at each split
`p`, and the full power of every inert prime and of `1 + i`; so `Norm g = N e^{-d_xy}`. Since `g | z_x - z_y != 0`,
`N e^{-d_xy} <= |z_x - z_y|^2 = 4N sin^2((delta_x - delta_y)/2) <= N (delta_x - delta_y)^2`. Both points lie on the arc, so
`|delta_x - delta_y| <= Theta`. Finally `|hat u_x - hat u_y| = |delta_x - delta_y|/|u|_inf >= e^{-d_xy/2}/(C e^{-W/4}) = (C'/C) e^{-eps_xy/2}`. QED

**Lemma 3.2 (discriminant of the angles).** For an actual cluster, `sum_{x<y} d_xy <= n^2 W/4`. Consequently

```text
sum_{x<y} eps_xy  <=  N_2 abar,      abar := W/(2(n-1)) + 2 ln C',      Gamma(hat u) <= abar,
prod_{x<y} |delta_x - delta_y|  >=  e^{-n^2 W/8}.
```

*Proof.* For each split `p`, `sum_{x<y} |a_{p,x} - a_{p,y}|` is convex in each `a_{p,x} in [0, e_p]`, hence maximal when every
`a_{p,x} in {0, e_p}`, where it is `e_p * #{0} * #{e_p} <= e_p n^2/4`; multiply by `log p` and sum, using
`sum_p e_p log p <= W`. Then `sum eps_xy <= n^2 W/4 - N_2 W/2 + 2 N_2 ln C' = n W/4 + 2 N_2 ln C' = N_2 abar`. By Lemma 3.1,
`ln^+(1/|hat u_x - hat u_y|) <= eps_xy/2`, so `Gamma(hat u) <= (1/N_2) sum eps_xy <= abar`. The product bound is Lemma 3.1
multiplied over pairs. QED

So the actual direction lies in the compact set `{|u|_inf = 1, Gamma(u) <= abar}` of §2.4. This is all the pair
inequalities say about positions, once averaged; `abar = W/(2(n-1)) + O(1)` is `Theta(log n)` for near-extremal clusters
(`round8/stability.md` (4.2)), and `<= (1/c + o(1)) ln n` for any cluster with `n >= c log R/loglog R` (Cor 3.5).

### 3.2 Separation from subspaces with many collision pairs

**Lemma 3.3 (separation; any `u` with bounded pair defect).** Let `u in R^d`, `|u|_inf = 1`, `Gamma(u) <= g`.
* (i) For every `lambda > 0`, at most `N_2 g/(2 lambda)` pairs `{x, y}` have `|u_x - u_y| < e^{-lambda}`.
* (ii) For every linear subspace `L` of `R^d` with `P(L) >= 1` collision pairs,
  `dist_inf(u, L) >= (1/2) exp(-N_2 g/(2 P(L)))`.

In particular, for an actual cluster, `dist_inf(hat u, L) >= (1/2) exp(-N_2 abar/(2 P(L)))` for every such `L`.

*Proof.* `sum_{x<y} ln^+(1/|u_x - u_y|) = N_2 Gamma(u)/2 <= N_2 g/2`. (i): each such pair contributes more than `lambda`.
(ii): let `y in L` with `|u - y|_inf = rho = dist(u, L)`. For a collision pair `{x, x'}` of `L`, `y_x = y_{x'}` (with
`y_n = u_n = 0`), so `|u_x - u_{x'}| <= 2 rho`. If `2 rho >= 1` the claim is trivial; otherwise each of the `P(L)` pairs
contributes `>= ln(1/(2 rho))`, so `P(L) ln(1/(2rho)) <= N_2 g/2`. The last sentence is Lemma 3.2. QED

This is the arithmetic separation lemma asked for in (2), in a form that needs no case analysis: R11 §8.2 asked for a
pair orthogonal to `L` (a collision pair of `L`) with small excess, argued by "Markov's inequality on the pair excesses"
as a heuristic. Lemma 3.2 makes the average excess a theorem, and Lemma 3.3(ii) is Markov's inequality applied to it.
(The same argument is `t1rigour.md` Lemma 6.6.) What remains is to count collision pairs.

**Theorem 3.4 (collision pairs of cheap spans at T = 1; any profile).** Assume `2j <= n`. For `v in R^n` let `v_(1) <= ... <= v_(n)` be
its sorted entries, `tau_j(v) = sum_{i=1}^j (v_(n+1-i) - v_(i))`, and for integer `v`: window `Win(v) = [v_(j), v_(n+1-j)]`,
window length `ell(v) = v_(n+1-j) - v_(j)`, and `e(v) = #{x : v_x notin Win(v)}`.
* (a) `S_w := {k in {-1,0,1}^n : sum k = 0, |k|_1 <= 2j}` is contained in `A`, for every profile. Hence `h_A(v) >= 2 tau_j(v)`.
* (b) For integer `v`, `tau_j(v) >= j ell(v) + e(v)`.
* (c) Let `F = (v_1, ..., v_s)` be integer directions with `d_F = sum_i h_A(v_i) <= m' - 1` (a cheap family), `E` the set of
  points `x` with `v_{i,x} notin Win(v_i)` for some `i`, `n' = n - |E|`, and `q_F = prod_i (ell(v_i) + 1)`. Then
  `sum_i (j ell(v_i) + e(v_i)) <= (m'-1)/2`, `|E| <= (m'-1)/2 - j sum_i ell(v_i)`, `q_F <= 2^{sum_i ell(v_i)}`, and if
  `n' >= q_F`,

  ```text
  P(L_F)  >=  n' (n' - q_F) / (2 q_F),          L_F = span_R{v_i} (read in u-coordinates).
  ```
* (d) If `w >= 10` and `m' = ceil(4w)`, then `sum_i ell(v_i) <= 4`, so `q_F <= 16` and `|E| <= (m'-1)/2`.

*Proof.* (a) `k in S_w` has `|k_x| <= 1` and `|<k, beta_p>| <= |k|_1 max_x |beta_{p,x}| <= 2j <= w`. If `P` and `Q` are the positions
of the `j` largest and the `j` smallest entries of `v` (disjoint since `2j <= n`), then `+-(1_P - 1_Q) in S_w` and
`<1_P - 1_Q, v> = tau_j(v)`. (b) For `i <= j`, `v_(n+1-i) - v_(i) = ell + (v_(n+1-i) - v_(n+1-j)) + (v_(j) - v_(i))` with both
brackets `>= 0`. An entry above the window occupies a sorted position `n + 1 - i` with `i <= j - 1` and has
`v_(n+1-i) - v_(n+1-j) >= 1` (integers); symmetrically below. Summing over `i` gives `j ell + e`. (c) The first two claims
are (a)-(b) summed, with `|E| <= sum_i e(v_i)`. A point `x notin E` has `v_{i,x}` in an integer interval of `ell(v_i) + 1`
values, so the *patterns* `(v_{i,x})_i` of points outside `E` take at most `q_F` values, and `ell + 1 <= 2^ell`. Every element
of `L_F` is `sum_i t_i v_i` modulo constants, so two points with the same pattern are a collision pair of `L_F`. With
class sizes `n_1, ..., n_q` (`q <= q_F`, `sum n_c = n'`), convexity gives `sum_c binom(n_c, 2) >= q binom(n'/q, 2) >= n'(n'-q_F)/(2q_F)`.
(d) `j >= w/2 - 1`, so `10 j >= 5w - 10 >= 4w > m' - 1` for `w >= 10`; hence `sum ell_i <= (m'-1)/(2j) < 5`. QED

(Theorem 3.4 is `t1rigour.md` Lemmas 6.1-6.2 and Theorems 6.3, 6.7, re-proved; I use it only through (c)-(d).)

**Corollary 3.5 (separation from all cheap spans; radii in the growth regime).** Let an actual cluster satisfy (R1)-(R2),
`w >= 10`, and `n > (m'-1)/2 + 16`. For every cheap family `F` at order `m' = ceil(4w)`,

```text
dist_inf(hat u, L_F)  >=  sigma_F := (1/2) exp(-Q_F abar/2),     Q_F = N_2/P_min(F) <= 16 n(n-1) / ( (n - (m'-1)/2)(n - (m'-1)/2 - 16) ),
```

so `Q_F <= 16 (1 + O(w/n))`. In particular:
* *near-extremal clusters* (`round8/stability.md` `H(epsilon)`): `abar <= (1/epsilon + o(1)) ln n`, so
  `sigma_F >= n^{-(8/epsilon)(1 + o(1))}`;
* *any cluster in the growth regime* `n >= c log R/loglog R` (fixed `c in (0, 2/3]`, `W = 2 log R`): `abar <= (1/c + o(1)) ln n`
  (as `W/(2(n-1)) <= (1 + o(1)) loglog R/c` and `loglog R = (1 + o(1)) ln n`), so `sigma_F >= n^{-(8/c)(1+o(1))}`.

No Plotkin-tightness, smoothness, prime count or squarefree hypothesis is used: only distinctness, the arc, and
`sum_p e_p log p <= W`.

*Proof.* Theorem 3.4(c)-(d) with `n' >= n - (m'-1)/2 > 16 >= q_F`, then Lemma 3.3(ii) with `g = abar`. For near-extremal clusters,
`H(epsilon)` gives `W <= (2/epsilon) n loglog R` and `loglog R = (1 + o(1)) log n` (`stability.md` (1.2) and the remark after it).
QED

**Corollary 3.6 (the one-cut models satisfy the collision Lojasiewicz inequality).** In the setting of Theorem 2.5, with
`n >= m + 3`, every subspace in `Sigma(S, L, m)` has at least `(n - m)(n - m - 2)/4` collision pairs. Hence, for all `u` with
`|u|_inf = 1`,

```text
kappa_H(u)  >=  (m!/(2d)^m) (2(n-1))^{-m} exp(-m Q_1 Gamma(u)/2),        Q_1 = 2n(n-1)/((n-m)(n-m-2)),
kappa_H(u)  >=  (m!/(2d)^m) (n-1)^{-m} (G(u)/2)^{m Q_1},
```

i.e. CL(A_{S,L}, m; theta, gamma) with `theta = m Q_1 = 2m (1 + O(m/n))` and `ln(1/gamma) = O(m ln n)`.

*Proof.* Generic points: of `Sigma^I_X`, the values on `([d] \ X) u {n}` all vanish, and `n - |X| > n - m/2`; of `Sigma^II_X`,
the values on `S \ X` are equal and those on `(R \ X) u {n}` vanish, two groups of total size `n - |X| > n - m/2`, giving
`>= 2 binom((n - m/2)/2, 2)` pairs; of `Sigma^III_Y`, the values on `(R \ Y) u {n}` vanish, and `s + |Y| < m`, so this group
has more than `n - m` points. Each count is `>= (n-m)(n-m-2)/4`. Then Lemma 3.3(ii) bounds `dist(u, Sigma)` from below, and
(2.1) gives the first display. For the second, if `dist(u, L) = rho` then the `P = P(L)` collision pairs have
`|u_x - u_y| <= 2 rho` and the others `<= 2`, so `G^{N_2} <= (2 rho)^P 2^{N_2 - P}`, i.e. `(G/2)^{N_2} <= rho^P`; as `G/2 <= 1`
and `N_2/P <= Q_1`, `rho >= (G/2)^{Q_1}`; insert in (2.1). QED

**Corollary 3.7 (the separation beats the smallest admissible near-strata radius).** By Lemma 5.0 below, at distance `eta`
from a cheap span `L_F` (transversal order `o_F >= m' - d_F`), `kappa_H <= 4^{m'} eta^{o_F}`; so any near-strata statement of
the form "`kappa_H > eps_0` at distance `>= rho_F` from `L_F`" must have `rho_F >= rho_F^min := (eps_0 4^{-m'})^{1/(m' - d_F)}` (the
R11 radius `(eps 2^m)^{1/(m - d_F)}` with the rigorous constant). Take `C <= 1`, `m' = ceil(4w)` and `W >= W_min(n) := 3n ln n - 14n`
(`two_thirds` (T3)). Then `sigma_F > rho_F^min` for **every** cheap family `F` whenever the linear-in-`W` inequality

```text
(m' - 2js - 2e)( Q(s,e) abar/2 + ln 2 )  <  W/4 + m' ln 4 - ln(d/(m'+1)),     Q(s,e) = 2^{s+1} N_2/((n-e)(n-e-2^s)),
```

holds for all `0 <= s <= b`, `0 <= e <= floor((m'-1)/2) - js`, at `W = W_min` and in its `W`-coefficient. *Scan*
(`sep_vs_radius.py`, step 1): this holds for every `n in [5185, 20000]` (`r = n`), `[6425, 20000]` (`r = n ln n`),
`[10280, 20000]` (`r = n^2`), and fails just below; the `W`-coefficient condition (condition B of `compat.py`) holds on the
coarse scan up to `4 * 10^6`. *Proof of the reduction to the scan:* Corollary 3.5 with `Q_F <= Q(s,e)`, and `d_F >= 2(js + e)`
(Theorem 3.4(c)); the inequality is `ln(1/sigma_F) < ln(1/rho_F^min)` with the worst admissible `(s, e)`. QED

This answers task item (2) in the form R11 §8.2 posed it: *every* cheap span has collision pairs of small excess
(Theorem 3.4 + Lemma 3.3), and the actual point is farther from `L_F` than the smallest radius any ZE-TC_eps near `L_F` can
have, from `n ~ 5 * 10^3` on (`r = n`). It does **not** make R11's ZE-TC_eps sufficient: that also needs the non-cheap
rational strata (§4.2, §5.3).

## 4. Assembly: the conditional theorem and the minimal open statement

Lemma 1.1 holds verbatim at every order `m'' in [ceil(4w), m_max]` (the Siegel count needs only `m'' <= m_max`, and
`W_1` is decreasing in the order, so (R2) at `ceil(4w)` implies it at `m''`). Write `Lambda_{A,m''}` for (2.2) at order
`m''`, and `eps_0(m'') = C d e^{-W/4}/(m''+1)`.

### 4.1 The theorem

**Theorem 4.1 (conditional contradiction; one cluster, any norm).** Let a cluster satisfy (R1)-(R2) with `n >= 5`, let
`m'' in [ceil(4w), m_max]`, let `A` be its own `T = 1` slab box, and `abar = W/(2(n-1)) + 2 ln C'`. If

```text
(MIN)        Lambda_{A,m''}( abar )  >  C d e^{-W/4} / (m'' + 1),
```

then the cluster does not exist.

*Proof.* `|hat u|_inf = 1` and `Gamma(hat u) <= abar` (Lemma 3.2), so `kappa_H(hat u) >= Lambda_{A,m''}(abar) > eps_0(m'')`
by (MIN); but `kappa_H(hat u) <= eps_0(m'')` (Lemma 1.1, Step 0). QED

**Corollary 4.2 (growth, conditional; all norms).** Fix `C > 0`. Suppose there is `n_1` such that every profile on
`n >= n_1` points (any `r`, any columns `beta_p in [-1,1]^n`) with `ceil(4w) <= m_max` satisfies (MIN) at `m'' = ceil(4w)`
for every `W > W_1`. Then every cluster with arc constant `C` and `W > W_1` has

```text
n  <  max( n_1, n_R(r) ),        n_R(r) := min{ n_0 : ceil(4 w(n, r)) <= m_max(n) for all n >= n_0 },
```

and `r <= W/ln 5`. With `n_R(r) = O(log r)` (J3; `t1rigour.md` Thm 2.1(e) gives `n - 1 >= 316 ln 4r` as a sufficient
condition, unrefereed), this is `n = O(log W) = O(loglog R)`, on **all** norms. (`W > W_1 = O(sqrt(n/log r))` holds for
`n >= n_2(C)` by `two_thirds` (T3): `3n log n <= W + (K_0 + 4 log^+ C) n`.)

*Proof.* If `n >= max(n_1, n_R(r))` then (R1) holds and Theorem 4.1 applies to the cluster's own profile. QED

Compared with [D] Thm 5.7 (`n = O(loglog R)` on squarefree norms under ZE), Corollary 4.2 needs neither the rank lemma
([D] Lemma 3.1) nor a nondegenerate subcluster ([D] Lemma 3.2): distinctness of the points enters through the pair
gaps, and the collision structure of the strata is built into (MIN). The price is that (MIN) is quantitative.

### 4.2 Which hypotheses imply (MIN)

**Theorem 4.3 (sufficient conditions).** In the setting of Theorem 4.1, each of the following implies (MIN).
* (a) **Collision Lojasiewicz.** CL(A, m''; theta, gamma) with `theta < d` and
  `ln(1/gamma) + theta ln C' + ln(C d/(m''+1))  <  (1 - theta/d) W/4`.
* (b) **Distance form with collision counts** (R11's ZE-TC_eps, made precise). There are a family `𝓛` of linear subspaces
  of `R^d` and `Q_0 >= 1` with `P(L) >= N_2/Q_0` for every `L in 𝓛`, such that `kappa_H(u) > eps_0(m'')` for every
  `|u|_inf = 1` with `dist_inf(u, L) >= sigma := (1/2) e^{-Q_0 abar/2}` for all `L in 𝓛`.
* (c) **Lojasiewicz form with exponent.** With `𝓛, Q_0` as in (b): `kappa_H(u) >= gamma' min_{L in 𝓛} dist_inf(u, L)^alpha`
  for all `|u|_inf = 1`, where `alpha Q_0 < d` and
  `ln(1/gamma') + alpha (ln 2 + Q_0 ln C') + ln(Cd/(m''+1))  <  (1 - alpha Q_0/d) W/4`.

*Proof.* (a) Prop 2.6(iii): `Lambda(abar) >= gamma e^{-theta abar/2} = gamma e^{-theta W/(4d) - theta ln C'}`, and the
inequality says this exceeds `eps_0`. (b) If `Gamma(u) <= abar` then `dist(u, L) >= sigma` for every `L in 𝓛` (Lemma 3.3(ii)),
so `kappa_H(u) > eps_0`. (c) implies (b): at distance `>= sigma` from all `L in 𝓛`, `kappa_H >= gamma' sigma^alpha`
`= gamma' 2^{-alpha} e^{-alpha Q_0 W/(4d) - alpha Q_0 ln C'}`, which exceeds `eps_0` by the inequality. QED

*Reading.* In (c), with `𝓛` = the cheap spans and `Q_0 = Q_F <= 16(1 + O(w/n))` (Corollary 3.5), the local exponent
allowed near `L_F` is `alpha < d/16`. R11's guess `rho_F ~ (eps 2^m)^{1/(m - d_F)}` corresponds to `alpha = m' - d_F` near `L_F`;
then (c) needs `(m' - d_F) Q_F < d`, which §5.2 shows is the binding scale condition (`n > (8 + o(1)) w`). Note that (b)-(c)
need **every** rational stratum of `TC_H` in `𝓛` with `P(L) >= N_2/Q_0`: cheap spans are not enough
(`t1rigour.md` Thm 6.11: projection-cheap subspaces `R^{S'}` lie in `TC_H` and are not cheap spans; they have
`>= binom(n - 1 - |S'|, 2)` collision pairs, so they are admissible members of `𝓛`). Whether every rational subspace of
`TC_H(A, ceil(4w))` has `>= c n^2` collision pairs is open (`t1rigour.md` (Q1)); (a) and (MIN) itself do not need it as
a separate statement.

### 4.3 The minimal open statement

For `C <= 1`, `two_thirds` (T3) gives `W >= 3 n ln n - K_0 n` with `K_0 < 13.91`, so the level `g = abar = W/(2d)` of
Theorem 4.1 satisfies `g >= 1.5 ln n - 7` and `C d e^{-W/4}/(m'+1) = (d/(m'+1)) e^{-d g/2}`. Hence:

> **Hypothesis MIN\*(A)** (for one profile's `T = 1` slab box at `m' = ceil(4w)`). For every `g >= 1.5 ln n - 7`,
>
> ```text
> Lambda_A(g)  =  min{ kappa_H(u) : |u|_inf = 1, Gamma(u) <= g }   >   (d/(m'+1)) e^{-d g/2} .
> ```

MIN\* for all profiles with `ceil(4w) <= m_max` on `n >= n_1` points implies, by Corollary 4.2, `n = O(loglog R)` for
clusters with `C <= 1` on all norms (for `C > 1` the right side becomes `(C d/(m'+1)) e^{-(d/2)(g - 2 ln C)}` and the
range follows from (T3) with `K_0 + 4 log^+ C` in place of `K_0`). For the growth question alone (`M = o(log R/loglog R)`) it suffices to have MIN\* at
levels `g <= (1/c + o(1)) ln n` for each fixed `c > 0` and `n >= n(c)` (Corollary 3.5), since only clusters with
`n >= c log R/loglog R` must be excluded.

**Why this is the minimal statement relative to the arithmetic used.** For any `u*` with `|u*|_inf = 1` and
`Gamma(u*) <= abar`, the numbers `eps*_xy = 2 ln^+(1/|u*_x - u*_y|)` are `>= 0`, satisfy `|u*_x - u*_y| >= e^{-eps*_xy/2}` and
`sum eps*_xy <= N_2 abar`. So `u*` is compatible with every consequence of Lemmas 3.1-3.2 that the argument uses (the
pair gaps and their Plotkin average). A proof that uses, about the position of `hat u`, only Step 0 and these pair data
must therefore exclude every such `u*`, which is (MIN). (One can add the individual bounds `eps_xy <= W/2 + 2 ln C'`,
i.e. restrict to `min_{x<y} |u_x - u_y| >= e^{-W/4}/C'`; this changes nothing below. Arguments that use more arithmetic,
such as character inequalities for relations `v` that are not pair differences, or the algebraicity of `e^{i delta}`
(ZE-alg), could need less than (MIN).) The near-extremal structure of `round8/stability.md` (smoothness, `Theta(n)`
primes, Plotkin tightness) gives no further information about positions beyond `abar = Theta(log n)`; its margin
statements (Prop 4.1(d)-(e) there) are the same average.

**Qualitative core.** `Lambda_A` is non-increasing, so MIN\* gives `Lambda_A(g) > 0` for every `g`, which by Prop 2.6(i) is
ZE-TC^braid(A, ceil(4w)). Conversely ZE-TC^braid gives `Lambda_A(g) > 0` but no size (Prop 2.6(ii), Prop 2.7). The **open content** of route 5 at
`T = 1` is therefore: ZE-TC^braid for regime boxes at `m' = ceil(4w)`, **plus** the size bound of MIN\*, whose shape is fixed
by Remark 2.8 and §5: exponent (in `G`) below `d`, constant `>= e^{-(3/8 - o(1)) n ln n}`.

## 5. Necessary conditions and obstructions (where (MIN) is false or cannot come from general principles)

### 5.1 Collision counts of strata are forced

**Lemma 5.0 (transversal vanishing; projection cost).** Let `v_1, ..., v_s in Z^d`, `pi(k) = (<k', v_i>)_i`, and let `rho` be such
that every function on `pi(A)` is a polynomial of degree `<= rho` (e.g. `rho = sum_i h_A(v_i)`). Put `L = span_R{v_i}` and
`o = m - rho`. Then for every `c in V(A, m)`, `y in L`, `v in R^d` and real `eta`,

```text
Phi_c(y + eta v) = sum_{i >= o} binom(m, i) eta^i sum_k c_k <k',y>^{m-i} <k',v>^i,
kappa_c(hat{y + eta v}) <= 4^m |eta|^o      if |y|_inf = 1, |v|_inf <= 1, |eta| <= 1/2 .
```

*Proof.* For `y = sum t_i v_i`, `<k', y> = <t, pi(k)>`, so `<k',y>^{m-i}` agrees on `A` with a polynomial of degree
`<= min(m - i, rho)` in `k'`; times `<k', v>^i` it has degree `< m` when `i < m - rho`, and `c` annihilates it. Each remaining
coefficient is at most `|c|_1 d^m` in modulus, and `sum_i binom(m,i) |eta|^i <= 2^m |eta|^o`; dividing by
`|y + eta v|_inf^m >= 2^{-m}` gives the bound. QED

(This is `zetc.md` Lemma 1.1(7) and `t1rigour.md` Props 5.1-5.2, re-proved; for `rho = sum h_A(v_i)` it is R11 §5.3's
"transversal order `m - d_F`". Every linear subspace `L` contained in `TC_H` has transversal order `o >= 1` in the same
sense, by Taylor's formula, since all `Phi_c` vanish on `L`.)

**Proposition 5.1 (what (MIN) forces).** Let `L` be a linear subspace of `R^d` contained in `TC_H(A, m)`, and suppose all
`Phi_c`, `c in V_H`, vanish to order `>= o` transversally along `L` as in Lemma 5.0 (`o >= 1` always).
* (i) If `P(L) = 0`, then `Lambda_A(g) = 0` for all `g >= g_L := min{Gamma(u) : u in L, |u|_inf = 1} < infinity`. So (MIN)
  fails at every level `abar >= g_L`.
* (ii) If `P(L) >= 1`, then there is `Gamma_L < infinity` (depending on `L` only) such that for every `g >= Gamma_L`,

  ```text
  ln(1/Lambda_A(g))  >=  o (N_2/(2 P(L))) (g - Gamma_L)  -  o ln(4n)  -  m ln 4 .
  ```

  Consequently, if MIN\*(A) holds (for all `g >= 1.5 ln n - 7`), then **`o N_2 / P(L) <= d`**; in particular every linear
  subspace of `TC_H` has `P(L) >= N_2/d = n/2` collision pairs, and every projection-cheap span (Lemma 5.0) has
  `P(L) >= (m - rho) n/2`.

*Proof.* (i) `L` is not covered by the finitely many proper subspaces `L n {u_x = u_y}`, so it contains points off `Br`;
for them `Gamma < infinity` and `kappa_H = 0`. (ii) Choose `y in L`, `|y|_inf = 1`, with `y_x != y_{x'}` for every pair that is
not a collision pair of `L` (generic `y`), and let `delta_y > 0` be the least such difference. Choose `v` with `|v|_inf <= 1`,
`v_n = 0`, and, inside each class of points on which `L` is constant, pairwise distinct values from `{0, 1/n, ..., (n-1)/n}`
(with `0` for the class of `n` taken by `n` itself). For `0 < eta <= delta_y/4` and `u = y + eta v`: collision pairs of `L` have
`|u_x - u_{x'}| >= eta/n`, the others `>= delta_y/2`, and `|u|_inf <= 3/2`. Hence
`Gamma(hat u) <= (2P(L)/N_2) ln(3n/(2 eta)) + 2 ln^+(3/delta_y)` (normalising by `|u|_inf in [1/2, 3/2]` costs the factor
`3/2`). Put `a_y = 2 ln^+(3/delta_y)`, `Gamma_L = a_y + 2 ln(6n/delta_y)`, and `eta = (3n/2) exp(-(N_2/(2P(L)))(g - a_y))`; since
`N_2/(2P(L)) >= 1/2`, for `g >= Gamma_L` this `eta` is `<= delta_y/4 <= 1/2`, and `Gamma(hat u) <= g`. Lemma 5.0 gives `Lambda_A(g) <= kappa_H(hat u) <= 4^m eta^o`, which is the display
(with `ln(3n/2) <= ln(4n)` and `g - a_y >= g - Gamma_L`). If `o N_2/P(L) > d`, the right side of the display
exceeds `(d/2) g + ln((m'+1)/d)` for all large `g`, contradicting MIN\*. QED

So (MIN) contains a **quantitative collision statement about every stratum of the tangent cone**: the strata must have at
least linearly many collision pairs, and more if they are tangent to higher order. This is the precise version of
`t1rigour.md`'s open question (Q1), and it explains Prop 2.7 (one collision pair, order `m - 2`). It is a property of the
profile's slab box alone; it can be tested on profiles without reference to any cluster.

### 5.2 The order must be small: (MIN) is false at large orders on Sylvester profiles

**Proposition 5.2.** Let `n = 2^k` and take the Sylvester-Hadamard profile (points `x in F_2^k`, columns `1[a.x = 0]`,
`a != 0`, `r = n - 1`, all `e_p = 1`; R11 §5.1). For every order `m''` with `2jk < m'' <= m_max`, the span `L` of the `k`
columns `1[x_i = 0]` lies in `TC_H(A, m'')`, has `P(L) = 0`, and contains a point `hat u` with `Gamma(hat u) <= 3`. Hence
`Lambda_{A,m''}(g) = 0` for all `g >= 3`, and (MIN) fails at every such order for every cluster with this profile (whose
`abar >= 1.5 ln n - 7 > 3` when `C <= 1` and `k >= 10`). The window `(2jk, m_max]` is nonempty exactly from `k = 16`
(exact: `2jk = 40896 < m_max = 42065` at `k = 16`, `26340 > 21034` at `k = 15`), while `ceil(4w) <= 2jk` for every `k >= 15`.

*Proof.* Each column `1[x_i = 0]` is a sign column with `e_p = 1`, so its width on `A` is `<= 2j`; the `k` columns cost `2jk < m''`,
so `L` lies in `TC_H` ((I2) Prop 2.2). The map `x -> (1[x_i = 0])_i` is injective, so `P(L) = 0`. Let `x_0` be the base point
and `y_x = sum_i 2^{-i} 1[x_i != x_{0,i}]` (in `L` modulo constants); its values are `t/n`, `t = 0, ..., n - 1`, each taken once,
with `y_{x_0} = 0`. Then `hat u = y n/(n-1)` has `Gamma(hat u) = (2/N_2) sum_{t=1}^{n-1} (n - t) ln((n-1)/t)`, which is `<= 3`
(computed: `2.9996` at `k = 16`; the sum tends to 3). The numbers are in `sylvester_gamma_out.txt` (exact big-integer `m_max`). QED

So the reduction of Theorem 4.1 can work only at small orders; at the recommended `m' = ceil(4w)` the obstruction of
Proposition 5.2 is absent (Theorem 3.4(d): every cheap family has at most 16 core classes). This is the quantitative
form of R11's remark that the ZE-braid counterexample "does not apply at this order".

### 5.3 Cheap spans at m' = ceil(4w) pass the necessary condition

**Proposition 5.3.** Let `m' = ceil(4w)`, `w >= 10`, `2j <= n`, and let `F` be a cheap family. Then
`m' N_2 / P(L_F) <= d` (so the necessary condition of Prop 5.1(ii) holds for `L_F` **whatever its true transversal order**,
since `o <= m'` along any line not in `TC_H`) as soon as

```text
max_{0 <= s <= b, 0 <= e <= floor((m'-1)/2) - js}   m' * 2^{s+1} N_2 / ( (n - e)(n - e - 2^s) )   <=   d .
```

*Scan* (`compat.py`, `compat_fine.py`; integer `j, m', b, e`, ratios in floating point with an exact rational recheck at the
threshold): the displayed condition holds for every `n >= n_A` in the scanned range `n <= 4 * 10^6`, with

| `r` | `n_A` (order-free condition) | `n_B` (proved orders only) |
| --- | --- | --- |
| `n` | 473,716 | 5,185 |
| `n ln n` | 564,084 | 6,398 |
| `n^2` | 947,350 | 10,280 |

Here `n_B` is the threshold for the weaker condition with the *proved* orders,
`max (m' - 2js - 2e) 2^{s+1} N_2/((n-e)(n-e-2^s)) <= d` (what `t1rigour.md` Prop 5.5 tests; its value `5185` is
reproduced). Asymptotically the two conditions read `n > (8 + o(1)) w` and `n > (64 + o(1)) w`; the binding families have
`s = 4` (16 classes) for the order-free condition and `s = 2, 3` for the other.

*Proof.* Theorem 3.4(c)-(d): `P(L_F) >= n'(n'-q_F)/(2q_F)` with `q_F <= 2^s`, `s = sum ell_i <= b`, `n' >= n - e`,
`e <= floor((m'-1)/2) - js`; the right side is decreasing in `q_F` and in `e`. QED

**Reading.** Cheap spans are harmless for (MIN) at the level of necessary conditions once `n >= n_A`; Proposition 5.2
shows that this fails at large orders; and the non-cheap rational strata (projection-cheap spans and anything else in
`TC_H`) must be checked separately against Prop 5.1. For the projection-cheap spans `R^{S'}` of `t1rigour.md` Thm 6.11,
`P >= binom(n - |S'| - 1, 2)` and `o <= m'`, so the condition holds when `m' n(n-1) <= d (n - |S'| - 1)(n - |S'| - 2)`, e.g.
when `|S'| <= n/2` and `n >= 4m' + 3`; but `|S'|` can be as large as `m' - j` there.

## 6. Computations

All scripts are in `nearstrata_checks/`; each writes `<name>_out.txt`. None is used in a proof except where a statement
says "scan" or "exact", and then only for the finite range stated.

| script | what | status |
| --- | --- | --- |
| `one_cut_gadgets.py` | Theorem 2.5: gadget zonotope inside `A_{S,L}`, degree `= nu_eta`, all generators `eta`-good, product identity `sum c_k <k',u>^m = m! prod <v,u>^alpha`, `|c|_1 <= 2^m`, part (b) strata; 3000 random instances, `d <= 7`, all four cases | exact (Fractions/integers); the moment check uses random powers of linear forms |
| `onecut_lp_check.py` | Theorem 2.5(a) against LP values of `kappa_V` (977 cases with `nu_eta >= m`; min ratio `kappa_V/bound = 11.95`) | evidence (floating LP) |
| `duality_check.py` | Prop 2.1 on 12 small cut boxes: primal (best uniform approximation) = dual (sup over `V`), relative gap `<= 1.1e-14` | evidence (floating LP) |
| `sylvester_gamma.py` | Prop 5.2: `Gamma` of the binary-fraction point (`<= 2.9996` at `k = 16`), `j`, `ceil(4w)`, `2jk`, `m_max` at `k = 15, 16, 17` | exact big-integer `m_max`; `Gamma` by 30-digit summation |
| `compat.py`, `compat_fine.py` | Prop 5.3 thresholds `n_A`, `n_B` for `r = n, n ln n, n^2` | integer parameters, float ratios, exact rational recheck at `n_A` |
| `sep_vs_radius.py` | Corollary 3.7 thresholds (`5185 / 6425 / 10280`) | integer parameters, float margins (smallest margin beyond the threshold `> 2000`) |

## 7. Open points

1. **MIN\* for regime profiles at `m' = ceil(4w)`** (§4.3). This is the whole remaining content of route 5 at `T = 1` along
   the pair-inequality path. Its qualitative core is ZE-TC^braid (not ZE-TC^rat), and its shape is fixed: an inequality
   `Lambda_A(g) >= gamma e^{-theta g/2}` with `theta < d` and `ln(1/gamma) <= (3/8 - o(1)) n ln n` suffices (Theorem 4.3(a)),
   while general degree/height arguments provably cannot give it (Prop 2.7).
2. **Collision counts of all strata (necessary).** Prop 5.1: MIN\* forces `P(L) >= o_L n/2` for every linear subspace `L` of
   `TC_H` (transversal order `o_L`). Cheap spans pass from `n_A ~ 4.7 * 10^5` whatever their orders, and from `n_B = 5185`
   with the proved orders (`r = n`). Open: (a) are there rational subspaces of `TC_H(A, ceil(4w))` with few collision pairs
   (fewer than `~ m' n/2`)? (b) a lower bound for the projection cost `rho_A(F)` in terms of the collision classes of `F`,
   which would extend Theorem 3.4 to projection-cheap spans (`t1rigour.md` (Q1)); (c) an *upper* bound on the true
   transversal order of cheap spans better than `m'` (it would lower `n_A` towards `n_B`); product certificates cannot
   give it (Remark 2.10).
3. **A non-zonotopal model with an effective bound.** Theorems 2.4-2.5 are the only effective Lojasiewicz inequalities I can
   prove, and their boxes are zonotopal. The faithful Hadamard models of `compute.md` (`M8'`: `n = 8`, `T = 1`, `L = 3`, `m = 11`,
   seven explicit basis forms with ZE-TC^braid proved by saturation) are the natural first test of an explicit
   Nullstellensatz certificate `prod (collision forms)^e = sum g_c Phi_c`, which by the argument of Prop 2.6 would give an
   explicit `CL`. At `n = 8` they cannot satisfy MIN\* (Prop 5.1 needs `P(L) >= n/2` with `o_L N_2/P(L) <= 7`), so this would be
   a test of mechanism only.
4. **More arithmetic.** (MIN) is minimal only relative to the averaged pair data (§4.3). Strata with few collision pairs
   could still be handled by character inequalities for relations `v` orthogonal to `L` that are not pair differences
   (`X_v` small), or by the algebraicity of `e^{i delta}` (ZE-alg). Nothing in this write-up uses them.
5. **Constants.** The thresholds (`n >= 2900` for (R1), `5185` for Cor 3.7, `4.7 * 10^5` for the order-free Prop 5.3, all at
   `r = n`) come from Hoeffding plus a union bound and from Theorem 3.4's crude class count; none is optimised.
