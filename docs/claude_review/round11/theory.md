# Round 11 -- theory: the zero-exclusion hypothesis ZE for the slab-box exact system

## 0. Summary of results

**ZE is not proved, and no growth improvement is proved.** What is proved:

1. **`V_H` is almost all of `V(A,m)`** (Lemmas 1.1-1.2). The dichotomy of round-10 Theorem 4.1 holds on *any*
   support `S` in `B_T` whose cost is `<= 3mW/8`, for heights up to `H = e^{mW/40}`; with that `H`,
   `codim(V_H in V(S,m)) <= (13.8 + o(1)) E_S/n + 27.8|S|/(mW)` in the growth regime. This answers the referee's worry that `Z_P`
   might be cut out by only half of `V(A,m)`.
2. **Angle (a), the tangent cone** (§2). The tangent cone of `Z_P` at `e` lies in the common real zeros of the
   degree-`m` forms `Phi_c` (the degree-`m` least space of `A`; Prop 2.1). Its toric part (spans of *cheap*
   families, `sum_j h_A(v_j) < m`) is exactly degenerate (Prop 2.2). Cube slices reduce to parallelepipeds, i.e. to
   characters (Prop 2.3). Lifting from the tangent cone is trivial because `|delta| <= C e^{-W/4}` (Lemma 2.4). **But
   the leading-order (tangent-cone) criterion cannot decide ZE near low-rank strata** (Prop 2.5): `Phi_c` is
   `O(eta)` at relative distance `eta` from a toric stratum, the criterion needs `eta >~ e^{-W/4} n^{-3/2}`, and the gap
   lemma only guarantees `eta >~ e^{(W/2 - X_v)/2}` for relations `v` orthogonal to the stratum; this is useless when
   all of them have `X_v > W`.
3. **A delta-free substitute for ZE** (§3, the main new result). Widening the slab of **one** light prime `p`
   (`log p <= W/8`) to `|<k, beta_p>| <= m/e_p` costs at most `mW/8` in the LP, which the slack of Theorem 2.1(d)
   pays for. On this widened box `A_p` the level obstruction of round-10 Prop 5.5(ii) disappears, and:
   * **Theorem 3.1 (unconditional).** For every cluster in the regime (`n - 1 >= 300 ln 4r`,
     `W >= max(27 log^+ C + 40, 80 ln n + 9)`) and every split `p` with `log p <= W/8`, no vector of `V_0(A_p)` has a
     unique `p`-adic minimiser of `v_p(c_k) + <k, a_p>`. This is a necessary condition on the **profile alone**: the
     positions `delta` do not enter.
   * **Theorem 3.3 (conditional).** Under Hypothesis SL(n) (every admissible full-rank, pair-separated profile has
     such a vector at eight or more columns), clusters on squarefree norms satisfy `n = O(log W) = O(loglog R)`: the
     conclusion of round-10 Theorem 5.7, with ZE replaced by SL.
   * SL needs no nondegeneracy, no Diophantine input and no residues; for a given profile it is a finite computation
     (Prop 3.4: equivalent to the absence of a `p`-adic separator polynomial; in corner form, to a torus-fixed point not
     being a base point of the fat-point system). Prop 3.5 lists the obstructions: too few `a_p`-levels, and `a'_p` in the
     rational span of a cheap family (duplicated or nearly sparse columns, low rank). They are the `p`-adic shadows of
     the toric strata, and pair separation keeps the weight on columns they do not obviously hit.
4. **Angle (c), arithmetic** (§4). The actual point `w in Q(i)^{n-1}` has `pi`-adic valuation vector exactly `+-a'_p` at
   every split place (Prop 4.1). Theorem 3.1 exploits one such place; the remaining arithmetic tools (arithmetic
   Bezout, multiplicity estimates) are still exponentially short.
5. **Angle (b)** gave only the reformulation in Remark 3.6 (ZE and SL are both "no non-toric base locus", near `e` and
   at infinity in direction `a'_p`), not a positivity theorem.

**Evidence** (§5; small models only): in every small instance tested, SL held iff a product (parallelepiped)
certificate existed, and it failed once products were excluded. These sizes are far below the supercritical regime
(`n >= 2137`), so they neither support nor refute SL there; they show that SL is not automatic.

| item | status |
| --- | --- |
| ZE (round 10) | open; the leading-order tangent-cone criterion is provably inconclusive near low-rank strata (Prop 2.5) |
| `V_H` vs `V(A,m)` | proved: relative codimension `O(1/n)` with `H = e^{mW/40}` (Lemma 1.1) |
| Theorem 3.1: clusters force failure of SL at every light prime | **proved** |
| Theorem 3.3: SL(n) implies `n = O(loglog R)` on squarefree norms | **proved (conditional)** |
| Hypothesis SL(n) | open; delta-free; decidable per profile; small cases decided by products |
| growth improvement for actual circles | **not proved** |

## 1. Setup and notation (recalled from round10/determinant.md)

Labels, as in earlier rounds: *Theorem / Proposition / Lemma* = proved here with all quantifiers, or quoted with
an exact source; *Exact* = integer or exact-rational computation; *Evidence* = finite or floating-point experiment
that no proof uses; *Hypothesis* = unproved. Scripts are in `round11/theory_checks/` (outputs `*_out.txt`).

**Cluster data.** `N >= 2`, `W = log N`, `R = sqrt N`; a cluster is `n` distinct `z_x in Z[i]`, `|z_x|^2 = N`, on an arc
of length `<= C sqrt R`; `delta_x` = angle of `z_x` minus the arc midpoint, `|delta_x| <= Theta/2`,
`Theta <= C N^{-1/4}`. For a split `p = pi conj(pi)`, `p^{e_p} || N`: `a_{p,x} = v_pi(z_x)`, `beta_p = (2a_p - e_p 1)/e_p`.
On `Z^n_0 = {k : sum k = 0}`: `<k, a_p> = (e_p/2) <k, beta_p>`. `r` = number of split primes.

**Slab box (round 10, Thm 2.1).** `n >= 64`, `T = n`, `w = T sqrt(8(n-1) ln 4r)`, `D = 2(n-1)T`,
`B_T = {k in Z^n_0 : |k_x| <= T, x < n}`, `A = {k in B_T : |<k,beta_q>| <= w for all q}`, `m` = largest integer with
`binom(m+n-2, n-1) <= (2T+1)^{n-1}/8`; `m >= 0.7 n(n-1)`; if `n - 1 >= 300 ln 4r` then `m >= 4w`.
`V(S, m)` = vectors `c in Q^S` with `sum_k c_k P(k) = 0` for all polynomials `P` of degree `< m`; `V_H(S)` = the
`Q`-span of the integer `c in V(S,m)` with `|c|_inf <= H`. `f_delta(k) = e^{i<k,delta>}`,
`F_c(z) = N^{D/2} sum_k c_k e^{i<k,delta>}` (a Gaussian integer when `|k|_1 <= D` on `supp c`).
`Z_P = {delta : sum_k c_k e^{i<k,delta>} = 0 for all c in V_H(A)}`.

Coordinates. Using `k' = (k_1, ..., k_{n-1})` (and `k_n = -sum k'`), `<k, delta> = <k', u>` with `u_x = delta_x - delta_n`,
and `<k, beta_q> = <k', b_q>` with `b_{q,x} = beta_{q,x} - beta_{q,n}`. A point `delta` is *nondegenerate* iff
`<v, u> != 0` for every nonzero `v in Z^{n-1}`. `w = e^{iu}` denotes the torus point; `e = (1, ..., 1)`.
`h_S(v) = max_S <k,v> - min_S <k,v>` (width). `E = binom(m+n-2, n-1) = dim Poly_{<m}` in `n-1` variables,
`E_S = dim Poly_{<m}|_S <= E`, `U_S = |S|`.

**Imported (refereed).**
* (R10-4.1) `round10/determinant.md` Thm 4.1 and its proof. The proof shows, for *any* finite `S` in `B_T` and any
  integer `c in V(S, m)`, `c != 0`, `|k|_1 <= D` on `S`: if `F_c(z) != 0` then

  ```text
  cost(c) := sum_{q split} log q * h_{supp c}(a_q)  >=  mW/2 - 2 ln|c|_1 - 2m ln(DC/2) + 2 ln m!.           (1.1)
  ```

  (Upper bound `|F_c| <= N^{D/2} |c|_1 (D Theta/2)^m / m!` from Taylor; lower bound
  `Norm F_c >= N^D prod_q q^{-h_{supp c}(a_q)}` from forced divisibility; compare.) Only `m`, `D`, `C` and the
  support enter; the slab structure is used only to bound `cost(c)`.
* (R10-5.5) `round10/determinant.md` Prop 5.5(i): for a split `p` and integer `c`, if the minimum of
  `v_p(c_k) + <k, a_p>` over `supp c` is attained at exactly one `k`, then `F_c(z) != 0` (indeed
  `v_pi(F_c(z)) = e_p D/2 + min`). More generally `F_c(z) != 0` whenever the face residue sum is nonzero mod `pi`.
  Prop 5.5(ii): this never happens when `h_{supp c}(a_p) < m` (level-sum lemma 5.4).
* (R10-3.1/3.2) Lemma 3.1 (squarefree clusters have exponent rank `>= n/4 - 2`) and Lemma 3.2 (a full-rank
  subcluster with affinely independent exponent vectors exists, and its positions are nondegenerate), with the
  `(pi/4)Z` gap of `outside.md` Lemma 0.1 as corrected by the referee.
* (Ref-R1) `referee_determinant.md` Prop R1: for the uncut box, `V(B_T, m)` is spanned by the forward-difference
  vectors `Delta^alpha` (`|alpha| >= m`), and the base locus is the sparse strata.

### 1.1 Two technical upgrades of Theorem 2.1 / 4.1

**Lemma 1.1 (Siegel with large height; codimension of `V_H`).** Let `S` be finite, `U = |S|`, `E_S = rank` of the
moment matrix of `S` (rows `k'^alpha`, `|alpha| < m`, entries bounded by `B = T^{m-1}` on `B_T`). For every `H >= 1`,

```text
dim V(S, m) - dim V_H(S)  <=  (U ln 2 + E_S ln(UB)) / ln(2H + 1).
```

*Proof.* Keep `E_S` independent rows `M_0` of the moment matrix (same kernel `V(S,m)`). The map `x -> M_0 x` sends
`{0..H}^U` into `[-UHB, UHB]^{E_S}`, so some fibre has `>= (H+1)^U / (2UHB+1)^{E_S}` points; their differences with a
fixed point of the fibre are distinct integer kernel vectors with `|.|_inf <= H`, all in `V_H(S)`. A subspace of
dimension `d` contains at most `(2H+1)^d` integer points of `[-H,H]^U` (it projects injectively to some `d`
coordinates). With `L = ln(2H+1)`, `ln(H+1) >= L - ln 2` and `ln(2UHB+1) <= L + ln(UB)`:
`d >= U - E_S - (U ln 2 + E_S ln(UB))/L`, and `dim V(S,m) = U - E_S`. QED

So `V_H` is *almost all* of `V(S, m)` once `ln H` is large compared with `ln(UB)`. For `S` in `B_T` and
`n >= 64`, `ln(UB) <= (n-1) ln(2n+1) + (m-1) ln n <= 1.03 m ln n`. Lemma 1.2 below allows `ln H = mW/40`; then

```text
codim( V_H(S) in V(S,m) )  <=  27.8 U/(mW) + 41.2 E_S ln n / W  <=  (13.8 + o(1)) E_S / n + 27.8 U/(mW),
```

the last step for `W >= (3 - o(1)) n ln n`, which holds for every cluster in the growth regime
(`two_thirds` gives `W >= 3n log n - O(n)`). So the base locus `Z_P` of
the exact system is cut out by all of `V(A,m)` except for at most `~14 E_S/n` conditions, a relative codimension
`O(1/n)`. This does not decide ZE, but it answers the referee's worry (`referee_determinant.md` §5) that `V_H` might
be only a half-dimensional subspace of `V(A,m)`: after enlarging `H` it is not.

**Lemma 1.2 (general dichotomy on any support).** Let `n >= 64`, let the cluster satisfy `W >= 27 log^+ C + 40`, let
`S` be a finite subset of `B_T` (so `|k|_1 <= D` on `S`), `m >= 0.7 n(n-1)` an integer and `1 <= H <= e^{mW/40}`. If
every integer `c in V(S, m)` with `|c|_inf <= H` has `cost(c) <= 3mW/8`, then `F_c(z) = 0`, i.e.
`sum_k c_k e^{i<k,delta>} = 0`, for every `c in V_H(S)` (the whole `Q`-span, with no bound on heights).

*Proof.* For a generator `c` (integer, `|c|_inf <= H`): `2 ln|c|_1 <= 2 ln(|S| H) <= 2(n-1) ln(2n+1) + mW/20`. By
`ln m! >= m ln(m/e)` and `De/(2m) <= e/0.7 < 3.884`, the right side of (1.1) is at least
`mW/2 - mW/20 - 2(n-1) ln(2n+1) - 2m ln(3.884) - 2m log^+ C`. It exceeds `3mW/8` as soon as
`W > (80/3) ((n-1) ln(2n+1)/m + ln 3.884 + log^+ C)`; since `(n-1) ln(2n+1)/m <= ln(2n+1)/(0.7 n) <= 0.109` for
`n >= 64`, the hypothesis `W >= 27 log^+ C + 40` suffices. So (1.1) fails and `F_c(z) = 0`. Extend by linearity to
the span. QED

(Round 10's Theorem 4.1 is the case `S = A`, `H = (|A| n^{m-1})^2`, `cost <= wW <= mW/4`. The point of Lemma 1.2 is
that the support need not be the slab box: any `S` whose cost stays below `3mW/8` works, and heights up to
`e^{mW/40}` are allowed. Any fixed fraction `< 1/2` in place of `3/8` works with other constants.)

## 2. Angle (a): the tangent cone, its toric part, and lifting

Throughout this section `S` is any finite subset of `Z^n_0` (written in the coordinates `k'`), `m >= 1`, and
`Phi_c(u) = sum_{k in S} c_k <k', u>^m` for `c in V(S, m)`, `u in R^{n-1}`. Put

```text
TC(S)   = {u in R^{n-1} : Phi_c(u) = 0 for all c in V(S, m)},
TC_H(S) = {u in R^{n-1} : Phi_c(u) = 0 for all c in V_H(S)}   (contains TC(S)).
```

**Proposition 2.1 (what the tangent cone is).**
* (i) `u in TC(S)` iff the function `k -> <k', u>^m` on `S` coincides on `S` with a polynomial of degree `< m`.
* (ii) For `c in V(S,m)` and real `t`, `sum_k c_k e^{i t <k',u>} = (i t)^m Phi_c(u)/m! + O(t^{m+1})`. Hence the
  tangent cone at `0` of the real-analytic set `Z_P` (and of the base locus) is contained in `TC_H(A)`.
* (iii) In the language of interpolation theory: `{Phi_c : c in V(S,m)}` is the degree-`m` part of the de Boor-Ron
  least space of `S`, and (i) says `<k',u>^m` lies in the degree-`m` part of the ideal of leading forms of the
  vanishing ideal `I(S)`.

*Proof.* (i) `V(S,m)` is the annihilator of `Poly_{<m}|_S` in `Q^S`, so a vector of `C^S` is orthogonal to all of
`V(S,m)` iff it lies in `Poly_{<m}|_S (x) C`. Apply this to `(<k',u>^m)_k`. (ii) Expand `e^{ix}`; the terms of
order `j < m` vanish because `<k',u>^j` has degree `j < m`. (iii) (i) restated: `<k',u>^m - g` vanishes on `S` for
some `g` of degree `< m` iff `<k',u>^m` is the top form of an element of `I(S)`; the least-space statement is the
dual description (de Boor-Ron), and is not used below. QED

**Proposition 2.2 (the toric part of the base locus is degenerate).** Call a family `v_1, ..., v_s in Z^{n-1}` (as
linear forms `k' -> <k', v_j>`) *cheap for* `(S, m)` if `sum_j h_S(v_j) < m`.
* (i) If `(v_j)` is cheap, then for every `u in span_R(v_1, ..., v_s)` the function `k -> e^{i<k',u>}` coincides on `S`
  with a polynomial of degree `< m`. So `span_R(v_j)` lies in the base locus of `V(S,m)` (hence in `Z_P` when
  `S = A`), and in `TC(S)`.
* (ii) If `V(S, m) != 0`, every cheap family has rank `<= n-2`. Consequently every point of `span_R(v_j)` is
  degenerate: there is a nonzero `v in Z^{n-1}` with `<v, u> = 0` for all `u` in the span.
* (iii) Examples (`S = A`, sign profile): `{e_x : x in X}` with `2T|X| < m` (the sparse strata of round-10 Prop 5.1),
  `{b_q : q in Q}` with `2w|Q| < m` (low-rank strata, Prop 5.2), and mixtures with `2T|X| + 2w|Q| < m`.

*Proof.* (i) `y_j = <k', v_j>` takes at most `h_S(v_j) + 1` integer values on `S`. A function of `y = (y_j)` on the
grid `prod_j {values of y_j}` is a polynomial in `y` of degree `<= sum_j h_S(v_j)` (tensor Lagrange interpolation),
and `y` is linear in `k'`. For `u = sum_j t_j v_j`, `e^{i<k',u>} = prod_j e^{i t_j y_j}` is a function of `y`.
(ii) If some `n-1` of the `v_j` were independent, `k' -> y` would be injective on `Z^{n-1}`, so every function on `S`
would be a function of `y`, hence a polynomial of degree `< m`, i.e. `Poly_{<m}|_S = Q^S` and `V(S,m) = 0`. So the
`v_j` span a rational subspace of dimension `<= n-2`, whose orthogonal complement contains a nonzero integer `v`.
(iii) `h_A(e_x) <= 2T` and `h_A(b_q) <= 2w` by construction. QED

So the base locus contains a union `Z_tor` of rational linear subspaces through `0`, all degenerate; ZE-real is the
statement that, near `0`, nothing else meets the nondegenerate set. (Translates `g + span` with `g` not in the span,
and non-linear components, are the unknown parts; README §8.3.)

**Proposition 2.3 (cube slices: the tangent-cone analogue of round-10 Prop 5.3).** Let `k_0 in Z^{n-1}`,
`v_1, ..., v_M in Z^{n-1}`, `Psi(t) = k_0 + sum_i t_i v_i`, and `X` a subset of `{0,1}^M` with `Psi(X)` contained in `S`.
For `u in R^{n-1}` put `alpha_i = <v_i, u>`.
* (i) If `u in TC(S)`, then `(<k_0,u> + <alpha, t>)^m` coincides on `X` with a polynomial in `t` of degree `< m`.
* (ii) If moreover `X` contains a subcube `{t : supp t in J}` with `|J| >= m`, then fewer than `m` of the `alpha_i`,
  `i in J`, are nonzero. In particular, for nondegenerate `u` (all `alpha_i != 0`), a subcube of dimension `>= m` in a
  slice certifies `u notin TC(S)`.

*Proof.* (i) Pull back the identity of Prop 2.1(i) along `Psi`; polynomials of degree `< m` in `k'` pull back to
polynomials of degree `< m` in `t`. (ii) On the subcube, functions have a unique multilinear representative of degree
`<= |J|`, obtained by `t_i^2 -> t_i`; this does not raise degree. The multilinear representative of
`(c_0 + <alpha,t>)^m` has coefficient `m! prod_{i in J'} alpha_i` at `t^{J'}` for every `J'` in `J` with `|J'| = m`, and
the representative of a polynomial of degree `< m` has no such terms. QED

A subcube of dimension `>= m` in a slice of `A` is a lattice parallelepiped with `>= m` unit sides, so (ii) is
round-10 Prop 5.3 at the level of the tangent cone, and it is capped in the same way: it needs `>= m` directions `v_i`
with total width `sum_i |<v_i, b_q>| <= 2w` in every slab. For a sign profile `<v, b_q>` is an even integer, so each
slab can be "spent" by at most `w` of the `v_i`. *Square-free reformulation (remark).* Without the subcube, (i) asks
whether the degree-`m` square-free form `e_m(alpha o t) = sum_{|J|=m} alpha^J t^J` lies in the degree-`m` part of the
leading-form ideal of `I(X)` in `C[t]/(t_i^2)`. Since `I(X)` is generated (in the function ring of the cube) by the slab
polynomials `prod_{v}(<t, u_q> - v)` (`u_q = (<v_i, b_q>)_i`, `v` over the allowed slab values), this ideal contains
`(t_i^2, <t,u_q>^{d_q})`, `d_q` = number of slab values; whether it is generated by them in degree `m` is a
Froeberg-type question for powers of linear forms in a square-free algebra. Rank-one elements of the inverse system
give back exactly the subcubes of (ii). I did not get further with this.

**Lemma 2.4 (lifting from the tangent cone; quantitative).** Let `c in V(S, m)`, `u in R^{n-1}`, `t in R`, and
`D_u = max_{k in S} |<k', u>|`. Then

```text
| sum_k c_k e^{i t <k',u>}  -  (it)^m Phi_c(u)/m! |  <=  |c|_1 (|t| D_u)^{m+1} / (m+1)!.
```

Hence if `|Phi_c(u)| > |c|_1 |t| D_u^{m+1}/(m+1)`, then `t u notin Z(c)`; in particular `tu notin Z_P` if `S = A` and
`c in V_H(A)`. For the actual cluster (`S = A`, `D_u <= (D/2)|u|_inf`), writing `u = |u|_inf * u_hat`, a sufficient
condition for `delta notin Z_P` is

```text
kappa(u_hat) := sup_{c in V_H(A)}  |Phi_c(u_hat)| / (|c|_1 (D/2)^m)   >   (D/2) |u|_inf / (m+1).                (2.1)
```

*Proof.* `|e^{ix} - sum_{j <= m} (ix)^j/j!| <= |x|^{m+1}/(m+1)!`; the terms `j < m` cancel against `c in V(S,m)`.
For (2.1) apply this with `t = |u|_inf`, `u = u_hat`, `D_{u_hat} <= (D/2)` (as `|k'|_1 <= (n-1)T = D/2`). QED

The right side of (2.1) is tiny: `|u|_inf <= C e^{-W/4}`. So (2.1) is a *very* weak requirement, a relative lower bound
`e^{-W/4} poly(n, C)` on the leading form. The next proposition shows where even this fails.

**Proposition 2.5 (where the tangent cone cannot decide).** Let `L` be a linear subspace of `R^{n-1}` contained in
`TC_H(A)` (for example a toric stratum of Prop 2.2), let `u_hat` satisfy `|u_hat|_inf = 1`, and let `eta` be its
Euclidean distance to `L`. Then

```text
kappa(u_hat)  <=  m eta / sqrt(n-1).
```

Consequently (2.1) can hold only if `eta > sqrt(n-1) (D/2) |u|_inf / (m(m+1)) >= |u|_inf / (10 n^{3/2})` (`n >= 64`).
For the actual cluster with nondegenerate `delta`, every nonzero integer `v` orthogonal to `L` gives
`eta >= sqrt 2 e^{-X_v/2} / (|v|_2 |u|_inf)`, where `X_v = sum_q log q |<v~, a_q>|` and `v~ = (v, -sum v) in Z^n_0`.
This lower bound exceeds `|u|_inf/(10 n^{3/2})` whenever `X_v < W + 2 ln(14 n^{3/2} / (C^2 |v|_2))` (recall
`|u|_inf <= C e^{-W/4}`), and, when `|u|_inf` has its maximal order `C e^{-W/4}`, essentially only then: roughly,
the available arithmetic separates `delta` from `L` at the scale needed by (2.1) iff some `v` orthogonal to `L` has
`X_v < W`.
* Sparse strata `L = span{e_x : x in X}`: take `v = e_y`, `y notin X`, i.e. the pair `(y, n)`; `X_v = d_{yn} <= W`. The
  available arithmetic is borderline, up to `O(log(nC))`.
* Low-rank strata: the vectors orthogonal to `L` are relations among several columns, and I know no argument
  excluding `X_v > W + O(log n)` for all of them. At a `delta` that is nondegenerate but within relative distance
  below `|u|_inf/(10 n^{3/2})` of such a stratum, (2.1) fails for every `c`, whatever the tangent cone is: the
  leading-order criterion is inconclusive there, and membership in `Z_P` is decided by higher-order terms.

*Proof.* Let `u_L` be the orthogonal projection of `u_hat` on `L`; `Phi_c(u_L) = 0` for `c in V_H(A)`. Points `y` of the
segment `[u_L, u_hat]` satisfy `|y|_2 <= |u_hat|_2 <= sqrt(n-1)`, and `|k'|_2 <= T sqrt(n-1)` on `B_T`, so
`|<k',y>| <= T(n-1) = D/2`. Since `grad Phi_c(y) = m sum_k c_k <k',y>^{m-1} k'`,
`|grad Phi_c(y)|_2 <= m |c|_1 (D/2)^{m-1} T sqrt(n-1)`, and by the mean value theorem
`|Phi_c(u_hat)| <= m |c|_1 (D/2)^{m-1} T sqrt(n-1) eta`. Divide by `|c|_1 (D/2)^m = |c|_1 (D/2)^{m-1} T (n-1)`.
For the consequence: `m < (n-1)(2n+1) < 3n^2` (because `binom(m+n-2,n-1) >= (m/(n-1))^{n-1}` must be
`< (2n+1)^{n-1}`), so `m(m+1) <= 9.1 n^4` and `sqrt(n-1) n(n-1)/(9.1 n^4) >= 1/(10 n^{3/2})` for `n >= 64`. For the
arithmetic bound: `eta >= |<v, u_hat - u_L>|/|v|_2 = |<v, u>|/(|v|_2 |u|_inf)`, `<v,u> = <v~, delta>`, and for
`<v~,delta> != 0` the corrected round-10 Prop 5.3(iii) gives `|<v~,delta>| >= 2 arcsin(e^{-X_v/2}/sqrt 2)
>= sqrt 2 e^{-X_v/2}`. Finally `d_{yn} = sum_q log q |a_{q,y} - a_{q,n}| <= sum_q e_q log q <= W`. QED

**Verdict on angle (a).** Proved: the tangent cone is the degree-`m` least space (2.1); its toric part is degenerate
(2.2); cube slices reduce to products (2.3); lifting is a triviality (2.4) because `|delta|` is so small. Not proved:
that `TC_H(A)` has no non-toric real points. And even a complete description of `TC_H(A)` would not give ZE, because
actual clusters may sit (in relative terms) exponentially close to low-rank strata (2.5). The tangent-cone route
therefore has to be combined with an analysis of the base locus *along* the strata, i.e. a multi-scale argument. I
did not find one. This is the reason for the change of angle in §3.

## 3. A delta-free substitute: valuation certificates on a widened slab box

### 3.1 The observation

Round-10 Prop 5.5 offers a nonvanishing certificate that does not involve `delta` at all: if the minimum of
`v_p(c_k) + <k, a_p>` over `supp c` is attained once, then `F_c(z) != 0` by an exact `pi`-adic valuation. It was
discarded because of Prop 5.5(ii): in the slab box every column has `h_A(a_p) <= e_p w < m`, and then the level-sum
lemma forbids a unique minimiser. But the LP of Theorem 4.1 has slack: its cost is `<= wW <= mW/4`, against a
budget of nearly `mW/2`. That slack pays for removing the slab of **one** light prime `p` (`log p <= W/8`), so that
`a_p` gets `>= m` levels on the support, while every other slab stays narrow. The level obstruction then disappears,
and Lemma 1.2 still forces `F_c(z) = 0`. So a single element of the corresponding linear system with a unique
`p`-adic minimiser contradicts the existence of the cluster. **The positions `delta` no longer enter**: the missing
statement becomes a property of the profile alone, decidable by linear algebra over `Z_(p)`.

### 3.2 Definitions

Fix a profile `P` on `n` points (columns `a_q in {0..e_q}^n`, `beta_q = (2a_q - e_q 1)/e_q`), a column `p`, and the
round-10 parameters `T = n`, `w`, `m`. Put `l_p(k) = <k, a_p>` (an integer on `Z^n_0`) and

```text
A_p = {k in B_T : |<k, beta_q>| <= w for every column q != p,   |<k, beta_p>| <= m/e_p},
H_0(S) = (|S| n^{m-1})^2 ,     V_0(S) := V_{H_0(S)}(S)   for S contained in A_p.
```

A vector `c in Q^S` has the **unique-minimiser property at `p`** (`UM_p`) if `min_{k in supp c} (v_p(c_k) + l_p(k))` is
attained at exactly one `k`.

> **Hypothesis SL(P, p) (separation at `p`).** `V_0(A_p)` contains a vector with `UM_p`.

`V_0(S)` grows with `S` (both the support and `H_0` grow), so `S = A_p` is the most generous choice. SL(P,p) depends
only on the profile; given the profile it is a finite problem of linear algebra over `Z_(p)` (Prop 3.4).

### 3.3 The theorem

**Theorem 3.1 (valuation route).** Let `n >= 64` and let `z_1, ..., z_n` be a cluster with arc constant `C`, `r` split
primes and profile `P`, such that

```text
n - 1 >= 300 ln(4r),     W >= 27 log^+ C + 40,     W >= 80 ln n + 9.                                   (3.1)
```

Let `p` be a split prime of `N` with `log p <= W/8`. Then:
* (i) `sum_{k} c_k e^{i<k, delta>} = 0` (equivalently `F_c(z) = 0`) for every `c in V_0(S)`, for every `S` contained in
  `A_p`;
* (ii) SL(P, p) fails. More precisely, no `c in V_0(A_p)` has `UM_p`, and for every integer `c in V_0(A_p)` the face
  residue sum of round-10 Prop 5.5(i) vanishes mod `pi`; symmetrically at `conj(pi)` with maximisers of
  `l_p(k) - v_p(c_k)`.

*Proof.* (i) Let `c in V(S, m)` be an integer vector with `|c|_inf <= H_0(S)`. For `q != p`,
`h_{supp c}(a_q) = (e_q/2) h_{supp c}(beta_q) <= e_q w`; for `p`, `h_{supp c}(a_p) <= (e_p/2)(2m/e_p) = m`. Hence
`cost(c) <= w sum_q e_q log q + m log p <= wW + mW/8 <= 3mW/8`, using `m >= 4w` (round-10 Thm 2.1(d), which needs
`n - 1 >= 300 ln 4r`). Also `ln H_0(S) = 2 ln|S| + 2(m-1) ln n <= 2(n-1) ln(2n+1) + 2(m-1) ln n <= mW/40` because
`W >= 80 ln n + 9` and `80 (n-1) ln(2n+1)/m <= 8.7` for `n >= 64`. Lemma 1.2 applies (`S` is in `B_T`, `m >= 0.7n(n-1)`,
`W >= 27 log^+ C + 40`) and gives (i).
(ii) If `c in V_0(A_p)` has `UM_p`, multiply it by a nonzero integer to make it integral; this adds `v_p(lambda)` to every
`v_p(c_k)`, so `UM_p` is kept. By round-10 Prop 5.5(i) (exact `pi`-adic valuation of `F_c(z)`, valid since `|k|_1 <= D` on
`B_T`), `F_c(z) != 0`, contradicting (i). The face-sum statement is round-10 Prop 5.5(iii) applied to (i). QED

**Corollary 3.2 (a new delta-free necessary condition).** If a cluster satisfies (3.1), then SL(P, p) fails for
every split prime `p` with `log p <= W/8`. At most seven split primes have `log p > W/8` (eight of them would give
`sum_q e_q log q > W`). So the profile of a cluster in the regime fails SL at all but at most seven of its columns.

**Theorem 3.3 (conditional growth, with SL in place of ZE).** Call a profile on `n` points *admissible* if its
exponent rank is `n - 1`, `n - 1 >= 300 ln(4r)`, and it is pair-separated:
`sum_q (log q / W) |a_{q,x} - a_{q,y}| >= 1/2 - 1/(2n)` for all `x != y`. Assume

> **Hypothesis SL(n).** Every admissible profile on `n` points satisfies SL(P, p) at eight or more columns `p`.

Let `n_1 >= 64` and assume SL(`n'`) for all `n' >= n_1`. Let `z_1, ..., z_n` be a normalised cluster on a squarefree
norm `N'` with arc constant `C'` and `W' = log N'`, with `W' >= 4n log^+ C'`, `W' >= 27 log^+ C' + 40` and
`W' >= 80 ln n + 9`. Then `n < 4 max(n_1, 300 ln(4r) + 1) + 4 <= 4 max(n_1, 300 ln(4W'/log 5) + 1) + 4`, so
`n = O(log W') = O(loglog R)`.

*Proof.* By round-10 Lemma 3.1 the exponent rank is `k >= n/4 - 2`; by Lemma 3.2 there is a subcluster `J` of
`n' = k + 1 >= n/4 - 1` points with affinely independent exponent vectors, so its profile `P_J` (the columns
restricted to `J`, the same `r` primes) has exponent rank `n' - 1`. `J` is a cluster on the same norm and arc, and
it is pair-separated: by (I4) of round 10 (`d_xy >= W'/2 - 2 log^+ C'`) and `W' >= 4n log^+ C'`, we get
`d_xy/W' >= 1/2 - 1/(2n) >= 1/2 - 1/(2n')`. Suppose `n' >= max(n_1, 300 ln(4r) + 1)`. Then `P_J` is admissible and
SL(`n'`) gives eight columns with SL; by Corollary 3.2 (applied to `J`, whose parameters satisfy (3.1) because
`n' <= n`) at most seven columns can satisfy SL. Contradiction. Hence `n/4 - 1 <= n' < max(n_1, 300 ln(4r) + 1)`,
and `r <= W'/log 5`. QED

**Comparison with ZE (round-10 Thm 5.7).** Theorem 3.3 has the same conclusion. Its hypothesis differs in kind:
* ZE quantifies over all nondegenerate *positions* `delta` near `0`; SL(n) quantifies only over *profiles*, and for
  each profile it asks for one vector in an explicit rational subspace;
* SL needs no nondegeneracy, no information about `delta`, and no Diophantine input;
* SL(P,p) is decidable for a given profile (a finite computation: enumerate the short kernel vectors spanning `V_0`,
  then a rank test over `Z_(p)` by Prop 3.4), while ZE for a given profile is a statement about a real-analytic set;
* the barrier is not an issue in the same way: Theorem 3.1 uses integrality and the exact `pi`-adic valuation of the
  multi-term forms `F_c(z)`, the input that `outside.md` §6 item 1 places outside the barrier. SL(n) itself is a
  statement about profiles; it concerns fakes only through their profiles, and a fake's profile satisfying SL says
  nothing about the fake (which has no `F_c`).

Neither hypothesis implies the other as far as I can see. Both say that the base locus of the slab-box linear
system has no non-toric part in a certain region: ZE near `e`, SL "at infinity in the direction `a_p`" (Prop 3.5,
Remark 3.6).

### 3.4 Equivalent forms

**Proposition 3.4 (`p`-adic separators).** Let `S` be finite, `V` a `Q`-subspace of `Q^S`, `l : S -> Z`, `k* in S`. Let
`Lambda = {(c_k p^{l(k)})_k : c in V (x) Q_p} cap Z_p^S` and `R = Lambda mod p`, an `F_p`-subspace of `F_p^S`.
* (i) Some `c in V` has `UM_p` with minimiser `k*` iff `e_{k*} in R`.
* (ii) For `V = V(S, m)` (the full space): `e_{k*} notin R` iff there is a polynomial `g in Q_p[k']` of degree `< m` with
  `v_p(g(k)) >= l(k) - l(k*)` for all `k in S` and `v_p(g(k*)) = 0` (a "`p`-adic separator").
* (iii) (Corner form.) Suppose `S` is contained in `k* + {sum_i t_i g_i : t in N^{n-1}}` for a lattice basis `(g_i)` with
  `l(g_i) >= 1` for all `i`. Then `c in V(S,m) cap Z^S` with `p` not dividing `c_{k*}` has `UM_p`. In generating-function
  form: such `c` exist iff some `Phi in Z[x_1, ..., x_{n-1}]` with exponents in `{t : k* + sum t_i g_i in S}` vanishes to
  order `>= m` at `x = (1, ..., 1)` and has `p` not dividing `Phi(0)`. Over `Q` (ignoring `p`): iff the torus-fixed
  point `x = 0` of the affine chart at the vertex `k*` is not a base point of the linear system "exponents in `S`,
  order `>= m` at `(1, ..., 1)`".

*Proof.* (i) If `c` has `UM_p` at `k*`, rescale by a power of `p` so that `v_p(c_{k*}) + l(k*) = 0`; then
`c' = (c_k p^{l(k)})` is in `Lambda`, `c'_{k*}` is a unit and `c'_k in pZ_p` for `k != k*`, so `c' mod p` is a nonzero
multiple of `e_{k*}`. Conversely, if `c' in Lambda` reduces to `e_{k*}`, write `c' = (c_k p^{l(k)})` with
`c in V (x) Q_p`; `V` has a basis of rational vectors, so a rational `c'' in V` close enough to `c` `p`-adically has
`(c''_k p^{l(k)}) == e_{k*} mod p`, i.e. `UM_p` at `k*`.
(ii) `Lambda` is a saturated `Z_p`-submodule of `Z_p^S`; its orthogonal `Lambda^perp = {y in Z_p^S : <y, Lambda> = 0}`
is saturated, and their reductions are mutually orthogonal of complementary dimensions, so `R^perp = Lambda^perp mod p`.
For `V = V(S,m)`, `V^perp = Poly_{<m}|_S`, hence
`Lambda^perp = {(g(k) p^{-l(k)})_k : g in Poly_{<m}(Q_p)} cap Z_p^S`. So `e_{k*} notin R` iff some element of
`Lambda^perp` has a unit at `k*`, which is the stated `g` after multiplying by `p^{l(k*)}`.
(iii) Every `k != k*` in `S` has `l(k) >= l(k*) + 1` and `v_p(c_k) >= 0`, so `k*` is the unique minimiser when
`v_p(c_{k*}) = 0`. In the coordinates `t`, `c` corresponds to `Phi(x) = sum_t c_{k*+sum t_i g_i} x^t`;
`sum_k c_k P(k) = 0` for all `P` of degree `< m` iff `Phi` vanishes to order `m` at `(1,...,1)` (apply the operators
`(x_i d/dx_i)` and evaluate at `1`), and `c_{k*} = Phi(0)`. QED

### 3.5 Obstructions to SL, and why they do not apply generically

**Proposition 3.5 (necessary conditions for SL(P, p)).** Let `S` be contained in `A_p`.
* (i) (Levels.) If `h_S(a_p) < m`, no `c in V(S, m)` has `UM_p`. [Round-10 Prop 5.5(ii).]
* (ii) (Toric strata at infinity.) Let `a'_p = (a_{p,x} - a_{p,n})_{x<n}`, so that `l_p(k) = <k', a'_p>`. If `a'_p` lies
  in the `Q`-span of a cheap family `v_1, ..., v_s` for `(S, m)` (Prop 2.2), then no `c in V(S,m)` has `UM_p`.
* (iii) More generally, if the base locus of `V(S,m)` over a finite extension `K` of `Q_p` contains a point `z` of the
  torus whose valuation vector is exactly `a'_p` (normalisation `v(p) = 1`), then no `c in V(S,m)` has `UM_p`.

*Proof.* (iii) If `c` had `UM_p`, then in `sum_k c_k z^{k'}` the term with the unique minimal valuation
`v_p(c_k) + <k', a'_p>` dominates, so the sum is nonzero, contradicting `z in` base locus. (ii) Write
`N a'_p = sum_j mu_j v_j` with integers `N >= 1`, `mu_j`. In `K = Q_p(varpi)`, `varpi^N = p`, put `z_j = varpi^{mu_j}`. The
proof of Prop 2.2(i) works over any field of characteristic `0`: the function `k -> prod_j z_j^{<k', v_j>}` coincides on
`S` with a polynomial of degree `< m` with coefficients in `K` (tensor Lagrange in `y_j = <k',v_j>`). Hence
`sum_k c_k prod_j z_j^{<k',v_j>} = 0` for `c in V(S,m)`. But `prod_j z_j^{<k',v_j>} = varpi^{N <k', a'_p>} = p^{l_p(k)}`,
i.e. this is the point `z = (p^{a'_{p,x}})_x` written in the torus, whose valuation vector is `a'_p`. Apply (iii). QED

So SL(P, p) fails whenever `a'_p` is a rational combination of few cheap directions: duplicated columns, a column
that is constant outside a set of fewer than `m/(2T) ~ 0.37 n` points not containing the base point `n` (then `a'_p`
is supported on that set, and `h(e_x) <= 2T`), a column in the `Q`-span of fewer than `m/w` other columns with
`e_q = 1` (`h_S(a'_q) <= e_q w`), and mixtures. Low-rank profiles (round-10 Prop 5.2) are of this type, consistently with the fact that their positions lie
in `Z_P` automatically. *Heuristic (not a proof):* if, as for the coordinate and slab directions, nonzero lattice
directions have width `>~ 2T` on `A_p` unless they are slab directions (width `<= e_q w`), cheap subspaces have
dimension `<= 0.37 n + O(sqrt(n/ln r))`; a balanced column `a'_p` that is not a short rational combination of other
columns and coordinate vectors then avoids them, so (ii) is not expected to obstruct admissible profiles.
Pair separation also forces balance where it matters: on a squarefree profile, summing
`d_xy >= W/2 - 2 log^+ C` over all pairs gives `sum_q (log q/W) s_q (n - s_q) >= (n(n-1)/4)(1 - 4 log^+ C/W)`,
where `s_q` is the size of the minority side of column `q`; since `s_q(n - s_q) <= n^2/4`, almost all of the weight
sits on nearly balanced columns, and a nearly balanced column has `~n/2 > 0.37 n` nonzero entries in `b_p`, so it is
not in the span of the coordinate directions alone.

**Remark 3.6 (what SL looks like geometrically).** By Prop 3.4(iii), in the corner form SL asks that a torus-fixed
point of the toric variety of `conv(S)` is not a base point of the fat-point system "order `m` at `e`". By Prop 3.5,
any base locus of that system which reaches the boundary in the direction `a_p` (tropically) kills SL. ZE and SL are thus the same kind of statement, "no non-toric base locus", in two different regions:
near `e` (ZE) and on the ray `a'_p` at infinity (SL). A proof of a Seshadri-type positivity for `L_A - mE` away from
the closures of the cheap subtori would give both; I did not find one.

## 4. Angle (c): arithmetic of the actual point

**Proposition 4.1 (the actual point is on every profile ray, `p`-adically).** Put `w_x = z_x conj(z_n)/N in Q(i)`
(`x < n`), so that `w = e^{iu}`. For every split `p = pi conj(pi)`:
`v_pi(w_x) = a_{p,x} - a_{p,n} = a'_{p,x}` and `v_{conj pi}(w_x) = -a'_{p,x}`; for inert `q`, `v_q(w_x) = 0`.

*Proof.* `v_pi(z_x) = a_{p,x}`, `v_pi(conj z_n) = e_p - a_{p,n}`, `v_pi(N) = e_p`. Inert primes divide `z_x` and `z_n` to the
same power `f` and `N` to `2f`. QED

So `w` is a `Q(i)`-point of the base locus `Z_P` (Thm 4.1) whose valuation vectors at the places above split primes
are exactly `+-a'_p`, while its archimedean distance to `e` is `<= C e^{-W/4}`. Theorem 3.1 is precisely the
exploitation of the non-archimedean information at one place: by Prop 3.5(iii), SL(P, p) is the statement that no
point with valuation vector `a'_p` (any residues) lies in the base locus of `V_0(A_p)` *in the strong sense that one
element of the system has a monomial initial form*. A weaker sufficient condition, also delta-free except through
residues, is:

> **SL-res(P, p).** The `pi`-adic initial forms of the integer vectors of `V_0(A_p)` (face residue sums, round-10
> Prop 5.5) have no common zero at the residue vector `rho in (F_p^*)^{n-1}` of the cluster.

Theorem 3.1(ii) says that an actual cluster in the regime makes `rho` a common zero for every light `p`, and also at
the places `conj(pi)`. With `p` small, `rho` takes at most `p - 1` values per coordinate; collisions push towards common
zeros (round-10 Prop 5.5(iii)), which is why SL (monomial initial forms) is the natural target and SL-res is not
pursued here.

**Why the other arithmetic tools do not close the gap (summary, no new claims).** (a) Arithmetic Bezout bounds the
height of an *isolated* point of `Z_P` by `deg^{n-2}` times the height of the equations, exponential in `n` against
`h(w) ~ W/2 = O(n log n)` (round-10 §5.5). (b) Multiplicity estimates (Philippon, Masser-Wuestholz) need vanishing on
many points of a group; Theorem 4.1 forces the powers `w^j` into `Z_P(A)` only while the inherited cost `j w W`
stays below about `mW/2`, i.e. for `j` up to about `m/(2w) ~ sqrt(n/ln r)`, far too few. (c) Liouville at the archimedean place is exactly what Theorem 4.1 already spends. The
new ingredient in §3 is a non-archimedean Liouville at a single prime, which costs nothing in the LP beyond widening
one slab.

## 5. Evidence (small models; not used in any proof)

**Boolean SL tests** (`bool_sl.py`, `bool_sl2.py`, `bool_sl2b.py`, `bool_sl3.py`; outputs `*_out.txt`). A box slice
`k = k* + sum_i t_i v_i` (`t in {0,1}^M`, `v_i = +-e_x` chosen with `l_p(v_i) >= 1`) turns SL in the corner form
(Prop 3.4(iii)) into the following question, modelled abstractly with random slab vectors `u_q` and `l(t) = |t|`: on `X = {t in {0,1}^M : |<t, u_q> - c_q| <= s for all q}` (with `0 in X`), is there `phi` orthogonal
to all multilinear polynomials of degree `< m` with `phi(0) != 0` (char 0), with `phi(0) = 1` and `phi` integral at `p`
(`p`-integral)? Method: rank of the 0/1 evaluation matrix modulo the prime `2^31 - 1` for char 0 (a lower bound for the
rational rank, so a "yes" is certain only up to that caveat), and Smith reduction over `Z/3^19` for `p = 3`.
Instances: `M = 8..12`, `m = 3..5`, `r = 3..8` random `{0, +-1}` slab vectors, `s = 1` (three slab values, so every slab
polynomial has degree `3 < m` when `m >= 4`). A *product certificate* at `0` is a family of disjoint blocks whose
unions all lie in `X`; one exists iff the number of blocks can reach `m`.

* Every instance with `dim V(X, m) > 0` and a product certificate of order `m` satisfied SL, both in char 0 and
  `3`-integrally (Smith exponents all `0`).
* The instances with `dim V(X, m) > 0` and **no** product certificate (largest block family `< m`) were few, with
  `dim V(X,m) <= 3`, and in all of them SL failed in char 0 (`bool_sl2_out.txt`: one instance; `bool_sl3_out.txt`:
  two instances among 46 product-free random systems, the other 44 having `V(X,m) = 0`). One `bool_sl2` instance that
  looked product-free under a block search restricted to blocks of size `<= 2` has a product certificate with
  blocks of size `<= 5` (`bool_sl2b_out.txt`). So in every instance where the block search was run, SL held iff a
  product certificate existed.

**Widened slab box of the Sylvester-Hadamard profile, `n = 8`** (`had_sl.py`, `had_prod.py`; output
`had_sl_out.txt`). `T = 1`, narrow slabs `|<k, beta_q>| <= 2` for the six columns `q != p`, `|<k, beta_p>| <= 2m`,
`m = 4`; `|A_p| = 145`, `E = 120`, `dim V(A_p, 4) >= 54` (ranks mod `2^31 - 1`). For 10 sampled `k*`, SL in char 0 and
`3`-integrally held exactly at the `k*` on level `l* = -2` (seven levels above), and failed at `l* = -1, 0`. At the two
successful `k*` that were examined, an explicit product certificate (four unit generators with `l_p >= 1`, all 16
subset sums in `A_p`) exists; at the two failing ones none exists.

**Fourier-uniform sets** (`fourier_sl.py`, exact rational arithmetic after averaging over `AGL(k,2)`; output
`fourier_sl_out.txt`). Coordinates `F_2^k`, all `2^k - 1` nonzero characters as slabs,
`X = {Y : |sum_{y in Y} (-1)^{q.y}| <= s for q != 0}`. For `(k, s) = (3, 1)`: SL at the empty set holds for `m = 3` and
fails for `m >= 4`; for `(k, s) = (4, 2)` (`|X| = 1170`): holds for `m = 4, 5`, fails for `m >= 6`. In both cases SL fails
exactly once `m` exceeds the number `2s + 1` of slab values, i.e. once the model is locally supercritical; but in
these instances the invariant part of `V(X, m)` is then `0` for counting reasons (`|X|` is below `dim Poly_{<m}`), so
they do not test the regime either. A Fourier-uniform model is expected to become supercritical only around
`M = 2^k >= 64`, beyond enumeration.

Reading: at these sizes SL is decided by product (parallelepiped) structure, which is exactly the capped mechanism.
The regime where the slab box is supercritical (`m >= 4w`, i.e. `n >= 2137`) is far out of reach, so these tests
neither support nor refute SL there. They do show that SL is not automatic: it can fail for admissible-looking
small slab systems once products are excluded.

**Scripts** (all in `theory_checks/`, each writes the `*_out.txt` of the same name when run as in the outputs):
`sl_test.py` (first lattice test in dimension `2-3` with exact `p`-saturation; the supports there are images of boxes,
so certificates are products), `bool_sl.py` (Boolean SL; char 0 and `3`-integral), `bool_sl2.py` / `bool_sl2b.py` /
`bool_sl3.py` (with product-certificate search), `had_sl.py` / `had_prod.py` (Hadamard `n = 8` widened box),
`fourier_sl.py` (Fourier-uniform sets, exact after symmetry reduction). All are evidence only.

## 6. What remains open

1. **ZE itself** (near `e`): not proved. §2 shows the tangent cone at `e` cannot decide it alone, because actual
   positions may lie exponentially close (relatively) to low-rank toric strata (Prop 2.5); a multi-scale analysis of
   the base locus along the strata would be needed.
2. **SL(n)** (Theorem 3.3): not proved. It is the cleanest form of the missing step found so far: delta-free,
   residue-free, one vector in an explicit rational subspace with a monomial `p`-adic initial form. Concrete
   sub-questions, in increasing difficulty:
   * (a) the full-space version SL-full(P,p) (`V(A_p, m)` instead of `V_0(A_p)`): decide it for structured admissible
     profiles (Sylvester-Hadamard, Paley) at moderate `n`, where the symmetry group of the profile acts on `A_p` and
     reduces the linear algebra to orbit sums;
   * (b) bridge from SL-full to SL: by Lemma 1.1 with `H = e^{mW/40}` (allowed in Theorem 3.1),
     `codim(V_H in V(A_p,m)) <= ~14 E/n`; a statement that the unique-minimiser vectors are not confined to a
     subspace of that codimension would suffice;
   * (c) the corner form (Prop 3.4(iii)): show that a torus-fixed point of `X_{conv(A_p)}` whose cone is spanned by
     directions with `l_p >= 1` is not a base point of "order `m` at `e`", for slab boxes in the regime `m >= 4w`. This is
     a Seshadri-type positivity statement for `L - mE` away from the cheap subtori (Remark 3.6).
3. **Non-squarefree norms.** Theorem 3.1 needs no squarefreeness; Theorem 3.3 uses it only through the rank
   bound of round-10 Lemma 3.1. A version of SL(n) for profiles of rank `>= c n` would extend the conditional bound.
4. **Small-model evidence** is dominated by products (§5); the supercritical regime (`n >= 2137`) is out of reach
   of direct computation, so any further test of SL must exploit symmetry (2(a)).
