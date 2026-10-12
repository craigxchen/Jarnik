# Lang–Waldschmidt and the growth rate: a linear form sees only the support weight of its character, and the skew double core is LW-immune

Round 8, target: decide whether the Lang–Waldschmidt conjecture (single linear forms in logarithms over
`Q(i)`) implies `M = o(log R/loglog R)` for actual circles. Scripts are in `round8/lw_checks/`. Each runs
from that directory with `python3 <script>`; its output is stored next to it as `*_out.txt`.

Conventions as in `round4/outside.md` and `round5/construct.md`:
* *Theorem*, *Proposition* and *Lemma* mean proved here with all quantifiers, or cited with the exact source.
* *Exact check* means integer or `Fraction` arithmetic. *Evidence* means a finite computation that suggests
  an asymptotic statement but does not prove it.
* `W = log N = 2 log R`. An arc of length `<= C sqrt R` has angular width `Delta <= C e^(-W/4)`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**No growth improvement follows from Lang–Waldschmidt through single linear forms. For each fixed `eps`, `LW(eps)`
(in the task's form, for all forms in the logarithms of the primes and of `i`) gives no growth improvement at all. Whether the full conjecture (with `eps -> 0`) helps is reduced to a counting
problem and is open; the heuristic says it does not.**

1. **LW sees only the support weight (Lemma 2.2, Corollary 2.3, new).** For every representation of a character
   form `sum_j v_j Log u_j + b Log i` by numbers `alpha_k in Q(i)^*`, `sum_k log H(alpha_k) >= w(supp v)`, the total
   weight of the primes on which the character is nonzero; grouping the primes by value attains this up to
   `O(K log beta_K)`. So an LW inequality can exclude a cluster through one form only if some character has
   support weight `< W/4 + log(1/C_0)` (Theorem 3.2), and it does so whenever a character of moderate complexity has
   support weight `< W/(4(1+eps)) - O(K log beta_K)` (Theorem 3.1). Liouville, by contrast, sees the full weight
   `||v||_w/2`.
2. **Theorem L of `flipped.md` is correct** (re-derived, identities re-run). It is the label-pair case of
   Theorem 3.1. It kills every refereed fake family of rounds 5-6 (`construct.md`, `lower.md`, `sharp.md`: support
   fraction `~ 1/(b+1)`), as the task suspected: under LW none of them is a genuine circle for large `M`.
3. **The skew double core is immune (Theorem 5.1, new).** For the Paley skew double core of `exact.md` §4 (the
   Plotkin-optimal profile, best known fake constant), every nonzero character is nonzero on at least `M/2` of the
   `2(M-1)` types: a two-basis uncertainty principle for `K, K^T`, proved modulo `q`. Exactly, the minimum is
   `M/2 + 1` (`M = 8, 12`). Hence every character has support weight `>= W/4 + (1/2) log P - 1/16`, and **no single
   LW inequality, for any `eps >= 0`, with naive heights and exponent exactly `1`, excludes it** (Corollary 5.2).
   The lead's multi-basis idea is therefore not needed for single forms; its stated form ("support `>= (1-1/t)W`")
   is false because pairs always have support about `W/2` (§7).
4. **Simultaneous fakes for fixed `eps` (Theorem 6.1, new).** Angle models on the double core satisfying every
   Liouville gap and `LW_J(eps, C)` (the task's max form) for every linear form in the logarithms of their primes and
   of `i` (all groupings, all coefficients, rational characters and non-character monomials included; homogeneous
   forms without `Log i` and Waldschmidt's product form only at the single-form level) exist at
   `W <= (max(8, 4/eps) + o(1)) M log M`. So no argument from these inputs at a fixed `eps` proves
   `M <= (min(1/4, eps/2) - eta) log R/loglog R`.
5. **Open.** The `1/eps` comes only from the union bound (the double core is exactly at the threshold, margin
   `O(log P)`). An `eps`-uniform barrier follows from either a `mu`-spread family (Proposition 7.2, e.g. if `t >= 3`
   random Paley bases are spread; Lemma 7.1 handles all small-support characters) or a character-level count for
   the double core (§8).

| question | answer | where |
|---|---|---|
| does a single LW form exclude all `O(M log M)` profiles? | **no**: the skew double core is immune for every `eps >= 0` | Thm 5.1, Cor 5.2 |
| does `LW(eps)`, fixed `eps`, give growth? | **no**: fakes at `max(8, 4/eps) M log M` satisfy it for all forms | Thm 6.1 |
| does full LW (`eps -> 0`) give growth? | open; reduced to (S_mu) or (N); heuristically no | §7, §8 |
| is Theorem L right? | yes; it excludes all refereed round-5/6 fakes | §4 |
| multi-basis claim `(1-1/t)W`? | false (pairs: `W/2`); what matters is `sigma > 1/4` | §7 |


---------------------------------------------------------------------------------------------

## 1. Setting and the hypotheses used

### 1.1 Angle models and character forms (recalled)

As in `outside.md` §1 and `construct.md` §1: distinct primes `p_j = 1 mod 4` (`j in J`, `|J| = r`), weights
`w_j = log p_j`, exponents `e_j`, rows `a_x in prod_j [0, e_j]` (`x = 1..M`), `W = sum_j e_j w_j`. A real vector
`phi in R^J` realises the profile if

```text
<2 a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x        (x = 1..M),   delta_x in [0, Delta],  k_x in Z.        (1.1)
```

For actual (normalised, primitive) clusters, `phi_j = arg pi_j` for a fixed Gaussian prime `pi_j | p_j`, and
`u_j := pi_j / conj(pi_j)` has `Log u_j = 2 i phi_j + 2 pi i m_j`. For `v in Z^J` write
`w(supp v) = sum_{j : v_j != 0} w_j` (the **support weight**) and `||v||_w = sum_j |v_j| w_j` (the **Liouville
weight**; it is `w(v)` of `outside.md`).

A **character** is a zero-sum `lambda in Q^M` with integral image `v(lambda) = sum_x lambda_x a_x in Z^J`; `L` is
the lattice of character images. By (1.1), `2 <v(lambda), phi> = lambda.delta + (pi/2) lambda.k`. Hence for actual
clusters and every `b in Z` the **character form**

```text
Lambda(lambda, b) := sum_j v_j Log u_j + b Log i  =  i ( lambda.delta + (pi/2)(lambda.k + b) + 2 pi n ),   n in Z,    (1.2)
```

and for integral `lambda` and the right `b`, `Lambda = i lambda.delta`, so `|Lambda| <= (||lambda||_1/2) Delta`
(`lambda` is zero-sum and `delta in [0, Delta]^M`). This upper bound is the **only** way an actual cluster forces a
linear form in the `Log u_j` to be small: by `outside.md` Theorem A.1(i), `<v, phi>` is undetermined by the
cluster for `v` outside `L tensor R`. Clusters below are normalised as in `two_thirds/proof.md` Lemma 1.1
(Gaussian gcd removed), so the content `D` of `flipped.md` is `0`; a nonzero content only shortens the arc.

### 1.2 The Lang–Waldschmidt conjecture over `Q(i)`

For `alpha in Q(i)^*`, `H(alpha)` is the naive height (largest absolute coefficient of the primitive minimal
polynomial of `alpha` over `Z`), `h(alpha)` the absolute logarithmic height, and `A(alpha) := max(H(alpha),
e^(|Log alpha|))`, `Log` the principal logarithm.

**Hypothesis `LW(eps, C)`** (the form fixed by the task). For all `n >= 1`, all `alpha_1, ..., alpha_n in Q(i)^*`
and all `b in Z^n` with `Lambda := b_1 Log alpha_1 + ... + b_n Log alpha_n != 0`,

```text
|Lambda|  >=  C^n  B^(-(n-1+eps))  prod_k A(alpha_k)^(-(1+eps)),        B = max_k |b_k|.                      (LW)
```

The **Lang–Waldschmidt conjecture** asserts: for every `eps > 0` there is `C(eps) > 0` with `LW(eps, C(eps))`.
Sources: S. Lang, *Elliptic Curves: Diophantine Analysis*, Springer 1978, introductions to Chapters X and XI;
M. Waldschmidt, *Diophantine Approximation on Linear Algebraic Groups*, Springer 2000, Chapter 1 (§1.2, the
"optimistic" archimedean conjecture, there with `B (|b_1| ... |b_n|)^(-1-eps)` in place of `B^(-(n-1+eps))` and
`C(eps)^(nD)` for degree-`D` fields); M. Waldschmidt, *Open Diophantine problems*, Moscow Math. J. 4 (2004),
Conjecture 2.5 (the case of rational integers `alpha_k`). I did not re-verify the section numbering of the 2000
book (no network access to it here).

*Max form versus product form.* Since `B (|b_1| ... |b_n|)^(-1-eps) >= B^(-(n-1) - n eps)`, the product form at `eps`
implies (LW) at `n eps`; the converse fails for unbalanced coefficients (e.g. `b = (1, 1, B)`), so the product form
is the stronger statement. Below: the exclusion results (§3, §4) use (LW), hence hold under either form for a bounded
number of logarithms; the single-form immunity results (Theorem 3.2, Corollary 5.2) hold for **both** forms, because
`B prod_k |b_k|^(-1-eps) <= B^(-eps)` as well; the simultaneous Theorem 6.1 is proved for (LW) only.

*Normalisation.* Any true instance of `LW(eps, C)` has `C <= 1`. Indeed, appending `k` dummy terms
`1 * Log 1 = 0` to `Lambda = Log 2` keeps `Lambda`, raises `n` by `k`, keeps `B = 1`, and multiplies the right side
by `C^k` (as `A(1) = 1`); for `C > 1` and `k` large this exceeds `log 2`. So below `0 < C(eps) <= 1`.

**Weaker forms.**
* *Naive-height Baker-type forms* (`flipped.md` §4, hypothesis `(H_kappa)`): for at most `n_0` logarithms,
  `|Lambda| >= c_* B^(-kappa') prod_k H(alpha_k)^(-kappa)`, with `kappa >= 1`. `LW(eps, C)` gives it with
  `kappa = 1 + eps`; larger `kappa` is weaker.
* *Proved Baker-type bounds* (Baker–Wüstholz 1993, Matveev 2000): `log(1/|Lambda|) <= c(n) prod_k h(alpha_k) log B`.
  `flipped.md` Prop B shows these are useless here (product of heights); nothing below changes that.

### 1.3 What it means for a fake circle to satisfy LW

An angle model fixes genuine primes but fake angles `phi_j`. Its **fake logarithms** are defined on the group
`Gamma_J = < i, u_j (j in J) >` (the unit-modulus numbers built from the primes of `N` and `i`): for
`U = i^m prod_j u_j^(e_j)` put `Log_f U := i (2 <e, phi> + m pi/2)` reduced to imaginary part in `(-pi, pi]`.
For actual circles `Log_f = Log`.

**Definition 1.1 (`LW_J(eps, C)`).** An angle model satisfies `LW_J(eps, C)` if (LW) holds with every `Log`
replaced by `Log_f`, for all `n >= 2`, all lists `alpha_1, ..., alpha_(n-1) in Gamma_J`, `alpha_n = i`, and all
`b in Z^n` (`b_n = 0` allowed) with `Lambda_f != 0`. These are the linear forms "in the logarithms of the primes of
`N` and of `i`" of the task, with an arbitrary grouping of the primes into the numbers `alpha_k` (grouping matters,
see §2.3). Two remarks on scope:
* Listing `i` costs one factor `C B^(-1) A(i)^(-(1+eps))` against a *homogeneous* form (no `Log i` at all). Every
  form the cluster controls is a character form (1.2), which carries the winding term `b Log i`; the forms of
  `flipped.md` (Prop 2.1-2.2) are of this kind. Homogeneous forms are covered by the single-form statements
  (Theorem 3.2, Corollary 5.2), which hold for every representation, but not by the simultaneous Theorem 6.1.
* Forms involving numbers outside `Gamma_J` mix in genuine numbers; Lemma 2.2 shows that for a cluster-controlled
  form no representation by numbers of `Q(i)^*` gives a better lower bound than the support weight, which a
  `Gamma_J`-representation attains up to `O(K log beta_K)` (Corollary 2.3).


---------------------------------------------------------------------------------------------

## 2. What a linear form in logarithms can see: the support weight

### 2.1 Heights in `Q(i)`

For `alpha in Q(i)^*` and a finite place `v` of `Q(i)` (a Gaussian prime up to units, of norm `Nv`) put

```text
P+(alpha) = sum_{ord_v alpha > 0} ord_v(alpha) log Nv,        P-(alpha) = sum_{ord_v alpha < 0} |ord_v(alpha)| log Nv.
```

**Lemma 2.1.** Let `alpha in Q(i)^*`.
1. `2 h(alpha) = max(P+(alpha), P-(alpha))`.
2. If `alpha in Q`, then `log H(alpha) = h(alpha) = max(P+, P-)/2`.
3. If `alpha notin Q`, then `e^(2h(alpha)) <= H(alpha) <= 2 e^(2h(alpha))`.
4. If `U = i^m gamma/conj(gamma)` with `gamma = prod_j pi_j^(e_j^+) conj(pi_j)^(e_j^-)` (`e in Z^J`, `e != 0`), then
   `||e||_w <= log H(U) <= log A(U) <= ||e||_w + log 2 + pi`.

*Proof.* (1) With the usual normalisation, `2h(alpha) = sum_v log^+ ||alpha||_v + 2 log^+ |alpha|`, where
`||alpha||_v = Nv^(-ord_v alpha)` at finite places. The finite part is `P-(alpha)`. The norm of the fractional
ideal gives `2 log|alpha| = P+ - P-`, so `2 log^+|alpha| = max(0, P+ - P-)`. Add.
(2) For `alpha = a/b` in lowest terms, `H = max(|a|,|b|)`, and `P+ = log N(a) = 2 log|a|`, `P- = 2 log|b|`.
(3) Write `alpha = (A + B i)/D` with `gcd(A, B, D) = 1`, `B != 0`, `D >= 1`. The minimal polynomial is
`(D^2 X^2 - 2AD X + A^2 + B^2)/g`, `g` the content. Its Mahler measure is
`(D^2/g) max(1,|alpha|)^2 = max(D^2, A^2 + B^2)/g = e^(2h)`, and `H = max(D^2, 2|A|D, A^2+B^2)/g` lies between
this and twice it, since `2|A|D <= A^2 + D^2`.
(4) `ord_{pi_j} U = e_j = -ord_{conj pi_j} U`, so `P+ = P- = ||e||_w`; apply (3) and `|Log U| <= pi`. QED

*Check.* `check_valuation.py` Part A verifies (1)-(3) exactly on 2,806 random `alpha` (PASS).

### 2.2 The valuation lemma

**Lemma 2.2 (valuation lemma).** Let `v in Z^J`, `v != 0`, and suppose

```text
b_1 Log alpha_1 + ... + b_n Log alpha_n  =  sum_j v_j Log u_j + b' Log i + 2 pi i n'            (alpha_k in Q(i)^*, b_k, b', n' in Z).
```

Then `sum_k log H(alpha_k) >= w(supp v)`.

*Proof.* Terms with `b_k = 0` can be dropped (`log H >= 0`). Exponentiating, `prod_k alpha_k^(b_k) = i^(b') prod_j
u_j^(v_j)`. Fix `j in supp v` and put `x_k = ord_{pi_j} alpha_k`, `y_k = ord_{conj pi_j} alpha_k`. Then
`sum_k b_k x_k = v_j` and `sum_k b_k y_k = -v_j`. Let `Q2` be the indices with `alpha_k notin Q` and `Q1` the
rational ones; for `k in Q1`, `x_k = y_k =: z_k`. By Lemma 2.1,

```text
log H(alpha_k) >= (P+ + P-)/2 >= (1/2) sum_j (|x_kj| + |y_kj|) w_j     (k in Q2),
log H(alpha_k) =  max(P+, P-)/2 >= (1/2) sum_j |z_kj| w_j               (k in Q1, as P+ + P- >= 2 sum_j |z_kj| w_j).
```

So it suffices that `tau_j := (1/2)[ sum_{Q2} (|x_k| + |y_k|) + sum_{Q1} |z_k| ] >= 1` for each `j in supp v`.
Subtracting the two valuation equations, `sum_{Q2} b_k (x_k - y_k) = 2 v_j != 0`, so some `k_1 in Q2` has
`x_{k_1} != y_{k_1}`. If `|x_{k_1}| + |y_{k_1}| >= 2`, then `tau_j >= 1`. Otherwise `|x_{k_1}| + |y_{k_1}| = 1`, so
`x_{k_1} + y_{k_1} = +-1`. Adding the two equations, `sum_{Q2} b_k (x_k + y_k) + 2 sum_{Q1} b_k z_k = 0`, and the
`k_1` term is `+-b_{k_1} != 0`. So some `k_2 != k_1` has `x_{k_2} + y_{k_2} != 0` (`k_2 in Q2`) or `z_{k_2} != 0`
(`k_2 in Q1`); either way it contributes at least `1/2` to `tau_j`, and `k_1` contributes `1/2`. QED

The bound is attained up to `O(1)` by grouping the primes by the value of `v` (Corollary 2.3(b)). A single number
can carry two primes at half price (`alpha = pi_1/pi_2` has `2h = max(w_1, w_2)`), but it then needs a rational or
conjugate partner, and the total is again `>= w_1 + w_2`. *Check:* `check_valuation.py` Part B finds, by brute
force over all representations with at most two numbers from a 624-element pool (exponents in `[-2,2]` at
`p = 5, 13`) and coefficients `|b| <= 4`, that the least `sum log H` is `>= w(supp v)` for six targets; the
minimum is within `0.18` of `w(supp v)`.

### 2.3 The LW cost of a character

**Corollary 2.3.** Let `Lambda = sum_j v_j Log u_j + b Log i` with `v != 0`, and let `0 < C <= 1`.
* (a) *(Lower bound.)* For every representation of `Lambda` as `sum_k b_k Log alpha_k` with `alpha_k in Q(i)^*`, the
  right side of (LW) is at most `B^(-eps) e^(-(1+eps) w(supp v))`.
* (b) *(Upper bound.)* Let `beta_1 < ... < beta_K` be the distinct values of `|v_j|` on `supp v`. Grouping the
  primes by value, `U_k = prod_{v_j = beta_k} u_j prod_{v_j = -beta_k} u_j^(-1)`, writes `Lambda` with `K + 1`
  logarithms (the last `Log i`) and coefficient bound `B <= (2K+1) beta_K + |Lambda|`, and (LW) then reads

  ```text
  log(1/|Lambda|) <= (1+eps) w(supp v) + (K+eps) log((2K+2) beta_K) + (K+1)[(1+eps)(pi + log 2) + log(1/C)].
  ```

*Proof.* (a) `C^n <= 1`, `B^(-(n-1+eps)) <= B^(-eps)`, `A >= H`, and Lemma 2.2. (b) `Lambda = sum_k beta_k Log U_k
+ b'' Log i` with `b'' in Z`; since every `|Log| <= pi`, `|b''| pi/2 <= |Lambda| + K beta_K pi`. Use Lemma 2.1(4)
for `A(U_k)` and `A(i) = e^(pi/2)`. QED

**Definition 2.4.** The *LW cost* `cost_eps(lambda, b)` of a character form is the infimum, over all
representations of `Lambda(lambda, b)` by numbers of `Q(i)^*`, of `log(1/RHS of (LW))`; so (LW) gives exactly
`log(1/|Lambda|) <= cost_eps`. By Corollary 2.3 (for `|Lambda| <= 1`, `C = C(eps)`),

```text
(1+eps) w(supp v)  <=  cost_eps(lambda, b)  <=  (1+eps) w(supp v) + (K+eps) log((2K+2) beta_K) + (K+1)(9 + log(1/C(eps))),
```

so **an LW inequality sees a character only through its support weight**, up to `K log beta_K`. By contrast
Liouville (the axis gap, `outside.md` Lemma 0.1) sees it through `||v||_w / 2`. LW beats Liouville exactly for
**concentrated** characters: large coefficients on a light support, `(1+eps) w(supp v) < ||v||_w/2`.


---------------------------------------------------------------------------------------------

## 3. The single-form criterion

**Theorem 3.1 (exclusion).** Assume `LW(eps, C(eps))`. Let an actual normalised cluster of `M` points lie on an arc of
length `<= C_0 sqrt R` of `x^2 + y^2 = R^2`, and let `lambda` be an integral character with `n = ||lambda||_1`,
image `v = v(lambda) != 0`, and `K`, `beta_K` as in Corollary 2.3(b). Then

```text
W/4  <=  (1+eps) w(supp v) + (K+eps) log((2K+2) beta_K) + (K+1)(9 + log(1/C(eps))) + log(n C_0/2).       (3.1)
```

*Proof.* Choose `b` in (1.2) so that `Lambda(lambda, b) = i lambda.delta`. It is nonzero, because
`prod_j u_j^(v_j) i^b = 1` forces `v = 0` (unique factorisation in `Z[i]`; `pi_j`, `conj pi_j` are distinct primes).
So `0 < |Lambda| <= (n/2) Delta <= (n C_0/2) e^(-W/4)`, and Corollary 2.3(b) gives the lower bound. QED

For a rational character apply (3.1) to an integral multiple. (3.1) with `K, beta_K` of size `M^O(1)` reads
`(1+eps) w(supp v) >= W/4 - O(K log M)`. Since `W >= 3 M log M - O(M)` and `K <= r`, **an LW inequality
excludes a cluster through a single form essentially exactly when the cluster has a character of moderate
complexity whose support weight is below `W/(4(1+eps))`** (Theorem 3.2 is the converse), i.e. a hyperplane
`lambda^perp` that contains columns of total weight more than about `(3/4 + eps/4) W` (for `e = 1`).

**Theorem 3.2 (no single form can exclude a spread profile).** Let a profile (with any `delta`, `k` and realising
`phi`) have `w(supp v) >= W/4 + x_0` for every nonzero character image `v in L`, with `x_0 > log(1/C_0)`. Then for
every `eps >= 0`, every `C <= 1`, every character `lambda`, every `b` and every representation of
`Lambda(lambda, b)` by numbers of `Q(i)^*`, the right side of (LW) is smaller than the cluster's upper bound
`(||lambda||_1/2) Delta`. The same holds for `(H_kappa)` with `kappa >= 1`, `c_* <= 1`, any `kappa'`.

*Proof.* By Corollary 2.3(a) the right side is `<= e^(-(1+eps) w(supp v)) <= e^(-W/4 - x_0) < C_0 e^(-W/4) <=
(||lambda||_1/2) Delta` (`||lambda||_1 >= 2`). For `(H_kappa)`: `c_* B^(-kappa') prod H^(-kappa) <= prod H^(-1)
<= e^(-w(supp v))` by Lemma 2.2. Rational characters: the support of `v(lambda)` equals that of an integral
multiple. QED

So "pick one linear form, bound it below by LW and above by the arc" — the method of `flipped.md` Theorem L —
can never exclude a profile all of whose characters have support weight `> W/4 + log(1/C_0)`, **uniformly in
`eps >= 0`**, even at exponent exactly `1` on naive heights.

**Corollary 3.3 (LW forces quarter-spreading).** Under the Lang–Waldschmidt conjecture, for every `eps > 0` there
is `M_1(eps, C_0)` such that every actual cluster with `M >= M_1` points on an arc `C_0 sqrt R` and
`W <= M^2` satisfies: every integral character `lambda` with `K + 1 <= W/(log M)^3`, `beta_K <= e^((log M)^2)` and
`||lambda||_1 <= e^((log M)^2)` has

```text
w(supp v(lambda)) >= W/(4(1+eps)) - W/log M.
```

*Proof.* In (3.1), `log(2K+2) <= 3 log M`, so the error terms are at most
`(W/(log M)^3)((log M)^2 + 3 log M + 9 + log(1/C(eps))) + (log M)^2 + log C_0 <= W/(2 log M)` for large `M`,
because `W >= 3 M log M - O(M)` (`two_thirds/proof.md`). QED

(Clusters with `W > M^2` already satisfy `M = O(sqrt(log R))`, so the restriction costs nothing for the growth
question.)

**Where the known profiles sit.** Write `sigma := min_v w(supp v)/W` over nonzero character images. For every
profile, the pair characters `e_x - e_y` have `w(supp) = ` weight of the columns separating `x, y`, and by the
cut bound `sum_{x<y} w(sep(x,y)) <= W M^2/4`, so `sigma <= M/(2(M-1))`. The LW window is therefore
`1/4 < sigma <= 1/2 + O(1/M)`.

| profile | lightest character | `sigma` | single-form LW |
|---|---|---|---|
| pure Hadamard core, `b` copies (walsh-extraction Thm B) | column `H(.,a)` | `1/(M-1)` | excludes (already excluded unconditionally) |
| core + one flipped copy per row, `b` copies (`construct.md`, `b = 8`; `lower.md`, `b = 16`; `sharp.md` P3, `b >= 33`) | column `H(.,a)`: block `a` plus all flipped columns | `~ 1/(b+1)` | **excludes** (Theorem 3.1, `K <= 3`, `beta_K = M`) |
| fully flipped core with flipped mass `F` | label pair `(H(.,a) - H(.,b))/2` | `~ F/(2W)` | excludes iff `F < W/2` (Theorem L) |
| **skew double core** (`exact.md` §4) | column of `K` | `(M/2 + 1)/(2(M-1))`, `>= M/(4(M-1))` proved | **immune** (§5) |
| `t` random row permutations of a Paley core | pairs: `1/2`; columns: `~ 1 - 1/t` | `<= 1/2`; small `M`: see §7 | unknown asymptotically (§7) |


---------------------------------------------------------------------------------------------

## 4. `flipped.md` Theorem L, re-checked

**Verdict: correct as stated.** I re-derived every step:
* Lemma 1.1, (1.2)-(1.4): the identity `sum_x H(x,a) s_xj = sigma_j [M delta_(a,a(j)) - 2 c_j(a)]` and the error
  `<= M Delta/2`; `d_j(a,b) in Z`. The exact identities are re-run: `check_two_log_identity.py`, copied to
  `lw_checks/rerun_flipped/`, passes all 2,600 exact identities with output identical to the stored one.
* Prop 2.2: `2 x (1.3)` gives `M Log u_ab - 4 Log w_ab = i pi (m_a - m_b) + 2 i eta + 2 pi i n`; with
  `b_3 = -2(m_a - m_b) - 4n`, `|Lambda_ab| = 2|eta| <= M Delta`; `|b_3| <= 2M + 9` when `M Delta <= pi/2`.
  Heights: `U_ab = G_a conj(G_b)` and `Q_ab` are conjugate-primitive of odd norm (disjoint primes), so
  `H(u_ab) <= 2 e^(W_a + W_b)`, `H(w_ab) <= 2 e^(F_ab)` (this is Lemma 2.1(4) above). Nonvanishing via the axis
  gap when some column of `a` or `b` has `|F_j| < M/4`.
* (L.1) is `(H_kappa)` plus logarithms. (L.2): the average of `W_a + W_b` over the `binom(M-1,2)` pairs is
  `<= 2W/(M-1)`; a row lies in `D_ab` for `P_x N_x <= (M-1)^2/4` pairs, so the average of `F_ab` is
  `<= F (M-1)/(2(M-2))`. The final inequality `(kappa delta/4) W <= E_M` and `W >= (M-1) log 5` are right.
* `LW => (H_(1+eps))`: with three logarithms, `kappa' = 2 + eps` in the task's form (`2 + 3 eps` in the product
  form), `c_* = C(eps)^3 e^(-3 pi (1+eps))` (as `A(alpha) <= e^pi H(alpha)` for `|alpha| = 1`).
  `M_0` depends on the (ineffective) `C(eps)`.

**Theorem L is the label-pair case of Theorem 3.1.** `Lambda_ab` is the character form of
`lambda = (H(.,a) - H(.,b))/2` (up to orientation signs), whose image is supported on the blocks `a, b` and on the
modified columns meeting `D_ab`; its support weight is `<= W_a + W_b + F_ab`, with equality in the one-flip-per-
column case. The single-label form of Prop 2.1 is the column character `H(.,a)`, support `W_a + W_F` in the fully
flipped class. Averaging over pairs, `F_ab ~ F/2`, gives the threshold `F < W/2`; single labels give `F < W/4`.
The phrase of the task "Hadamard orthogonality produces integral characters with coefficient `~M` on a low-height
block" is exactly this: `v(H(.,a)) = M/2` on block `a`, so `K <= 3` and `beta_K = M`.

**Consequences for the refereed fakes (checked against their definitions).** The profiles of `construct.md`
(`b = 8` unflipped copies per label plus one flipped copy per row; the residue design of its Lemma 4.1 uses
auxiliary primes that are residue data, not prime factors of `N`), of `lower.md` (`b = 16` copies plus capacity-two
twins) and of `round6_sharp/sharp.md` (P3-char, `b >= 33`) have modification mass `F/W -> 1/(b+1)`, i.e. `1/9`,
`1/17`, `<= 1/34`, and every label keeps unmodified columns. So under `(H_kappa)` with `kappa < (b+1)/2`, and in
particular under the Lang–Waldschmidt conjecture, **none of them is realised by an actual circle once
`M >= M_0(C, eps)`**. The task's statement "the fakes of `construct.md` have flip mass `W/9`, so under LW their
genuine versions do not exist" is correct in this sense (for large `M`, with an ineffective `M_0`). The same holds
for every Hadamard-core profile in the "Prop-P barrier class" of `walsh-extraction.md` and for any modification
pattern with `F <= (1/2 - delta) W`.

**What Theorem L does not reach.** (i) The referee of `flipped.md` showed that for the canonical Paley assignment
the pair inequalities force `W_F <= W/2 + O(M log C)`, so under LW the canonical fully flipped families must have
`W_F/W -> 1/2`: the *critical band*. (ii) The skew double core of `exact.md` §4 is a Hadamard core in which the
second copy of label `a` is the first copy modified on `{0, a}` (§5.1), so `F = W` and Theorem L is silent. §5
shows that this is not a gap of the proof: no single form excludes it.


---------------------------------------------------------------------------------------------

## 5. The skew double core is immune to every single LW form

### 5.1 The profile

Let `q = 3 mod 4` be prime, `M = q + 1`, `chi` the Legendre symbol mod `q` (`chi(0) = 0`), and index rows and columns
by `{0} cup F_q` (`0` is the border). With the Jacobsthal matrix `Q(x,y) = chi(y - x)` (skew, as `chi(-1) = -1`) put

```text
S = [[0, -1^T], [1, Q]],      K = I + S        (skew-Hadamard: K + K^T = 2I, K K^T = M I, K(x,0) = 1).
```

The **skew double core** (`exact.md` Theorem 4.1) has the `r = 2q` cut types, for `a in F_q`,

```text
s_a^(1) = K(., a),        s_a^(2)(x) = K^T(x,a) K^T(x,0)   (= s_a^(1) with the entries in rows 0 and a negated),
```

one prime `p_j = 1 mod 4` per type, `e = 1`, rows `a_xj = (1 + s_j(x))/2`. (The bracket is checked exactly in
`check_dcore_families.py` Part D.) So it is a Hadamard core with one unmodified and one modified copy of every
label, modification sets `{0, a}`. For an integral zero-sum `c`, `v(c)_j = <c, s_j>/2`, so
`supp v(c) = {j : <c, s_j> != 0}`. Every rational character `lambda` has `(M/2) lambda` integral
(`exact.md` Theorem 4.1(ii)), with the same support. Put `s(c) := |supp v(c)|`.

### 5.2 A two-basis uncertainty principle

**Theorem 5.1.** For every nonzero zero-sum `c in Z^M`, `s(c) >= (q+1)/2 = M/2`. Consequently, if all `p_j` lie in
`[P, P e^(1/(4r))]`, every nonzero character image `v` satisfies

```text
w(supp v)  >=  (M/2) log P  >=  W/4 + (1/2) log P - 1/16.
```

*Proof.* We may assume `c` primitive. For `a in F_q` put `Sigma_a := sum_{y in F_q} chi(a - y) c_y`. Since
`K(0,a) = -1` and `K(y,a) = [y = a] + chi(a - y)`, and `K(a,0) = K(0,0) = 1`, `K(0,y) = -1`,
`K(a,y) = [a = y] + chi(y - a)` for `y in F_q`,

```text
u_a := <c, s_a^(1)> = (c_a - c_0) + Sigma_a,          v_a := <c, s_a^(2)> = -(c_a - c_0) + Sigma_a.          (5.1)
```

So `u_a = v_a = 0` iff `c_a = c_0` and `Sigma_a = 0`. Let `Z_u, Z_v` be the zero sets in `F_q`; then
`s(c) = 2q - |Z_u| - |Z_v| = 2q - |Z_u cup Z_v| - |Z_u cap Z_v|`. Put `d = (q-1)/2`. Since `chi(z) = z^d mod q`
(also for `z = 0`), `Sigma_a = R(a) mod q` with

```text
R(X) = sum_{y in F_q} c_y (X - y)^d,     coefficient of X^(d-k) = (-1)^k binom(d,k) p_k,     p_k = sum_y c_y y^k,
```

and `binom(d,k) != 0 mod q` for `k <= d < q`.

*Case 1: `R != 0` in `F_q[X]`.* `Z_u cap Z_v` lies in the root set of `R`, so `|Z_u cap Z_v| <= d`. With
`|Z_u cup Z_v| <= q`, `s(c) >= 2q - q - d = (q+1)/2`.

*Case 2: `R = 0` in `F_q[X]`.* Then `p_k = 0 mod q` for `k = 0..d`, and by (5.1) `u_a = v_a = c_a - c_0 mod q`, so
`Z_u cup Z_v` lies in `E_0 = {a : c_a = c_0 mod q}`. Put `f_y = c_y - c_0`. As `sum_{y in F_q} y^k = 0 mod q` for
`0 <= k <= q - 2` (for `k = 0` the sum is `q`), `sum_y f_y y^k = p_k = 0 mod q` for `k = 0..d`. If `f != 0 mod q`,
the Vandermonde matrix `(y^k)_{k <= d, y in supp f}` has full column rank when `|supp f| <= d + 1`, so
`|supp (f mod q)| >= d + 2`, `|E_0| <= q - d - 2`, and `s(c) >= 2q - 2|E_0| >= q + 3`. If `f = 0 mod q`, then
`0 = sum_x c_x = c_0 + sum_y c_y = (q+1) c_0 = c_0 mod q`, so `c = 0 mod q`, contradicting primitivity.

*Weights.* `w(supp v) >= (M/2) min_j w_j` and `W <= 2q max_j w_j`, with `max - min <= 1/(4r) = 1/(8q)`. QED

This is an uncertainty principle for the pair of bases `K, K^T`: `u + v = 2 Sigma` is a polynomial of degree
`<= d = (q-1)/2` modulo `q`, and `u - v = 2(c - c_0)`. It is the two-basis analogue of `flipped.md` Lemma 6.2.

*Exact data* (`check_dcore_support.py`, `check_dcore_families.py`, `hyperplane_census.py`):
* the identities (5.1) on random `c`, `q = 7, 11, 19, 23`;
* **exact minimum** of `s(c)` over all `c != 0` (all hyperplanes spanned by types, `3,003` for `q = 7` and `646,646`
  for `q = 11`, maximisers re-verified in exact arithmetic): `min s = 5 = M/2 + 1` (`M = 8`) and `7 = M/2 + 1`
  (`M = 12`); exhaustive ternary `c` give the same minima, and small-support `c` at `q = 19` give `>= 18`;
* the columns `c = K(., a)` have `s = M/2 + 1` for `q = 7, 11, 19, 23, 43, 47`; label pairs `(K(.,a)-K(.,b))/2`
  have `M/2 + 3`; row pairs have `M`;
* census of light characters (primitive, up to sign, determined by their zero sets): `M = 8`: `s = 5, 6, 7, 8` occur `14, 28, 140, 525` times;
  `M = 12`: `s = 7, 8, 9, 10` occur `22, 0, 220, 1771` times. The lightest ones are exactly the `2q` columns.

So the bound of Theorem 5.1 is off by one, and the double core sits at `sigma = (M/2 + 1)/(2(M-1)) -> 1/4`: at the
edge of the LW window, but on the safe side by a margin of order `log P`.

### 5.3 Single-form immunity, uniformly in `eps`

**Corollary 5.2.** Consider any angle model on the skew double core with primes in `[P, P e^(1/(4r))]`,
`P > e^(1/8) C_0^(-2)`, and arc constant `C_0`. For every `eps >= 0`, every `C <= 1`, every character `lambda`, every
`b` and every representation of `Lambda(lambda, b)` by numbers of `Q(i)^*`, the (LW) lower bound is smaller than the
arc's upper bound `(||lambda||_1/2) Delta`. The same holds for `(H_kappa)` with any `kappa >= 1`, `c_* <= 1`.

*Proof.* Theorem 3.2 with `x_0 = (1/2) log P - 1/16 > log(1/C_0)`. QED

In particular the label forms (`H(.,a)`, support `(M/2 + 1) w`) and the label-pair forms of Theorem L (support
`(M/2 + 3) w`) of the double core have LW costs at least `(1+eps)(W/4 + (3/2) log P - 1/16)` and
`(1+eps)(W/4 + (7/2) log P - 1/16)`, above the arc's `W/4 - log(n C_0/2)`. **The single-form method cannot exclude a profile of size
`W = Theta(M log M)` that satisfies every Liouville gap, every residue axiom of the round-5 class
(`exact.md` Theorem 5.1, if correct) and has `M` points on an arc `C_0 sqrt R`.** Theorem 6.1 gives an independent
proof of the existence of such angle models (without residue data), together with simultaneous LW.


---------------------------------------------------------------------------------------------

## 6. Simultaneous LW: fake circles for each fixed `eps`

Single-form immunity says that no *one* comparison works. A fake circle that satisfies all LW inequalities
*simultaneously* (together with all Liouville gaps) rules out every argument from these inputs. For each fixed
`eps` this is achieved on the skew double core, at a size whose constant degrades like `1/eps`.

**Theorem 6.1.** Let `0 < eps <= 1` and `0 < C_0 <= 1`. There is `M_2(eps, C_0)` such that for every prime
`q = 3 mod 4` with `M = q + 1 >= M_2` and every `P` with

```text
P >= P_1 := max{ (40 M/C_0)^4 ,  (64 M)^(2/eps) },
```

there is an angle model on the skew double core of order `M` (one prime `= 1 mod 4` per cut, all in some
`[P', P' e^(1/(8q))]` with `P <= P' <= 2P`, `e = 1`, windings `k = e_1`) such that:
1. it has `M` distinct points on an arc of length `<= C_0 sqrt R`, and `W <= 2(M-1) log(2P) + 1/4`;
2. every monomial `v in Z^J \ {0}` satisfies the gaps of `outside.md` Lemma 0.1;
3. it satisfies `LW_J(eps, C)` (Definition 1.1) for every `C <= 1`, in particular for `C = C(eps)`.

With `P = P_1`: `W <= (max(8, 4/eps) + o(1)) M log M + 8 M log(1/C_0)`.

*Proof.* Write `r = 2q`, `J` = the cuts.

*Primes.* By the prime number theorem for the progression `1 mod 4`, `[P, 2P]` contains at least `12 q^2` such
primes for `P >= P_1` and `M` large; cutting `[log P, log 2P]` into `< 8q` intervals of length `1/(8q)`, one of
them holds `r` primes. So `w_min >= log P'` and `w_max - w_min <= 1/(8q)`.

*Margins.* For an integral character `c` with `n = ||c||_1`: `||v(c)||_w = sum_j w_j |<c,s_j>|/2 >= (w_min/2) Phi(c)`
with `Phi(c) = sum_j |<c, s_j>| >= 2(M-1) + n` (`exact.md` Theorem 4.1(iii); I re-derived its proof — the defect
identity, the parity of `u_a, v_a`, the cases `c_0 = 0`, `c_0 != 0` with `n = 4` and `n >= 6`, and amplitude
`>= 2` — and `check_dcore_families.py` Part B spot-checks it on random `c` for `q <= 47`). Hence

```text
||v(c)||_w - W/2 >= (n/2) log P' - 1/8,              w(supp v) - W/4 >= (1/2) log P' - 1/16   (Theorem 5.1).   (6.1)
```

For a rational non-integral character `lambda = lambda_0 1 + c` (`lambda_0 in (2/M)Z \ Z`), apply (6.1) to
`(M/2) lambda`: `||v(lambda)||_w >= W/M + (||lambda||_1/2) log P' - 1/(4M)`, and the support is unchanged.

*Random choice.* Draw `delta` uniformly from `[0, Delta]^M`, `Delta = C_0 e^(-W/4)`; given `delta`, draw the
realisation `phi` from the Haar measure of the realisation torus (`outside.md` Theorem A.1(i)). For a character
`lambda`, `2 <v(lambda), phi> = lambda.delta + (pi/2) lambda.k` depends on `delta` only; for `v` outside
`span_R L`, `<v, phi> mod 2 pi` is Haar-uniform for every `delta` (Theorem A.1(ii), proof).

*(a) Lemma 0.1 gaps for characters.* For integral `c`, the gap requires
`dist(c.delta/2, (pi/4) Z) >= arcsin(e^(-||v||_w/2))`, which fails with probability at most
`4 e^(-||v||_w/2)(n + pi/Delta)` (`exact.md` Lemma 5.0: density `<= 2/Delta`, at most `2n Delta/pi + 2` windows).
By (6.1) this is `<= 4.3 (n + pi/C_0) P'^(-n/4)`. There are at most `(2M)^n` vectors `c` with `||c||_1 = n`, and
`2 M P'^(-1/4) <= C_0/20`, so the total over `n >= 2` is `<= 4.3 sum_n (n + pi/C_0)(C_0/20)^n < 1/16`. Rational
non-integral `lambda` with `||lambda||_1 Delta <= pi/M` satisfy the gap deterministically: the winding gives
`dist(<v,phi>, (pi/4)Z) >= (pi/4)(2/M) - ||lambda||_1 Delta/4 >= pi/(4M) > e^(-W/(2M)) >= e^(-||v||_w/2)`. Those with
`||lambda||_1 > pi/(M Delta)` have `||v||_w >= (||lambda||_1/2) log P'` astronomically large, and the same window
count (density `<= M/(2 Delta)`) gives total probability `< 1/16` for large `M`.

*(b) Lemma 0.1 gaps for non-characters.* By Theorem A.1(ii), given (a), the bad realisations have Haar
proportion `<= 1.33(Pi - 1)` with `Pi = prod_j (1 + 2/(sqrt(p_j) - 1)) <= exp(2r/(sqrt(P') - 1)) <= 1 + 1/(200M)`
(as `sqrt(P') >= 1600 M^2`); so `< 1/8`.

*(c) LW_J.* Fix `K >= 1`, `e_1, ..., e_K in Z^J \ {0}`, `c in (Z\{0})^K`, and put `v = sum_k c_k e_k`,
`h' = sum_k ||e_k||_w`, `B_c = max_k |c_k|`. Torsion entries of a list only lower the right side of (LW) (they have
`A >= 1` and raise `n`), and `C^n <= 1`, `A(U_k) >= e^(||e_k||_w)` (Lemma 2.1(4)), `B >= B_c`. Since
`Lambda_f = i (2 <v, phi> + (pi/2) m)` for some `m in Z`, it suffices that

```text
dist(2 <v, phi>, (pi/2) Z)  >=  L(e, c) := B_c^(-(K+eps)) e^(-(1+eps) h')       for every such (e, c) with v != 0.    (6.2)
```

(If `v = 0`, `Lambda_f in i (pi/2) Z`, so `Lambda_f = 0` or `|Lambda_f| >= pi/2 > 1`.) The failure probability of
(6.2) is at most `(4/pi) L(e,c)` if `v` is not in `span_R L` (four windows on the circle), and at most
`4 M L(e,c)/Delta` if `v = v(lambda)` is a character image (density of `lambda.delta` `<= 1/(Delta ||lambda||_inf)`,
at most `2 ||lambda||_1 Delta/pi + 2` windows, `||lambda||_1 <= M ||lambda||_inf`, `||lambda||_inf >= 2/M`).

Sums. `#{c in (Z\{0})^K : B_c = B} <= K 2^K B^(K-1)`, so `sum_c B_c^(-(K+eps)) <= K 2^K (1 + 1/eps)`. For `s > 0` let
`Z(s) = sum_{e != 0} e^(-s ||e||_w) = prod_j (1 + 2/(p_j^s - 1)) - 1 <= 4r/(P'^s - 1)` (when `2r <= P'^s - 1`).
* Non-character part: `<= (4/pi)(1 + 1/eps) sum_K K (2 Z(1))^K <= (16/pi)(1 + 1/eps) Z(1) < 1/8`.
* Character part: every such `(e, c)` has `h' >= w(supp v) >= W/4 + x_0`, `x_0 = (1/2) log P' - 1/16 >= 0`
  (each `j in supp v` lies in some `supp e_k`; Theorem 5.1). Hence
  `(1+eps) h' >= (W/4 + x_0) + (eps/2)(W/4) + (eps/2) h'`, and the part is at most

  ```text
  (4M/Delta) e^(-W/4 - x_0 - eps W/8) (1 + 1/eps) sum_K K (2 Z(eps/2))^K  <=  (4M/C_0)(2/eps) e^(-eps W/8) * 4 Z(eps/2)/(1 - 2Z(eps/2))^2.
  ```

  With `P'^(eps/2) >= (64M)` we get `2 Z(eps/2) <= 1/4`, and `e^(-eps W/8) <= P'^(-eps q/4)`; the bound is `< 1/8`
  for `M >= M_2(eps, C_0)`.

The failure probabilities of (a)-(c) add up to less than `1`, so a good `(delta, phi)` exists. The points are
distinct (every pair of rows is separated by `M` cuts, `exact.md` Theorem 4.1(i)) and lie on the arc by
construction. `W = sum_j log p_j <= 2q (log P' + 1/(8q)) <= 2(M-1) log(2P) + 1/4`. QED

**Remarks.**
1. *Residue data.* Adding the point-level residue data of `construct.md` Lemma 4.1 (auxiliary primes of norm
   `<= L ~ 2M log M`, not dividing `N`) multiplies the Lemma 0.1 window by `L^(n/2)` exactly as in `exact.md`
   Theorem 5.1; it changes `(40M/C_0)^4` to about `(40 M L/C_0)^4` and does not touch (c). The constant becomes
   `max(12, 4/eps)`. I did not re-prove `exact.md` Theorem 5.1 (unrefereed); the residue part follows its proof
   verbatim.
2. *What Theorem 6.1 rules out.* For fixed `eps`, no argument whose inputs are the exponent profile, the
   Liouville gaps of all monomials, Haar-generic properties of the angles and `LW_J(eps, C(eps))` for all forms
   can prove `M <= (min(1/4, eps/2) - eta) log R/loglog R`. For `eps >= 1/2` the LW hypothesis costs nothing
   beyond the round-5 constant.
3. *What it does not rule out.* The constant `4/eps` comes from summing over all `K`-tuples of exponent vectors
   of height about `W/4`: there are about `exp((W/4) log r / log P)` of them, and only the factor `e^(-eps W/4)` of
   (LW) is available to beat this count, because the support-weight margin of the double core is only `O(log P)`.
   So an argument that uses `LW(eps)` with `eps -> 0` is **not** excluded by Theorem 6.1; see §8.


---------------------------------------------------------------------------------------------

## 7. Multi-basis profiles and the lead's suspicion

**Correction of the claim.** The lead suspected that profiles built from `t >= 2` independent Hadamard bases have
"every nonzero integral character of support weight `>= (1 - 1/t) W - small`". This is false for `t >= 3`, and for
every profile: a pair character `e_x - e_y` has support equal to the columns separating `x` and `y`, and the cut
bound gives `min_{x<y} w(sep(x,y)) <= W M/(2(M-1))`. For `t` row-permuted Hadamard bases with equal weights, every
pair has support exactly `t M/2` labels, `sigma = M/(2(M-1))`. The value `(1 - 1/t) W` is that of the *columns* of
one basis. What LW needs (Theorem 3.2) is `sigma > 1/4`, and the relevant window is `1/4 < sigma <= 1/2 + O(1/M)`.

**What is provable for `t` permuted Paley bases.** Let `H_i = Pi_i H` (`i = 1..t`) be row permutations of the
normalised Paley-I matrix of order `M = q + 1`, each contributing its `q` nonconstant columns with equal weights.

*Lemma 7.1.* Every nonzero zero-sum `c` with `|supp c| <= (1/4 - mu) M` has `|supp H_i^T c| >= (1/4 + mu) M - 1`
in every basis, hence `sigma(c) >= 1/4 + mu - O(1/M)`.

*Proof.* `|supp H_i^T c| = |supp H^T (Pi_i^T c)|`, and `flipped.md` Lemma 6.2 (stated for `H n`; the same proof
applies verbatim to `H^T`, since `(H^T c)_a` is, modulo `q`, a polynomial of degree `<= (q-1)/2` in `a` plus
`-c_a`) gives `>= (q+1)/2 - |supp c|` nonzero entries among the finite labels. QED

So only characters of large support can be light, e.g. columns of one basis (`sigma ~ 1/t + (1 - 1/t)(1 - O(M^(-1/2)))`
heuristically). Showing that no large-support character is sparse in two independent bases at once is a
"two-basis spark" problem for random permutations; I did not prove it.

*Data* (`hyperplane_census.py`, exact enumeration of all hyperplanes spanned by types, re-verified exactly):

| profile | `M` | types | min `s(c)` | `sigma = s/#types` |
|---|---|---|---|---|
| single Paley core | 8 | 7 | 1 | 0.143 |
| skew double core | 8 | 14 | 5 | 0.357 |
| skew double core | 12 | 22 | 7 | 0.318 |
| 2 random bases (seed 1) | 8 | 14 | 2 (two coinciding columns) | 0.143 |
| 2 random bases (seeds 1, 2, 3) | 12 | 22 | 4, 4, 7 | 0.18, 0.18, 0.32 |
| 3 random bases (seed 1) | 8 | 21 | 6 | 0.286 |

At these sizes random bases are *worse* than the skew double core: coincidences (two equal columns, or a 6-row
character sparse in both bases) are frequent. This is evidence only; it says nothing about large `M`, where such
coincidences should become rare, but it shows that "independent bases" is not automatically spread and that a proof
must control all characters, not just columns.

**Proposition 7.2 (spreading with a margin would give an `eps`-uniform barrier).** Suppose a family of profiles of
orders `M -> infinity` has `r = O(M)` types with one prime each in `[P', P' e^(1/(8r))]`, `e = 1`, satisfies the
Liouville margins `||v(c)||_w - W/2 >= (n/2) log P' - O(1)` for all integral characters (rational characters
being handled as in Theorem 6.1, e.g. by an injective winding), and is **`mu`-spread**: `w(supp v) >= (1/4 + mu) W`
for every nonzero `v in L`, with `mu > 0` fixed. Put `eta = 4mu/(1+4mu)`. Then the conclusions of Theorem 6.1 hold
with `P >= max{(40M/C_0)^4, M^(3/eta)}` for every `eps in (0,1]` and `M >= M_3(eps, mu, C_0)`; in particular
`W = O_mu(M log M)` with a constant independent of `eps`.

*Proof.* Only the character part of (c) changes. For a character tuple, `h' >= (1/4 + mu) W`, so
`W/4 <= h'/(1+4mu)` and `W/4 - (1+eps) h' <= -eta h'`. The character part is then at most
`(4M/C_0)(1 + 1/eps) sum_K K (2 Z(eta))^K <= (4M/C_0)(2/eps) * 4 Z(eta) * (16/9)` once `2 Z(eta) <= 1/4`, and
`Z(eta) <= 4r/(P'^eta - 1)`. With `P'^eta >= M^3` and `r = O(M)` this is `O(M^2/(eps C_0 M^3)) < 1/8` for
`M >= M_3`. Parts (a), (b) and the non-character part of (c) are as before. QED

So the lead's route would work for `eps`-uniform statements **if** some multi-basis family were proved `mu`-spread.
The skew double core is `0`-spread exactly (margin `O(log P)` only, Theorem 5.1), which is why Theorem 6.1 loses
`1/eps`.


---------------------------------------------------------------------------------------------

## 8. The decision, and what is left open

**(A) is refuted for single forms.** A conditional growth theorem "LW `=>` `M = o(log R/loglog R)`" proved by
finding, in every profile with `W = O(M log M)`, *one* linear form that LW forbids, does not exist: by
Corollary 2.3(a) and Theorem 3.2 a form is forbidden only through a character of support weight
`< W/4 + log(1/C_0)`, and the skew double core (Theorem 5.1) has none, for any `eps >= 0`, while it carries angle models with every Liouville gap at
`W = Theta(M log M)` (Theorem 6.1 for each fixed `eps`; `exact.md` Theorem 5.1 with residue data). This also covers
the weaker naive-height forms `(H_kappa)` (any `kappa >= 1`) and, a fortiori, Baker-type bounds.

**(B) holds for each fixed `eps`, with a constant that degrades.** Theorem 6.1: for fixed `eps`, fake circles
satisfying `LW_J(eps)` for every form, every Liouville gap and Haar-generic angle properties exist at
`W <= (max(8, 4/eps) + o(1)) M log M`. So no argument from `LW(eps)` at a fixed `eps` gives a growth improvement.

**What is not decided.** Whether the *full* conjecture (all `eps`, an argument free to take `eps -> 0` with `M`)
implies a growth improvement. Any such argument must use many forms at once (a covering of the `delta`-box by the
slabs `|lambda.delta| < LW-bound(lambda)`), since single forms never suffice (Corollary 5.2). The precise obstacle
in Theorem 6.1 is a counting problem: the union bound sums over all `K`-tuples of exponent vectors whose span meets
the character lattice, about `exp((W/4) log r/log P)` of them at the critical height, and only `e^(-eps W/4)` is
available against this count because the double core's support margin is `O(log P)`, not `Omega(W)`. Two ways to
close the gap, both open:
1. **(S_mu)** prove that some family satisfying the round-5 margins is `mu`-spread for a fixed `mu > 0`
   (Proposition 7.2), e.g. `t >= 3` random row permutations of a Paley core; Lemma 7.1 settles all characters with
   `|supp c| <= (1/4 - mu) M`, and the remaining case is a two-basis spark statement for random permutations;
2. **(N)** replace the union bound over representations by one over characters for the double core: bound the
   supremum of the (LW) right side over all representations of `Lambda(lambda, .)` by a function of `lambda` that
   decays fast enough in `||lambda||_inf`, and show that the number of characters with `||lambda||_inf <= T` and
   support excess `w(supp v) - W/4 <= x` is at most `e^(O(M) + x/2) T^(O(1))`. Theorem 6.1 instead sums over all
   representations, which overcounts by up to `e^((W/4) log r/log P)`.

*Heuristic (labelled as such).* The light characters of the double core are few and rigid: at `M = 8, 12` the
lightest are exactly the `2q` columns (excess `(3/2) log P`), the next level has `O(q^2)` members (label pairs,
excess `(7/2) log P`), and every character with `|supp c| <= d` has `s(c) >= M - 2|supp c|` (Lemma 7.1 applied
to `K` and `K^T`). A Khintchine–Groshev count weighted by the actual light characters (rather than by all
decompositions) then predicts an expected number of LW violations of order `M^O(1) P^(-3/2)` for random `delta`,
which tends to `0` with no help from `eps`. If this is right, (B) holds `eps`-uniformly and the
Lang–Waldschmidt conjecture gives **no** growth improvement for actual circles; the rate `log R/loglog R` would then
be beyond every single-form transcendence hypothesis, and only multi-form (simultaneous approximation /
Subspace-type) statements could improve it.

**Relation to the barrier map.** Before this round the refereed fake families (construct, lower, sharp) were all
`sigma ~ 1/(b+1)` single cores, and Theorem L showed that LW kills them. The remaining LW-proof of concept was
therefore "LW might kill every `O(M log M)` fake". Theorem 5.1 shows that the Plotkin-optimal double core of
`exact.md` sits exactly on the safe side of the LW threshold `sigma = 1/4` by an additive `O(log P)`: the critical
band of the `flipped.md` referee report is not just a gap in Theorem L but is occupied by the best known fake.


---------------------------------------------------------------------------------------------

## 9. Files (`round8/lw_checks/`)

| script | what it checks | result |
|---|---|---|
| `dcore.py` | library: Paley skew-Hadamard `K` (asserts skewness, `K K^T = M I`, `K(x,0) = 1`), double-core cuts, permuted bases, exact rank/nullspace | used by all |
| `check_dcore_support.py` | (5.1) on random `c` (`q = 7..23`); exhaustive minima of `s(c)` (ternary, and entries `<= 2`); exact minimum via all hyperplanes (`q = 7, 11`) | min `s = M/2 + 1`; bound `M/2` holds |
| `check_dcore_families.py` | `s` of columns, label pairs, row pairs (`q <= 47`); Liouville margin `Phi >= 2(M-1) + n` spot check; second family = first modified on `{0,a}` | `M/2+1`, `M/2+3`, `M`; PASS; PASS |
| `hyperplane_census.py q profile seed slack` | census of light characters (primitive, determined by their zero set) for `dcore`, `rand2`, `rand3`, `single`; outputs `hyperplane_census_*_out.txt` | table of §7, counts of §5.2 |
| `check_valuation.py` | Lemma 2.1 exactly on 2,806 random `alpha`; Lemma 2.2 by brute force over representations | PASS |
| `rerun_flipped/` | copy of `round4/flipped_checks/`; `check_two_log_identity.py` re-run | 2,600 identities PASS, output identical |

Every script runs from `lw_checks/` with `python3 <script>`; outputs are the `*_out.txt` files. Floating point
(SVD) is used only to *find* candidate hyperplanes; every reported minimiser is recomputed with `Fraction`
arithmetic. The repository `/home/user/Jarnik` was not modified.

