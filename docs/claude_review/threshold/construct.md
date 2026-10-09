# Constructing `ord(c) > 2t` is impossible, and forced divisibility reaches ratio exactly 1, never more

Task: construct a measure `c` on a shell `S_{s,t}` with `ord(c) > 2t` (problem (B) of
`../tau_problem.md`), or quantify route (C): several shells, `N`-dependent coefficients, and forced
Gaussian divisibility. Scripts and logs are in `scratchpad/tau/construct/` (§5). Other agents' files in
`scratchpad/tau/` (`upper.md`, `verify.md`, `upper_checks/`, ...) are read but not imported.

Conventions:

* *Theorem*, *Proposition*, *Lemma* and *Corollary* mean proved here, with all quantifiers.
* *Data* means an exact finite computation. Data is never used as a proof of a general statement.
* A modular rank gives a rigorous **upper** bound for `reg`, because the rank mod q is at most the
  rank over Q. Lower bounds come from explicit exact measures.

---------------------------------------------------------------------------------------------

## 0. Verdict

**(B) cannot be done.** No measure on any shell or slice, for any `M`, has `ord(c) > 2t`. There are
two independent proofs.

1. `upper.md` (concurrent agent) proves `ord(c) <= (2 - 1/ceil(M/2)) t`. I re-checked every step of
   its proof (§1.1) and reproduced its numbers with my own code (§1.2).
2. Theorem B below gives an independent self-contained proof of the weaker bound
   `ord(c) <= (2 - 2^{2-M}) t`, which is still below `2t` (§2.4, Corollary B2).

Every family the task suggested is therefore capped. These include alternating orbits, other isotypic
components, products and sums of Vandermondes, Schur and Jacobi–Trudi or generalised alternants,
determinants of other character sets, and hyperoctahedral compositions. Exact data for each family are
in §1.3.

**(C), forced divisibility, also cannot prove the uniform theorem, but it closes the gap exactly to the
boundary.** The main results are these.

* **Theorem C1 (forced divisor = tropical bound).** For any form `F` and any `M` points of norm `N`,
  the Gaussian divisibility of `F(z)` that holds for every choice of points with a given valuation
  profile is the tropical bound. At a split prime `p` it is governed by the width `h_A(b_p)` of the
  support `A` of `F` in the direction `b_p`. Here `b_p` is the centred `pi`-adic exponent vector of the
  points (definitions in §2.1). "Products of pair differences inherit `prod g_ij`" is the special case of
  root binomials.
* **Theorem C2 (exact criterion).** With this divisor, the degree `D` cancels. The Liouville step
  forces `F(z) = 0` exactly when

  ```text
  sum_p omega_p h_A(beta_p)  <  ord_Y F - O(1/log N)      (omega_p = e_p log p / log N)
  ```

  Every profile can be replaced by a **sign profile**, i.e. a measure on `{+-1}^M`, with the same pair
  separations and no smaller left side. So the ratio needed is

  ```text
  ord(c) / Lambda(A) > 1,
  ```

  where `Lambda(A)` is the value of an explicit linear programme over half-separating sign measures
  (§2.2). The pair inequality is the case `F = z_i - z_j`. The project's character inequalities are the
  case `A = {c, -c}`.
* **Theorem B (mean-width bound; new).** Let `S` be a finite subset of `Z^n`, and let
  `w_P(S) = max_S x(P) - min_S x(P)` for `P` a subset of `[n]`. Then

  ```text
  reg(S)  <=  2^{1-n} * sum_{P subset [n]} w_P(S).
  ```

  Equality holds for boxes. In homogeneous form, every auxiliary measure satisfies
  `ord(c) <= Omega(supp c) := E_eps h_{supp c}(eps)`, with `eps` uniform in `{+-1}^M`.
* **Corollary C3 (the cube profile defeats every auxiliary polynomial).** In the **cube profile**,
  each prime of `N` splits the cluster by an independent fair coin. Equivalently, the sign vectors
  `eps_p` run uniformly over `{+-1}^M`. This profile satisfies every pair inequality with equality,
  `d_ij = W/2`. At the exponent level, it also satisfies the forced-divisibility Liouville inequality
  of **every** `F`:
  every `M`, every degree, any number of shells, coefficients polynomial in `N`. So
  `Lambda(A) >= Omega(A) >= ord`, and the needed ratio `> 1` is never reached.
* **Proposition C4 (coefficients with arithmetic content).** Arbitrary `N`-dependent Gaussian-integer
  coefficients can carry any prime content. Then a uniformly random block-bijective assignment of the
  cube profile satisfies the inequality in expectation. So no `F` gains a margin proportional to
  `log N`.
* **Quantified answer to the task's question (§2.7).** The ratio needed is `ord / Omega > 1`.
  * Forced divisibility raises the Vandermonde from `ord/(2t) = (2s'+1)/(2s'+2)` to
    `ord/Omega = 1` exactly. The same holds for every product of pair differences. This is the Cilleruelo–Córdoba
    endpoint, now critical.
  * Theorem B says nothing exceeds `1`.
  * Against the best half-separating adversary, the Vandermonde ratio is `1 - 1/(2 ceil(M/2))`.

**Limits of these statements.**

* The barrier concerns the method: one auxiliary form, Liouville, and valuation-forced
  divisibility. It does not say that clusters with cube-like profiles exist on short arcs.
* Constant-size effects in the equality case of Theorem B (`ord = Omega`, for arcs `C < 1`) are
  settled only conditionally:
  * Conjecture E: the equality cases are exactly root-binomial products. Checked on 1662/1662 exact
    equality cases.
  * Conjecture C_bal: a balanced-profile version with slack. Checked on all tested sets.

  See §2.5 and §2.8.

---------------------------------------------------------------------------------------------

## 1. Problem (B): `ord(c) > 2t` is impossible

### 1.1 Re-verification of `upper.md`

`upper.md` claims `tau(M) = tau*(M) = 2 - 1/ceil(M/2)`. I checked every step it needs for the upper
bound.

* **Lemma 1.1 (order = regularity).**
  * If `c != 0` annihilates `Poly_{<m}|_A`, then `Poly_{<m}` does not restrict onto `Q^A`. Hence
    `reg(A) >= m`, and so `ord(c) <= reg(supp c)`.
  * Conversely, a kernel vector attains `reg(A)`.
  * Correct.
* **Lemma 2.1 (peeling).** Let `P = P_1 + L_1 P'`. Then `P = g` on `A_1`, and `P = P_1 + (g - P_1) = g`
  on `A \ Z(L_1)`. Its degree is at most `max(reg A_1, reg A' + 1)`. Induction with the index shift
  gives `reg(A) <= max_i (reg A_i + i - 1)`. Correct.
* **Lemma 2.2 (slices of `Q^M(p,n)`).** `(c, k')` lies in `Q^M(p,n)` iff `k'` lies in
  `Q^{M-1}(p-c, n)` (for `c >= 0`) or in `Q^{M-1}(p, n-|c|)` (for `c < 0`), with `-n <= c <= p`.
  Correct.
* **Lemma 2.3 (sorted peeling).** `N(r_(i)) >= i` and `N(v) <= K - v + 1` together give
  `r_(i) + i - 1 <= K`. Correct.
* **Theorem 3.1, case `j = 1`.**
  * The bounds are `r^+_c <= (p-c) + (m-1)n` and `r^-_c <= p + (m-1)n`.
  * `N(v) <= (V - v + 1) + n = K - v + 1`, with `K = p + mn`.
  * Correct. The second sub-case is not even needed.
* **Theorem 3.1, case `j = m`.** This is `k -> -k` together with `l^M_m(p,n) = l^M_1(n,p)`. Correct.
* **Theorem 3.1, case `2 <= j <= m-1`.**
  * `alpha = (m-1)/(m-j)` and `beta = (m-1)/(j-1)`, so `1/alpha + 1/beta = 1`.
  * The induction is used with `j' = j` on the positive slices and `j' = j-1` on the negative slices.
    Both are admissible indices for `M - 1`.
  * Identity (3.1): `V+/alpha + V-/beta = m p/(m+1-j) + m n/j = l^M_j(p,n)`. I re-derived it by hand.
  * Every slice value is at most `min(V+, V-)`. This uses monotonicity of `l^{M-1}` in both arguments,
    applied with both indices.
  * The count is `N(v) <= (x+ + 1) + x- = K - v + 1`.
  * Correct.
* **Corollary 3.2.**
  * `l^{2k}_k(p,n) = (2k-1)(p+n)/k`.
  * The mean of `l^{2k+1}_k` and `l^{2k+1}_{k+1}` is `(2k+1)(p+n)/(k+1)`.
  * Correct.
* **Lemma 1.2 (slices, not shells, are the right objects).** This was confirmed independently by
  `verify.md` and `referee_conditional.md`. I use it as given. Shells are subsets of slices and `reg` is
  monotone, so it does not affect the shell question asked here.

Conclusion: for every `M`, every shell `S_{s,t}` and every `c != 0` on it,
`ord(c) <= (2 - 1/ceil(M/2)) t < 2t`. The bound is attained by the Vandermonde measures (`upper.md`
Prop. 4.1, Thm 5.1). The same holds for full slices, and hence for any number of shells combined.

### 1.2 Independent numerics

`construct/check_barrier.py`, log `check_barrier.log`, uses my own Hilbert-function code
(`construct/clib.py`). That code works in the binomial basis, modulo `2097143` and `1999993`, and the two
primes always agreed.

* It computed `reg Q^M(p,n)` for:
  * `M = 3`, `t <= 10`;
  * `M = 4`, `t <= 8`;
  * `M = 5`, `t <= 6`;
  * `M = 6`, `t <= 5`.
* All 156 values satisfy `reg <= phi_M(p,n)`, which is `upper.md`'s Theorem 3.1.
* The maxima of `reg/t` are `3/2, 3/2, 5/3` for `M = 3, 4, 5`, equal to `2 - 1/ceil(M/2)`.
* For `M = 6` the maximum over `t <= 5` is `3/2`. The value `5/3` first occurs at `t = 9`.
* These are upper bounds by the modular-rank principle. Equality with the known product constructions
  of `upper.md` makes them exact.

### 1.3 The families the task proposed, one by one

Each item gives the exact reason it is capped, plus data where computed.

* **Alternating orbits with non-consecutive entries** (determinants with rows `z^a`, `zbar^b` chosen
  freely).
  * For distinct `v`, `Alt(delta_v)` has generating function `a_v(w) = V(w) s_lambda(w) w^mu`
    (bialternant formula).
  * `s_lambda(1, ..., 1)` is the number of semistandard tableaux, which is positive. So
    `ord = binom(M,2)` exactly.
  * The shell index is `t = sum |v_i|`, which is minimal for consecutive `v`.
  * So every alternating orbit has ratio at most the consecutive (Vandermonde) value
    `2 - 1/ceil(M/2)`. This is `(2s'+1)/(s'+1)` for `M = 2s'+1`.
* **Other isotypic components (hook, two-row Specht), Young-subgroup alternating sums, sums over set
  partitions, Schur or Jacobi–Trudi combinations (generalised alternants), hyperoctahedral
  compositions.**
  * All of these are measures on slices, so §1.1 caps them at `2 - 1/ceil(M/2)`.
  * Data: on `S^5(6,3)` I constructed independently the generalised alternant
    `c = sum_i gamma_i Alt(delta_{v_i})` with `gamma = (-2, 4, -2, -6)`
    (`construct/genalt.py`).
  * Its support has 480 points. Its order is **14**, certified by exact vanishing of all
    `binom(17,4)` moments of degree `<= 13`. Its ratio is `14/9 = 1.556 < 5/3`.
  * It beats products of alternants (13), as `upper.md` Remark 5.4 says. It does not beat the
    barrier.
* **Symmetric cut polytopes** `{sum k = s, k(B) <= U_|B|}` (`construct/scan_sym*.py`, log
  `scan_sym.log`). An exhaustive scan of profiles with `|U_j| <= 5, 4, 3` for `M = 3, 4, 5`
  (44, 70 and 37 distinct sets) computed exact `reg`. The maximal `reg / Lambda`, with `Lambda` the
  best symmetric half-separating adversary of §2.2, is `3/4, 3/4, 5/6`. This is
  `(2 - 1/ceil(M/2))/2`. Since `Lambda >= 2 w_bal`, it shows that no scanned set beats `upper.md`'s
  barrier.
* **Multi-shell combinations with `N`-polynomial coefficients.** These are forms on `X_M`. Their order
  is the minimum over slices `A_{s,D}` (`upper.md` Lemma 1.2), and §1.1 applies to each slice.
* **Coefficients depending on `N` non-polynomially.** On one circle,
  `F = sum_j lambda_j(N) G_j` with `G_j` homogeneous of degree `D_j` restricts to
  `sum_j lambda_j R^{D_j} g_j(theta)`.
  * Its order along the diagonal is the order of a measure on `A_{s, max D_j}`. Only one parity of
    `D` occurs for each `s`.
  * Its Liouville size is at least `max_j |lambda_j| R^{D_j}`.
  * So the effective degree is at least `max D_j`, and the slice barrier applies unchanged.

---------------------------------------------------------------------------------------------

## 2. Route (C): Liouville with forced Gaussian divisibility

### 2.1 The forced divisor is the tropical bound

Notation:

* `N >= 1`.
* `z_1, .., z_M` in `Z[i]` with `|z_j|^2 = N`, not necessarily distinct.
* `u_j = x_j + i y_j` and `v_j = x_j - i y_j`.
* A form of degree `D`, reduced on `X_M`, is
  `F = sum_{k in A} c_k m_k` with `m_k = (u_1 v_1)^{(D-|k|_1)/2} prod_j u_j^{k_j^+} v_j^{k_j^-}`.
  Here `A` is a finite subset of `{k : |k|_1 <= D, |k|_1 = D (mod 2)}` and `c_k` lies in `Z[i]`.
* `F(z)` means `F` evaluated at `u_j = z_j`, `v_j = conj(z_j)`.

**Theorem C1.**

* *(Split primes.)* Let `p = pi * conj(pi)` be a split prime with `p^e || N`. Put
  `a_j = v_pi(z_j)` (so `0 <= a_j <= e`) and `b_j = a_j - e/2`. Then
  ```text
  v_pi(F(z))       >=  min_{k in A} ( v_pi(c_k)       + eD/2 + <k,b> ),
  v_conjpi(F(z))   >=  min_{k in A} ( v_conjpi(c_k)   + eD/2 - <k,b> ).
  ```
* *(Inert primes and `1+i`.)* For an inert `q` with `q^{2f} || N`:
  `v_q(F(z)) >= min_k v_q(c_k) + fD`. For `1+i` with `2^g || N`:
  `v_{1+i}(F(z)) >= min_k v_{1+i}(c_k) + gD`.
* *(Consequence.)* Suppose every `c_k` is coprime to `N` and `F(z) != 0`. Then
  ```text
  Norm F(z)  >=  N^D * prod_{p split} p^{ - h_A(b_p) },     h_A(b) := max_{k in A} <k,b> - min_{k in A} <k,b>.   (2.1)
  ```

*Proof.*

* *Split primes.*
  * Conjugation swaps `pi` and `conj(pi)`, so `v_pi(conj z_j) = v_conjpi(z_j)`.
  * Also `v_pi(z_j) + v_conjpi(z_j) = v_p(N) = e`. Hence `v_pi(conj z_j) = e - a_j`.
  * Then
    ```text
    v_pi(m_k(z)) = e(D - |k|_1)/2 + sum_{k_j>0} k_j a_j + sum_{k_j<0} |k_j| (e - a_j) = eD/2 + <k, b>.
    ```
    The second equality uses `sum_j k_j a_j - (e/2) sum_j k_j = <k,b>` and
    `(e/2)(sum k^- - sum k^+) + e sum k^- = (e/2)|k|_1`.
  * The ultrametric inequality gives the first bound. The `conj(pi)` bound is symmetric.
* *Inert primes and `1+i`.*
  * `q` is fixed by conjugation up to a unit, so `v_q(z_j) = v_q(conj z_j) = f`. Every `m_k(z)` then has
    `v_q = fD`.
  * The same holds for `1+i`, since `1-i` is an associate of `1+i`.
* *Consequence.* `Norm F(z)` is the product over Gaussian primes of `Norm(prime)^{v}`.
  * At a split `p` the exponent sum is at least `eD + min_k <k,b> - max_k <k,b> = eD - h_A(b)`.
  * The inert part contributes at least `q^{2fD}` and the ramified part at least `2^{gD}`.
  * The product of these is `N^D prod_p p^{-h_A(b_p)}`.

`[]`

**Sharpness and special cases.**

* (i) For `F = u_i - u_j` the bound at `p` is `p^{e - |a_i - a_j|}` per conjugate pair. Their product
  is `Norm(g_ij)`, where `g_ij = gcd(z_i, z_j)`. So "products of pair differences inherit
  `prod g_ij`" is exactly C1 for root binomials.
* (ii) Data (`construct/gauss_check.py`, log `gauss_check.log`):
  * Points: all 96 lattice points of `N = 2 * 3^2 * 5^2 * 13 * 17 * 29`.
  * Forms: 352 random forms (Vandermonde alternants, random measures, and coefficients with random
    `pi`- and `conj(pi)`-content), at random tuples, `M <= 4`.
  * Every bound of C1 held, computed in exact Gaussian arithmetic.
  * The split-prime bound was attained with equality in 1322 of 1408 prime checks. The remaining
    cases are accidental residue cancellations.
* So C1 is the complete valuation-forced divisibility. Anything further depends on residues of the
  points, which is the larger-sieve input of `auxiliary.md`, not a forced divisor.

### 2.2 The criterion, and why the degree disappears

Fix `F` with content-free coefficients and `m := ord_Y F >= 1`. By `upper.md` Lemma 1.2 this is the
minimum over slices of the order of `c` restricted to the slice. Liouville (`tau_problem.md` step 1;
proved in `verify.md` §2) gives a constant `K = K(F)` with the following property.

* For points in an arc of length at most `C N^{1/4}`:
  `|F(z)| <= K R^D (C R^{-1/2})^m`.
* Explicitly, `K <= sum_k |c_k| (|k|_1)^m e^{|k|_1}/m!` for arcs of angle `<= 1`. This follows from
  Taylor's formula along the diagonal, slice by slice.
* So `Norm F(z) <= K^2 C^{2m} N^{D - m/2}`.

Combining with (2.1): if `F(z) != 0` then

```text
sum_{p split} log p * h_A(b_p)  >=  (m/2) log N - 2m log C - 2 log K.                       (2.2)
```

**The degree `D` has cancelled.** Normalise:

* `W = log N`;
* `omega_p = e_p log p / W`, so `sum omega_p <= 1`;
* `beta_p = (2/e_p) b_p`, which lies in `[-1,1]^M`.

Since `h_A(b_p) = (e_p/2) h_A(beta_p)`, (2.2) reads

```text
sum_p omega_p h_A(beta_p)  >=  m - (4m log C + 4 log K)/W.                                   (2.3)
```

**Theorem C2 (exact criterion).**

* (a) For every cluster tuple with `R >= R_0`, the forced-divisibility Liouville argument proves
  `F(z) = 0` exactly when its profile `(omega_p, beta_p)` violates (2.3).
* (b) The pair inequality is (2.3) for `F = u_i - u_j`, with `m = D = 1`, `K = 1` and
  `A = {e_i, e_j}`:
  ```text
  sum_p omega_p |beta_{p,i} - beta_{p,j}|  >=  1 - 4 log C / W,
  ```
  which is `d_ij >= W/2 - 2 log C` (the [P] inequality of `auxiliary.md`).
* (c) Let `A = {c, -c}` with `F = G - conj(G)` (order 1), where `G = prod z^{c+} conj(z)^{c-}` and
  `sum c = 0`. Then (2.3) is the project's character inequality. Since `h_A(b) = 2|<c,b>|`, it reads
  ```text
  sum_p e_p log p |<c, beta_p>|  >=  W/2 - O(log(C |c|_1)).
  ```
  The left side is the log-norm of `G / gcd(G, conj G)`. The inequality holds whenever `G` is
  nonreal.

*Proof.* (a) is (2.2) and (2.3) together with Theorem C1. (b) and (c) are direct substitutions;
`|z_i - z_j|` is at most the arc length. `[]`

**Lemma C2' (sign profiles suffice).** For every profile `(omega_p, beta_p)` there is a measure `nu`
on `{+-1}^M` with the following properties.

* `nu` has total mass `sum_p omega_p`.
* `sum_eps nu(eps) |eps_i - eps_j| = sum_p omega_p |beta_{p,i} - beta_{p,j}|` for all `i, j`.
* `sum_eps nu(eps) h_A(eps) >= sum_p omega_p h_A(beta_p)`.

*Proof.*

* For `beta` in `[-1,1]^M` and `tau` uniform on `[-1,1]`, put `eps_tau = 2 * 1[beta > tau] - 1`.
* Coordinatewise `E_tau eps_tau = beta`. Also `E_tau |eps_{tau,i} - eps_{tau,j}| = |beta_i - beta_j|`.
* `h_A` is sublinear (a max minus a min of linear forms). Hence
  `h_A(beta) = h_A(E eps_tau) <= E_tau h_A(eps_tau)`.
* Let `nu` be the sum over `p` of `omega_p` times the law of `eps_tau` for `beta_p`. `[]`

So, up to `O(1/W)` errors, the method with `F` works only if `m > Lambda(A)`, where

```text
Lambda(A) := max { sum_eps nu(eps) h_A(eps) :  nu >= 0,  sum nu <= 1,  sum_eps nu(eps) |eps_i - eps_j| >= 1  for all i < j }.
```

`Lambda` is a linear programme. Conversely, `m > Lambda(A)` would make (2.3) fail for every admissible
profile once `W` is large, by continuity of the LP value in its right-hand side.

In cut language: on one slice, `h_A(eps) = 2 w_P(A)` with `P = {eps = +1}`, and the constraint says
that every pair is separated with `nu`-probability at least `1/2`.

### 2.3 What ratio is needed

Without divisibility the requirement is `ord > 2D`. With the tropical divisor it is

```text
ord(c)  >  Lambda(supp c)  >=  Omega(supp c) := 2^{-M} sum_{eps in {+-1}^M} h_{supp c}(eps).          (2.4)
```

The second inequality holds because the uniform measure on `{+-1}^M` is feasible:
`E|eps_i - eps_j| = 1`. On `A_{s,D}` every width satisfies `h(eps) <= 2D`, so `Omega <= (2 - 2^{2-M}) D`.
So divisibility genuinely lowers the requirement: the target moves from `2D` down to the cube mean
width. Theorem B shows that it is still never met.

### 2.4 Theorem B (mean-width bound)

For finite `S` in `Z^n` and `P` a subset of `[n]`, let `x(P) = sum_{i in P} x_i` and
`w_P(S) = max_{x in S} x(P) - min_{x in S} x(P)`. In particular `w_empty = 0`. Put

```text
MW_n(S) := 2^{1-n} sum_{P subset [n]} w_P(S) = 2 E_P w_P(S)        (P uniform among the 2^n subsets).
```

**Theorem B.** For every `n >= 0` and every finite nonempty `S` in `Z^n`: `reg(S) <= MW_n(S)`.
Equality holds for every box `prod_i [0, a_i]`, for which `reg = sum a_i = MW_n`.

*Proof.* Induction on `n`.

* *Base `n = 0`.* `S` is a point, so `reg = 0 = MW_0`.
* *Setup for `n >= 1`.* Slice by `x_1 = c`. Each nonempty slice `S_c` is an affine copy of a set
  `S'_c` in `Z^{n-1}` (coordinates `2..n`), so `reg(S_c) = reg(S'_c)`. By induction,
  ```text
  reg(S_c) <= g_c := MW_{n-1}(S'_c) = 2 E_{P'} l_c(P'),
  ```
  where `P'` is uniform among subsets of `{2..n}` and `l_c(P') = w_{P'}(S_c)`.
* *Splitting `MW_n`.* Splitting on whether `1` lies in `P`,
  ```text
  K := MW_n(S) = E_{P'} [ w_{P'}(S) + w_{P' cup {1}}(S) ].
  ```
* *Key inequality.* Let `c, c'` index nonempty slices, and fix `P'`.
  * Let `[a_1, a_2]` be the range of `x(P')` on `S_c`, and `[b_1, b_2]` its range on `S_{c'}`.
  * Since `x(P' cup {1}) = x_1 + x(P')`, we have `w_{P' cup 1}(S) >= (c' + b_2) - (c + a_1)`.
  * Also `w_{P'}(S) >= a_2 - b_1`.
  * Adding,
    ```text
    w_{P'}(S) + w_{P' cup 1}(S) >= (c' - c) + l_c(P') + l_{c'}(P').
    ```
  * Averaging over `P'`:
    ```text
    K >= (c' - c) + (g_c + g_{c'})/2.                 (2.5)
    ```
    With `c = c'`, this gives `g_c <= K` for every slice.
* *Counting.* For real `v`, let `C_v` be the set of slices with `g_c >= v`. If `C_v` has two elements,
  (2.5) applied to its extreme elements gives `max C_v - min C_v <= K - v`. Hence
  `N(v) := |C_v| <= K - v + 1`. If `|C_v| = 1` the same holds, because `v <= g_c <= K`.
* *Conclusion.* Order the slices by decreasing `g`. The peeling lemma (§1.1, Lemma 2.1) gives
  ```text
  reg(S) <= max_i ( reg(S_(i)) + i - 1 ) <= max_i ( g_(i) + i - 1 ).
  ```
  Since `N(g_(i)) >= i`, the counting step gives `g_(i) + i - 1 <= K`.

`[]`

**Homogeneous form (Corollary B1).** Let `A` be finite in the slice `{k in Z^M : sum k = s}`. Then
every `c != 0` on `A` has

```text
ord(c) <= reg(A) <= Omega(A) = 2^{-M} sum_eps h_A(eps) = 2^{1-M} sum_{P subset [M]} w_P(A).
```

*Proof.*

* Drop `k_M`, an affine bijection onto its image `x(A)` in `Z^{M-1}`.
* For `P` not containing `M`, `w_P(A) = w_P(x(A))`.
* For `P` containing `M`, `k(P) = s - k(P^c)`, so `w_P(A) = w_{P^c}(x(A))`.
* Hence `2^{1-M} sum_{P subset [M]} w_P(A) = MW_{M-1}(x(A))`.
* Finally `h_A(eps) = 2 w_{eps=+1}(A)` on a slice, and Lemma 1.1 of `upper.md` gives
  `ord <= reg`. `[]`

**Corollary B2 (second barrier for (B)).** For every `M`: `tau*(M) <= 2 - 2^{2-M} < 2`.

*Proof.* On `A_{s,D}`, every `k` satisfies `k(P) <= (D+s)/2` and `k(P) >= (s-D)/2`, so
`w_P <= D`. The two trivial cuts have width `0`. Hence
`Omega(A_{s,D}) <= 2(1 - 2^{1-M}) D`. `[]`

This is independent of `upper.md`. It is weaker for `M >= 4`, equal at `M = 2, 3`, and it coincides
with the "crude" constant mentioned in `upper.md` §3.

**Remarks.**

* **(Equality cases.)** Every product of root binomials `w_i - w_j` attains equality in Corollary B1.
  Its support is the Minkowski sum of the segments `[e_i, e_j]`, so `Omega = sum_edges 1`, which is the
  order. This covers Vandermonde measures and products of Vandermondes over set partitions.
  * Data (`construct/tight.py`, log `tight.log`): 1662 measures with `ord = Omega(supp)`. These are
    generic elements of exact rational kernels on random sets, `n <= 4`.
  * **All 1662** are, up to a monomial, products of root binomials `(x_i - 1)` and `(x_i - x_j)`.
    This was found by exact division.
  * **Conjecture E:** equality in Theorem B for a measure holds only for root-binomial products.
* **(Data.)** `construct/verify_B.py` (log `verify_B.log`) checked 6400 random sets in dimension
  `n <= 4`. In every case `reg <= (bound produced by the proof) <= MW`, and equality held in about 28%
  of them. `conjB.py`, `cbal.py` and `quantify.py` add further checks.
* **(Novelty.)**
  * I grepped `research/docs` for `mean width`, `Rademacher`, `regularity index`,
    `Hilbert function`, `Newton polytope`, `cut-width` and `annihilator`.
  * The closest note is `sl2_invariant_average_cut_width.md`. It is an additive analogue: the average
    cut half-width of the **Newton support of a translation-invariant (`SL(2)`) polynomial** of degree
    `E` is at least `E/4`, with equality for bracket products. Its proof uses a different grid lemma.
  * That note bounds the Newton support of a polynomial. Theorem B bounds the support of a vanishing
    measure in terms of its order of vanishing. I did not find Theorem B there.
  * The group-algebra "annihilator capacity" framework
    (`gaussian_annihilator_capacity.md`, `annihilator_capacity_criticality.md`) works dually, with
    linear combinations of the cluster's points. It found the same endpoint-criticality for its own
    width functional.

### 2.5 The cube profile defeats every auxiliary polynomial

**Corollary C3.** Fix `M`, a support `A` (any number of slices), and `c` of order `m >= 1` on it. Take
`N = prod_{eps in {+-1}^M} p_eps`, with `2^M` distinct primes `p_eps = 1 (mod 4)` in an interval
`[X, (1+theta)X]`. Assign to `p_eps` the sign vector `eps`, and
let

```text
z_j = prod_eps pi_eps^{[eps_j = +1]} conj(pi_eps)^{[eps_j = -1]}.
```

Then the following hold.

* (i) `|z_j|^2 = N`, and for all `i != j`
  ```text
  |d_ij - W/2| <= 2^{M-1} log(1+theta).
  ```
* (ii) Let `F` be any form on `X_M` with `ord_Y F = m` and content-free coefficients. Then
  ```text
  | sum_p log p * h_A(b_p) - (W/2) Omega(A) | <= 2^{M-1} log(1+theta) max_eps h_A(eps),     Omega(A) >= m.
  ```
* Here `W = log N` tends to infinity with `X`, and `theta` may tend to `0` (for example
  `theta = 1/log X`, which still leaves many primes `= 1 (mod 4)` in the interval).
  So (2.2) is satisfied with slack `(W/2)(Omega(A) - m) >= 0`, up to `o(1)` and the constant
  `2m log C + 2 log K`.

*Proof.*

* (i) `e_p = 1` and `b_p = eps_p / 2`. Exactly half of the sign vectors separate `i` and `j`.
  Moreover each `log p_eps` lies within `log(1+theta)` of the mean of the logarithms. Both errors are
  half a sum of `2^M` such deviations.
* (ii) Use `h_A(b_p) = h_A(eps_p)/2`. For `Omega(A) >= m`: if `s_0` is a slice where `c` has order
  `m`, then `Omega(A) >= Omega(A_{s_0}) >= ord(c_{s_0}) = m`. The first inequality holds because `h` is
  monotone in `A`; the second is Corollary B1. `[]`

**Consequences.**

* **(a) Exponent level.** No auxiliary form, of any degree and with any number of shells, can force
  vanishing on clusters whose valuation profile is (close to) the cube profile.
  * At the same time these profiles satisfy every pair inequality.
  * By Theorem C2(c) they also satisfy the character inequalities.
  * In LP terms, `Lambda(A) >= Omega(A) >= ord(c)` for every `c`. So the strict inequality (2.4)
    needed by the method is never achieved.
* **(b) Constant level.**
  * *Non-tight `F`, `Omega(A) > m`.* Let `C > 0` be arbitrary, even `C < 1`. Mix the cube profile with
    the uniform-balanced profile `mu_bal` (`|P| = floor(M/2)`) at weight
    `eta = 4 log(1/C) / (W (2 pi_bal - 1))`. Here `pi_bal = ceil(M/2)/(2 ceil(M/2) - 1)` is the
    probability that a balanced cut separates a pair.
    - The pair inequalities with `C` then hold.
    - The payoff stays `>= Omega(A) - eta Omega(A)`.
    - So (2.3) holds for `W >= W_0(F, C)`.
  * *Tight `F`, `Omega(A) = m`.* This is the equality case of Theorem B: the Vandermonde, and products
    of pair differences.
    - The cube profile meets (2.3) only up to a bounded constant.
    - If Conjecture E holds, tight `F` are root-binomial products. Their inequality is literally the
      sum of the pair inequalities, which every admissible profile satisfies.
    - If Conjecture C_bal (§2.8) holds, `mu_bal` itself satisfies (2.3) with slack proportional to `W`.
    - Either conjecture therefore makes the barrier robust for every `C > 0`.

### 2.6 Coefficients with arithmetic content (`N`-dependent coefficients)

Now let the coefficients `c_k` lie in `Z[i]` and depend on `N` arbitrarily. They may then carry prime
content at the primes of `N`.

**Proposition C4.**

* *Hypotheses.*
  * Take `F` as in §2.1 with arbitrary `c_k` in `Z[i]`, support `A`, and `ord_Y F = m`.
  * Let `L_k = log Norm(c_k)`.
  * Let `N` be squarefree, with prime factors split and grouped into blocks. Each block has `2^M`
    primes in an interval `[X_t, (1+theta)X_t]`. The pair inequalities then hold up to
    `2^{M-1} log(1+theta)` per block, as in Corollary C3(i).
  * The expectation bound below is exact and needs no such condition.
  * Within each block, assign the sign vectors by independent uniformly random bijections, and form
    `z^sigma` as in Corollary C3.
* *Conclusion.* Put `alpha_k(p) = v_pi(c_k)` and `alpha'_k(p) = v_conjpi(c_k)`, and
  ```text
  T_p(b) = min_k (alpha_k(p) + <k,b>) + min_k (alpha'_k(p) - <k,b>),
  ```
  which is the split-prime forced exponent minus `eD = D` (Theorem C1, `e = 1`). Then
  ```text
  E_sigma  sum_p log p * T_p(b_p)  <=  max_k L_k - (W/2) Omega(A)  <=  max_k L_k - (m/2) W.
  ```
* *Consequence.* Liouville gives `Norm F(z) <= K_0^2 |A|^2 e^{max L_k} C^{2m} N^{D-m/2}`, with `K_0`
  depending only on `(A, m)`. So some assignment satisfies the forced-divisibility inequality up to a
  constant. Content gains nothing proportional to `log N`.

*Proof.*

* *Per-prime inequality.* Let `u, v` be argmin and argmax of `<k, eps>` on `A`, chosen by a fixed rule
  that depends only on `eps`. Then
  ```text
  T_p(b) + T_p(-b) <= (alpha_u + alpha'_u + alpha_v + alpha'_v)(p) - 2(<v,b> - <u,b>).
  ```
  For the bound on `T_p(b)`, use `u` in the first min and `v` in the second. For `T_p(-b)`, swap the
  roles.
* *Expectation at one prime.* Each prime's sign vector is uniform and its law is symmetric under
  `eps -> -eps`. With `b = eps/2`,
  ```text
  E T_p <= (1/2) E_eps [ beta_{u(eps)}(p) + beta_{v(eps)}(p) ] - (1/2) E_eps h_A(eps),
  ```
  where `beta_k = alpha_k + alpha'_k`.
* *Summing over `p`.* Multiply by `log p` and sum. By linearity,
  `sum_p log p E_eps beta_{u(eps)}(p) = E_eps sum_p log p beta_{u(eps)}(p) <= max_k L_k`. The
  same holds for `v`.
* Finally `(1/2) sum_p log p E h_A = (W/2) Omega(A) >= (W/2) m`.
* Data: the per-prime inequality held in 20000 random exact instances (`quantify.py`). `[]`

### 2.7 The numbers (exact; `construct/quantify.py`, log `quantify.log`)

The method needs `ord > Lambda` (the LP of §2.2). Its necessary consequence is `ord > Omega`, the
cube payoff.

| measure | ord | `Omega` (cube) | ord/`Omega` | `2E_bal w` | `Lambda` (exact LP) | ord/`Lambda` |
|---|---|---|---|---|---|---|
| Vandermonde `M=3`, `(-1,0,1)` | 3 | 3 | **1** | 4 | 4 | 3/4 |
| Vandermonde `M=4`, `(-1,0,1,2)` | 6 | 6 | **1** | 8 | 8 | 3/4 |
| Vandermonde `M=5`, `(-2..2)` | 10 | 10 | **1** | 12 | 12 | 5/6 |
| Vandermonde `M=6`, `(-2..3)` | 15 | 15 | **1** | 18 | 18 | 5/6 |
| Vandermonde `M=7`, `(-3..3)` | 21 | 21 | **1** | 24 | 24 | 7/8 |
| generalised alternant on `S^5(6,3)` | 14 | 65/4 | 0.862 | 18 | 18 | 7/9 |
| full slice `Q^4(2,2)` (`reg = 6`) | 6 | 7 | 0.857 | 8 | 8 | 3/4 |
| full slice `Q^5(3,3)` (`reg = 10`) | 10 | 45/4 | 0.889 | 12 | 12 | 5/6 |
| full slice `Q^6(2,2)` (`reg = 6`) | 6 | 31/4 | 0.774 | 8 | 8 | 3/4 |
| spike profile `M=7`, `n=8` (peeling bound) | `<= 111` | 287/2 | `<= 0.774` | 128 | 152 | `<= 0.730` |

How to read the table.

* Forced divisibility moves the Vandermonde from `ord/(2t) = (2s'+1)/(2s'+2) < 1` to
  `ord/Omega = 1` exactly. This is the Cilleruelo–Córdoba endpoint, made critical.
* Against the best admissible profile, the Vandermonde ratio is `1 - 1/(2 ceil(M/2))`. This equals
  `tau(M)/2` of `upper.md`.
* Nothing exceeds `1` (Theorem B).
* The "spike" polytope `{sum x = 0, x(B) <= n max(|B|, M - |B|)}` is narrow on balanced cuts and wide
  on unbalanced ones. It is the natural candidate against balanced adversaries. Even there the
  rigorous peeling bound for `reg` is far below `Omega` and `Lambda`. Peeling bounds give
  `reg/n = 13, 13.5, 14, 13.75, 13.8, 14, 13.86, 13.88` for `n = 1..8` (`spike_peel.log`).

### 2.8 Further evidence, the balanced conjecture, and a failed shortcut

* **Conjecture C_bal.** For every `c`,
  ```text
  ord(c) * pi_bal <= E_bal w_P(supp c),
  ```
  where `P` is a uniform balanced cut and `pi_bal = ceil(M/2)/(2 ceil(M/2) - 1)`.
  * It is tight for root-binomial products.
  * It contains `upper.md`'s barrier, since `w_P(Q^M(p,n)) = p + n` for every cut.
  * It would make `mu_bal` (pair separation `pi_bal > 1/2`) defeat every `F` with slack proportional
    to `W`. It would also give `ord/Lambda <= 1 - 1/(2 ceil(M/2))`.
  * Data:
    - random exact sets, `M = 4, 5, 6` (`cbal.log`): no violation, minimum margin 0 at root
      segments;
    - symmetric cut polytopes, `M = 6..9` (`cbal_sym.log`): rigorous peeling upper bounds for `reg`,
      no violation, maximum ratio `1.000` at alternant-type profiles.
* **A shortcut that fails.** My first proof attempt for Theorem B was a different induction:
  restrict to `w_i = w_j` if `f` does not vanish there, and factor `w_i - w_j` otherwise.
  * It needs a "good pair" with `E[w_P ; P separates i,j] >= E[w_P ; P does not]`.
  * After symmetrisation, the summed form of this condition is
    `sum_b binom(M,b) psi(b) ((M - 2b)^2 - M) <= 0`, where `psi(b)` is the average width of cuts of
    size `b`.
  * This holds for all width profiles when `M <= 6`. The proof uses the seminorm identities
    `2 e_i = (e_i + e_j) + (e_i + e_l) - (e_j + e_l)` and
    `2 e_i = 1_{ijk} + 1_{ilm} - 1_{jklm}`.
  * It fails at `M = 7` for the spike profile `psi = (6, 5, 4)`: the sum is `+252`.
  * Random sets never violate it (`bm_gpl.log`). The slicing proof of §2.4 avoids the issue
    entirely.
  * "Balance monotonicity" (balanced cuts are widest on average) is false: 312 of 3000 random
    small sets violate it at `M = 4`, and 37 of 3000 at `M = 7` (`bm_gpl.log`). The summed good-pair
    inequality held on all random sets tested, for `M = 7` and `M = 9`.

---------------------------------------------------------------------------------------------

## 3. What this leaves, and what the theorem would need

1. **Single-place auxiliary polynomials are exhausted, with or without forced divisibility.**
   * Without divisibility the exponent is `(k-1)/(2k-1) < 1/2`.
   * With the tropical divisor the method sits exactly at the endpoint. Its inequalities are all
     satisfied by the cube profile, which also satisfies all pair and character inequalities.
2. **The question is now realizability, not auxiliary constructions.** A proof must show that
   cube-like valuation profiles cannot occur on short arcs. This needs an input that ties the
   archimedean phases to the prime splitting.
   * In the notes the relevant model is the complete Walsh/Hadamard profile, with columns equal to all
     sign vectors. The repeated-Walsh growth results (items 311–335) use phase information beyond
     valuations.
   * The Paley obstruction (item 334) and `annihilator_capacity_criticality.md` reach the same
     endpoint-criticality from the character and group-algebra sides.
3. **Inputs not covered by the barrier.**
   * Residue information (congruences of units, the larger sieve). It is not forced divisibility.
   * Exact phase relations, e.g. using `arg pi_p`.
   * Several auxiliary forms jointly in a determinant. `auxiliary.md` Prop. A reduces these to
     Vandermonde times a symmetric factor.
   * Constant-level effects in the equality case of Theorem B. These are closed under Conjecture E or
     C_bal.

---------------------------------------------------------------------------------------------

## 4. Answers to the task's explicit items

* **Construct `c` with `ord(c) > 2t`.** Impossible, for every `M` and every shell: §1.1, and
  Corollary B2 independently.
* **Isotypic components, products and sums of Vandermondes, Schur or Jacobi–Trudi, other character
  determinants, hyperoctahedral compositions.** All capped at `2 - 1/ceil(M/2)` (§1.3). Exact data: a
  generalised alternant with order 14 and `t = 9`; symmetric cut polytopes for `M <= 5`.
* **Multi-shell with `N`-dependent coefficients.**
  * Polynomial in `N`: these are forms, so the slice barrier applies.
  * Arbitrary: effective degree, §1.3.
  * With arithmetic content and forced divisibility: Proposition C4.
* **Forced divisibility, "ratio below 2 plus forced divisor".**
  * Quantified exactly in §2.2–2.3: `D` cancels, and the needed condition is
    `ord(c) > Lambda(supp c) >= Omega(supp c)`.
  * Achievable: `ord = Omega` exactly, by the Vandermonde and by every product of pair differences.
  * Never `ord > Omega`, by Theorem B.
  * So the forced divisor supplies exactly the missing `R^{D/(2 ceil(M/2))}` for the Vandermonde, and
    nothing beyond the boundary for any `F`.

---------------------------------------------------------------------------------------------

## 5. Files (`scratchpad/tau/construct/`)

| file | content |
|---|---|
| `clib.py` | Hilbert function and `reg` mod two primes (binomial basis); cut widths; exact rational simplex for `Lambda`; symmetric adversary |
| `gens.py` | point sets: `Q^M(p,n)`, orbits, cut polytopes |
| `check_barrier.py` → `check_barrier.log` | independent check of `upper.md` Thm 3.1 values |
| `verify_B.py` → `verify_B.log` | Theorem B: random sets, the proof's recursive bound, boxes, root products |
| `tight.py` → `tight.log` | equality cases of Theorem B: 1662 found, all root-binomial products |
| `gauss_check.py` → `gauss_check.log` | Theorem C1 on actual Gaussian integers (exact) |
| `genalt.py` | generalised alternant on `S^5(6,3)`: order 14, `Omega = 65/4` |
| `quantify.py` → `quantify.log` | the table of §2.7; Proposition C4 per-prime check |
| `peeldp.py`, `spike_peel.log` | rigorous sorted-peeling upper bounds for symmetric cut polytopes |
| `scan_sym.py`, `scan_sym2.py` → `scan_sym.log` | exhaustive symmetric cut-polytope scans `M <= 5` |
| `scan_peel.py`, `cbal_sym.py` → `scan_peel.log`, `cbal_sym.log` | symmetric profiles `M <= 9` (peeling bounds) |
| `cbal.py` → `cbal.log`; `conjB.py` | Conjectures C_bal and B on random exact sets |
| `bm_test.py`, `sum_test.py` → `bm_gpl.log` | balance monotonicity (false) and good-pair statistics |
| `spike.py` | spike polytope lattice points and exact `reg` for `M = 5, 6`, `n <= 2` |
