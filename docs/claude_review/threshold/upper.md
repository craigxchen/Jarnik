# The auxiliary-polynomial threshold is exactly `tau(M) = 2 - 1/ceil(M/2)`, so it is always below 2

Task: prove `tau(M) <= 2` (a barrier), or compute `tau(M)` exactly. Problem statement:
`../tau_problem.md`.

## 0. Verdict

**Barrier theorem, with the exact value.** For every `M >= 2`,

```text
tau(M)  =  tau*(M)  =  2 - 1/ceil(M/2)          (= 1, 3/2, 3/2, 5/3, 5/3, 7/4, 7/4, ... for M = 2, 3, 4, ...)
```

The two quantities are these.

* `tau(M)` is the lead's shell quantity from `tau_problem.md`.
* `tau*(M)` is the quantity the reduction actually needs: the pseudo-effective threshold of
  `pi^*H - tE` on `Bl_Y X_M`. It is defined through full slices of the cross-polytope (§1).

The supremum is attained, by the Vandermonde (alternating-orbit) measures. Consequences:

1. For every `M` and every form `F` of degree `D` with `F` not identically zero on `X_M`,
   `ord_Y F <= (2 - 1/ceil(M/2)) D < 2D`. So step 1 of the claimed reduction (Liouville, which
   needs `m > 2D`) can never be applied. This holds whether `F` is built from one shell, several
   shells, or anything else.
2. The Ru–Vojta constant of `conditional.md` §6 satisfies
   `beta_M = beta(H,E) <= tau*(M) <= 2 - 1/ceil(M/2) < 2` for every `M`. The open sub-question
   "is `beta(H,E) > 2` on some `Bl_Y X_M`?" therefore has the answer **no**. The twisted
   constants `b + beta(H - bE, E)` are also `<= tau*(M) < 2`.
3. The exponent reached by single-place auxiliary polynomials on `X_M` is exactly
   `alpha_M = 1 - 1/tau*(M) = (k-1)/(2k-1)` with `k = ceil(M/2)`. This is `1/3, 2/5, 3/7, ...`,
   the Cilleruelo–Córdoba sequence `1/2 - 1/(4k-2)`. The method cannot pass this value for any
   fixed `M`. It reaches `1/2` only in the limit, never at it.

**Finer exact result.** For each pair `(p,n)` (the "positive mass" and "negative mass" of a
slice; see §1) the maximal vanishing order is bounded by the explicit concave function

```text
phi_M(p,n) = (M-1) * min_{1<=j<=M-1} ( p/(M-j) + n/j ).
```

This bound is asymptotically sharp for every `(p,n)`. It is attained exactly at the
`M` "Vandermonde points" `u_j = (T_{M-j}, T_{j-1})`, where `T_x = x(x+1)/2`. For `M = 3, 4`
the maximal order is known in closed form for every `(p,n)` (Theorem 5.3).

**A correction to the claimed reduction.** "The N-weights separate shells" is false. On `X_M`
the basis element `N^{(D-t)/2} [z^k]` (reduced monomial `[z^k] = prod u_j^{k_j^+} v_j^{k_j^-}`)
restricts near `Y` to `u^{(D+s)/2} v^{(D-s)/2} w^k`. Equivalently, in Laurent form it is
`N^{(D-s)/2} z^k`, not `N^{(D-t)/2} z^k`. The weight depends on `s` only, not on `t`. So shells with the same `s` interact. An explicit example on
`X_2` is given in Remark 1.3. The correct per-`s` support is the full slice
`A_{s,D} = {k : sum k = s, |k|_1 <= D}`. This is the set used in `conditional.md` Prop. 6.2,
and the referee confirmed it independently (`referee_conditional.md`, check R1). Since every
shell lies in a slice, `tau(M) <= tau*(M)`. Both equal `2 - 1/ceil(M/2)`, so the error does not
change any value. But the pseudo-effective threshold is `tau*`, not `tau` by definition.

**Method.** The tool is (iii) from the task list: explicit interpolation, done by a hyperplane
peeling lemma. Cut the slice `k_1 = c`, sort the slices, and induct on `M`. The induction only
closes when it carries the whole family of `M-1` linear bounds `l_j`, not just the symmetric
bound. The lower bound is the Vandermonde/alternant divisibility argument together with products.
By-products:

* the crude version of the induction already gives `tau*(M) <= 2 - 2^{2-M}`;
* lines in either factor of the Cremona-graph model give the BDPP-type curve bounds `l_1` and
  `l_{M-1}` (Remark 3.3).

**Novelty check.** I grepped `research/docs/` and `growth/` for "Hilbert function",
"interpolat", "Seshadri", "Cremona", "root polytope", "pseudo", "threshold", "Bl_Y" and
"tau(M)". `conditional.md` (§6) and `referee_checks/tau_threshold.py` define and compute the
same threshold numerically, up to `M=5, N=6`. They prove no upper bound, and they leave
`beta > 2` open. I found no peeling or interpolation bound in the notes. I did no literature
search. The peeling lemma (Lemma 2.1) is a standard fact about regularity indices of finite
point sets.

Everything labelled Theorem, Lemma, Proposition or Corollary is proved below. The numerics in §7
are evidence and consistency checks only.

---------------------------------------------------------------------------------------------

## 1. Setting, the correct reduction, and the duality `order = regularity`

**Notation.**

* `M >= 1`.
* For integers `p, n >= 0` let
  ```text
  Q^M(p,n) = { k in Z^M : sum_i k_i = p - n,  sum_i |k_i| <= p + n }        ("slice", or "ball")
  S^M(p,n) = { k in Z^M : sum_i k_i^+ = p,     sum_i k_i^- = n }             ("shell")
  ```
  Here `k^+ = max(k,0)` and `k^- = max(-k,0)`. Then `S^M(p,n)` is a subset of `Q^M(p,n)`, and
  `Q^M(p,n)` is the disjoint union of the shells `S^M(p-i, n-i)` for `0 <= i <= min(p,n)`.
* In the problem's notation: `s = p - n`, `t = p + n`, and `S_{s,t} = S^M(p,n)`.
* The full slice of degree `D` and weight `s` is `A_{s,D} = Q^M((D+s)/2, (D-s)/2)`. This is
  nonempty iff `|s| <= D` and `s = D (mod 2)`.
* For a finite nonempty `A` in `Q^M`, the **regularity index** is
  `reg(A) = min{ d >= 0 : every function A -> Q is the restriction of a polynomial of degree <= d }`.
  Lagrange interpolation gives `reg(A) <= |A| - 1`.
* Two elementary properties of `reg`:
  - it is monotone: `A' subset A` implies `reg(A') <= reg(A)`, because one can extend by 0;
  - it is invariant under affine bijections, which preserve degree.

  The same `reg` is obtained over any field of characteristic 0, because ranks of rational
  matrices do not change under field extension.
* For `c : A -> C` not identically zero, `ord(c)` is the largest `m` such that
  `sum_k c_k P(k) = 0` for all polynomials `P` of degree `< m`. This is the problem's
  definition.

**Lemma 1.1 (order and regularity).**

* (a) `ord(c)` equals the order of vanishing at `w = (1, .., 1)` of the Laurent polynomial
  `f_c(w) = sum_k c_k w^k` on `(C^*)^M`.
* (b) `max{ ord(c) : c != 0 supported on A } = reg(A)`.

*Proof.*

* (a) Put `w = e^y`. Then `f_c(e^y) = sum_{d>=0} (1/d!) sum_k c_k (k.y)^d`. The map `exp` is a
  local biholomorphism at `0`, so the order at `w = 1` is the least `d` such that
  `sum_k c_k (k.y)^d` is not the zero polynomial in `y`. Expanding `(k.y)^d` shows that this
  polynomial vanishes identically iff `sum_k c_k k^alpha = 0` for all `|alpha| = d`.
* (b) Suppose `ord(c) >= m` and `c != 0`. Then `c` annihilates `Poly_{<m}|_A`, so
  `Poly_{<m}|_A` is not all of `Q^A`, and hence `reg(A) >= m`. Conversely, if
  `d = reg(A) >= 1`, then `Poly_{<=d-1}|_A` is a proper subspace of `Q^A`. A nonzero `c`
  orthogonal to it has `ord(c) >= d`. If `reg(A) = 0`, any `c != 0` has `ord(c) >= 0`. `[]`

**Lemma 1.2 (correct reduction to slices).** Work over `Q(i)` with `u_j = x_j + i y_j` and
`v_j = x_j - i y_j`, so that `X_M = {u_1 v_1 = ... = u_M v_M}` and
`Y = {u_j = u_1, v_j = v_1}`. Let `F` be a form of degree `D`.

1. On `X_M`, `F` has a unique expansion `F = sum_k c_k N^{(D - |k|_1)/2} [z^k]`. Here
   `[z^k] = prod_{k_j>0} u_j^{k_j} prod_{k_j<0} v_j^{-k_j}` and `N = u_1 v_1`. The index `k`
   runs over `|k|_1 <= D` with `|k|_1 = D (mod 2)`. Uniqueness holds because `X_M` is
   projectively normal and these monomials form a basis of `H^0(O(D))` (`conditional.md`
   Prop. 6.2, step 1).
2. Near a point of `Y` with `uv != 0`, write `u_j = u w_j`, `v_j = v/w_j`, `w_1 = 1`. Then
   ```text
   F = sum_s u^{(D+s)/2} v^{(D-s)/2} G_s(w),     G_s(w) = sum_{k in A_{s,D}} c_k w^k.
   ```
3. Hence `ord_Y F = min_{s : G_s != 0} ord_{w=1} G_s`, and
   ```text
   max { ord_Y F : deg F = D, F not 0 on X_M }  =  max_s reg(A_{s,D}).
   ```
4. The pseudo-effective threshold of `pi^*H - tE` on `Bl_Y X_M` is
   ```text
   tau*(M) := sup_{D,s} reg(A_{s,D})/D  =  sup_{p+n>=1} reg Q^M(p,n)/(p+n).
   ```
   The two expressions agree via `D = p+n` and `s = p-n`.

*Proof.*

* For item 2, note `N^{(D-t)/2}[z^k] = (uv)^{(D-t)/2} u^{|k^+|} v^{|k^-|} w^k`. Since
  `|k^+| = (t+s)/2` and `|k^-| = (t-s)/2`, the exponents of `u` and `v` are `(D+s)/2` and
  `(D-s)/2`. These depend on `s = sum k` only.
* For item 3, use `(v/u, w_2, .., w_M)` as local coordinates near a point of `Y` with
  `uv != 0`. Then `Y = {w = 1}`. The Taylor coefficient of `(w-1)^alpha` in `F/u^D` is
  `sum_s rho^{(D-s)/2} (d^alpha G_s)(1)/alpha!`, with `rho = v/u`. The powers of `rho` are
  distinct for distinct `s`, so this coefficient vanishes identically in `rho` iff every
  `(d^alpha G_s)(1)` does. Since `G_s` is homogeneous of degree `s` in `w`, its order at
  `(1, .., 1)` on `(C^*)^M` equals the order of `G_s(1, w_2, .., w_M)` at `1`.
* The maximum in item 3 is attained by `F` with a single nonzero `G_s`. Apply Lemma 1.1.
* For item 4, `Y` lies in the smooth locus of `X_M` (`conditional.md` Lemma 2.1). So sections
  of `D pi^*H - mE` are the forms of degree `D` vanishing to order `>= m` along `Y`. Max-orders
  are superadditive in `D`, because products of sections are sections, so the limit equals the
  supremum. `[]`

The lead's quantity is `tau(M) = sup reg(S_{s,t})/t`. By Lemma 1.1 and monotonicity of `reg`,
`tau(M) <= tau*(M)`.

**Remark 1.3 (the shells do not decouple).** Take `M = 2` and `D = 2`, and let
`F = u_1 v_2 + v_1 u_2 - 2 u_1 v_1` (in the circle model,
`z_1 zbar_2 + zbar_1 z_2 - 2|z_1|^2`).

* On `X_2`, `F = -|z_1 - z_2|^2`, so `ord_Y F = 2`.
* Its `t = 2` shell part `u_1 v_2 + v_1 u_2` and its `t = 0` shell part `-2N` do not vanish
  on `Y` at all. On `Y` they equal `2N` and `-2N`.
* So the order of `F` is not the minimum of the orders of its shell components. Here
  `G_0(w) = w_1/w_2 + w_2/w_1 - 2`, which has order 2. (Check (D) in `certify.py`.)

**Remark 1.4 (the `C^*`-weights).** The diagonal torus acts transitively on `Y minus {uv = 0}`.
So the generic order along `Y` is also the order at every real point of `Y`. The real points
of `Y` are exactly the possible cluster centres. The Liouville step needs this order uniformly
in the arc position, which is why `ord_Y` is the relevant quantity.

## 2. The peeling lemma and the slice structure

**Lemma 2.1 (peeling).**

* *Hypotheses.* Let `A` in `Q^d` be finite. Let `L_1, .., L_r` be affine-linear functions with
  `A` contained in `Z(L_1) cup ... cup Z(L_r)`. Put
  `A_i = A cap Z(L_i) minus (Z(L_1) cup ... cup Z(L_{i-1}))`.
* *Conclusion.*
  ```text
  reg(A) <= max_{i : A_i nonempty} ( reg(A_i) + i - 1 ).
  ```

*Proof.* Induct on `r`. Let `A' = A_2 cup ... cup A_r = A minus Z(L_1)`; it is covered by
`L_2, .., L_r`. Let `g : A -> Q` be given.

1. Choose `P_1` of degree `<= reg(A_1)` with `P_1 = g` on `A_1`.
2. `L_1` does not vanish on `A'`. Choose `P'` of degree `<= reg(A')` with
   `P' = (g - P_1)/L_1` on `A'`.
3. Then `P = P_1 + L_1 P'` equals `g` on `A_1` and on `A'`, and
   `deg P <= max(reg A_1, reg A' + 1)`.

Now apply the induction hypothesis to `A'`. `[]`

**Lemma 2.2 (slices of the slice).**

* *Hypotheses.* `M >= 2` and `p, n >= 0`.
* *Conclusion.* `Q^M(p,n)` is the disjoint union, over `c = -n, .., p`, of the nonempty
  hyperplane sections `{k_1 = c}`. The section at `c` is
  ```text
  {c} x Q^{M-1}(p - c, n)      (0 <= c <= p),         {c} x Q^{M-1}(p, n - |c|)      (-n <= c <= -1).
  ```
  Its regularity index is that of the `(M-1)`-dimensional set.

*Proof.* `(c, k')` lies in `Q^M(p,n)` iff `sum k' = p - n - c` and `|k'|_1 <= p + n - |c|`.

* For `c >= 0` this reads `Q^{M-1}(p-c, n)`.
* For `c < 0` it reads `Q^{M-1}(p, n-|c|)`.
* If `c > p`, then `|k'|_1 >= |sum k'| = c + n - p`, which forces `|c| + |k'|_1 > p + n`. The
  case `c < -n` is symmetric.
* `Q^1(p', n') = {p' - n'}` is never empty. `[]`

**Lemma 2.3 (sorted peeling).** Let `r_1, .., r_nu` be the slice regularities, sorted so that
`r_(1) >= ... >= r_(nu)`. For real `v` put `N(v) = #{i : r_i >= v}`.

* Lemma 2.1, with the slices taken in sorted order, gives
  `reg <= max_i (r_(i) + i - 1)`.
* If `N(v) <= K - v + 1` for every `v` in `{r_1, .., r_N}`, then
  `max_i (r_(i) + i - 1) <= K`.

*Proof.* `N(r_(i)) >= i`. `[]`

## 3. The upper bound

For `M >= 2` and `1 <= j <= M-1` put

```text
l^M_j(p,n) = (M-1) * ( p/(M-j) + n/j ),        phi_M(p,n) = min_{1<=j<=M-1} l^M_j(p,n).
```

**Theorem 3.1.** For all `M >= 2`, all integers `p, n >= 0` and all `1 <= j <= M-1`,

```text
reg Q^M(p,n)  <=  l^M_j(p,n).
```

Hence `reg Q^M(p,n) <= phi_M(p,n)`.

*Proof.* Induct on `M`.

*Base case `M = 2`.* `Q^2(p,n)` consists of `p+n+1` points on a line, so
`reg = p + n = l^2_1(p,n)`.

*Inductive step, `M >= 3`.* Put `m = M - 1 >= 2`, and assume the theorem for `M - 1`, i.e. for
all indices `1 <= j' <= m-1`. Write the slice regularities of Lemma 2.2 as

* `r^+_c = reg Q^{M-1}(p-c, n)` for `0 <= c <= p`;
* `r^-_c = reg Q^{M-1}(p, n-c)` for `1 <= c <= n`.

`N(v)` denotes the number of slices with regularity `>= v`. By Lemma 2.3 it suffices to show
`N(v) <= K - v + 1` for every slice value `v`, where `K = l^M_j(p,n)`.

*Case `j = 1`.* Here `K = p + m n`.

* By the hypothesis with `j' = 1`: `r^+_c <= (p - c) + (m-1)n` and
  `r^-_c <= p + (m-1)(n-c) <= p + (m-1)n =: V`. So every slice value is `<= V`.
* For `v <= V`, at most `min(p+1, V - v + 1)` values `c` in `[0,p]` satisfy
  `r^+_c >= v`. At most `n` values `c` in `[1,n]` satisfy `r^-_c >= v`.
* If `v >= (m-1)n`: `N(v) <= V - v + 1 + n = K - v + 1`.
* If `v < (m-1)n`: `N(v) <= p + n + 1 < K - v + 1`.

*Case `j = m`.* The map `k -> -k` is an affine bijection `Q^M(p,n) -> Q^M(n,p)`. Also
`l^M_m(p,n) = m p + n = l^M_1(n,p)`. So this case is the case `j = 1` at `(n,p)`, which was
proved above in the same step.

*Case `2 <= j <= m-1`* (so `M >= 4`). Put

```text
alpha = (m-1)/(m-j),   beta = (m-1)/(j-1),   so 1/alpha + 1/beta = 1;
V+ = l^{M-1}_j(p,n),   V- = l^{M-1}_{j-1}(p,n).
```

The hypothesis with the two admissible indices `j' = j` and `j' = j-1` gives

```text
r^+_c <= l^{M-1}_j(p-c, n)     = V+ - c*alpha     (0 <= c <= p),
r^-_c <= l^{M-1}_{j-1}(p, n-c) = V- - c*beta      (1 <= c <= n).
```

Both `l`'s are non-decreasing in each argument. So every slice value is
`<= min(V+, V-)`.

The identity below is a direct computation, also checked in (A) of `certify.py`:

```text
V+/alpha + V-/beta = (p + n(m-j)/j) + (p(j-1)/(m-j+1) + n) = m p/(m+1-j) + m n/j = l^M_j(p,n) = K.   (3.1)
```

So `K` is a convex combination of `V+` and `V-`, and `K >= min(V+, V-)`.

Fix a slice value `v`, so `v <= min(V+, V-)`. Put `x+ = (V+ - v)/alpha >= 0` and
`x- = (V- - v)/beta >= 0`.

* The values `c` in `[0,p]` with `r^+_c >= v` all satisfy `c <= x+`, so there are at most
  `x+ + 1` of them.
* The values `c` in `[1,n]` with `r^-_c >= v` all satisfy `c <= x-`, so there are at most `x-`
  of them.
* By (3.1) and `1/alpha + 1/beta = 1`, `N(v) <= x+ + x- + 1 = K - v + 1`. `[]`

**Corollary 3.2.** For all `M >= 2` and `p + n >= 1`:

```text
reg S^M(p,n) <= reg Q^M(p,n) <= (2 - 1/ceil(M/2)) (p + n).
```

Hence `tau(M) <= tau*(M) <= 2 - 1/ceil(M/2) < 2`.

*Proof.*

* If `M = 2k` is even, then `l^M_k(p,n) = ((2k-1)/k)(p+n)`.
* If `M = 2k+1` is odd, then
  `(l^M_k + l^M_{k+1})(p,n)/2 = ((2k+1)/(k+1))(p+n)`, and the minimum is at most the mean.
* In both cases the factor is `2 - 1/ceil(M/2)`. `[]`

The crude version of the induction carries only the symmetric bound `tau'(p+n)` for `M-1`.
It gives `reg Q^M(p,n) <= tau'(p+n) + (2 - tau') min(p,n)`, hence the recursion
`tau_M <= 1 + tau_{M-1}/2`, i.e. `tau*(M) <= 2 - 2^{2-M}`. That is already a barrier, but
the exact value needs the whole family `l_j`. The middle facets `2 <= j <= M-2` come from
pairing the bound `l_j` on the positive side with `l_{j-1}` on the negative side, which is what
makes `1/alpha + 1/beta = 1`.

**Remark 3.3 (curves: BDPP (i) for the outer facets).** Use the slice model (§6.2): Laurent
polynomials on `Q^M(p,n)` are restrictions of bidegree-`(p,n)` forms to the closure `Gamma_M`
of `{(w, w^{-1})}` in `P^{M-1} x P^{M-1}`. Let `C` be a general line through `e = (1,..,1)` in
the first factor, `w = 1 + lambda v`.

* Its image in the second factor is `[prod_{i != j}(1 + lambda v_i)]_j`, of degree `M - 1`.
* A Laurent polynomial `f` on `Q^M(p,n)` restricts to `G(lambda)/prod_i(1 + lambda v_i)^n`,
  with `deg G <= p + (M-1)n`.
* Lines through `e` cover the torus, so `f|_C != 0` for a general `C`.
* Hence `ord_e f <= p + (M-1)n = l^M_1(p,n)`. Lines in the second factor give `l^M_{M-1}`.

So the outer facets have covering-curve proofs. For `2 <= j <= M-2`, BDPP duality predicts
movable classes with `(H_1.Gamma, H_2.Gamma)/E.Gamma = ((M-1)/(M-j), (M-1)/j)`. I did not
construct them. Lemma 2.1 replaces them.

On the smooth model `X^#` of `Bl_Y X_M`, the supremum of `E.Gamma/H.Gamma` over movable classes
`Gamma` is exactly `1/tau*(M) = k/(2k-1)`, with `k = ceil(M/2)` (this uses Theorem 5.1 below).

* It is `<= 1/tau*` because `H - tau* E` is `Q`-effective, and movable classes are
  nonnegative on effective divisors.
* It is `>= 1/t` for every `t > tau*` by BDPP, because `H - tE` is then not
  pseudo-effective.

Since `k/(2k-1) > 1/2`, movable classes with contact ratio `> 1/2` exist for every `M`. So the
counting heuristic of `tau_problem.md` ("no movable rational curves with contact ratio
`>= 1/2` for `M >= 5`") can at most concern families of rational curves, not movable classes.
The heuristic itself is not tested here.

## 4. The lower bound: alternating orbits, products, and the cone tiling

**Proposition 4.1 (alternating orbits).**

* *Hypotheses.* Let `a, b >= 0` with `r = a + b + 1 <= M`. Let
  `v = (-a, .., -1, 0, 1, .., b, 0, .., 0)` in `Z^M`: `r` distinct entries, then `M - r`
  zeros. Put
  `c = sum_{sigma in S_r} sgn(sigma) delta_{sigma v}`, with `S_r` permuting the first `r`
  coordinates.
* *Conclusion.* `c != 0`, `c` is supported on the shell `S^M(T_b, T_a)`, and
  `ord(c) >= C(r,2)`.

*Proof.*

* The entries of `v` in the first `r` places are distinct, so the `r!` points `sigma v` are
  distinct and `c != 0`.
* For a polynomial `P`, `sum_k c_k P(k) = (Alt P)(v)`, where
  `Alt P = sum_sigma sgn(sigma) P o sigma`.
* `Alt P` is alternating in the first `r` variables, so it is divisible by
  `prod_{i<j<=r}(k_i - k_j)`, which has degree `C(r,2)`.
* So `Alt P = 0` whenever `deg P < C(r,2)`. `[]`

Exact integer check (B) in `certify.py`: all moments vanish below `C(r,2)` and the pairing with
the Vandermonde is nonzero. This is done for all `(a,b)` with `M <= 5`, and for `r = M = 6, 7`
via power sums.

**Lemma 4.2 (products and monotonicity).**

* (i) If `f_i` is supported in `Q^M(p_i, n_i)` with `ord_1 f_i = m_i` (`i = 1, 2`), then
  `f_1 f_2 != 0` is supported in `Q^M(p_1+p_2, n_1+n_2)` and has order `m_1 + m_2`.
* (ii) If `p' <= p` and `n' <= n`, then
  `Q^M(p', n') + (p-p') e_1 - (n-n') e_2` is contained in `Q^M(p,n)`. So
  `reg Q^M(p,n)` is non-decreasing in `p` and in `n`.

*Proof.*

* (i) The lowest homogeneous parts (in `log w`) multiply, and `C[y]` is a domain. Supports
  add, and `|k + k'|_1 <= |k|_1 + |k'|_1`.
* (ii) Check the two defining conditions. Multiplying by a monomial does not change the order
  at `1`. `[]`

**Lemma 4.3 (cone tiling).**

* *Hypotheses.* For `1 <= j <= M` put `u_j = (T_{M-j}, T_{j-1})`. These are the costs of
  Proposition 4.1 with `r = M` and `(a,b) = (j-1, M-j)`.
* *Conclusion.*
  1. `l^M_j(u_j) = l^M_j(u_{j+1}) = C(M,2)` for `1 <= j <= M-1`.
  2. The slopes `T_{j-1}/T_{M-j}` strictly increase from `0` (at `j = 1`) to `infinity`
     (at `j = M`). So the cones `C_j = R_{>=0} u_j + R_{>=0} u_{j+1}` (`1 <= j <= M-1`)
     tile the quadrant.

*Proof.*

* `l^M_j(u_j) = (M-1)((M-j+1)/2 + (j-1)/2) = C(M,2)`.
* `l^M_j(u_{j+1}) = (M-1)((M-j-1)/2 + (j+1)/2) = C(M,2)`.
* The slopes are clear. `[]`

**Theorem 4.4 (exactness).**

* (a) For every `1 <= j <= M`:
  `reg S^M(u_j) = reg Q^M(u_j) = C(M,2) = phi_M(u_j)`.
* (b) For every integer pair `(p,n)`:
  ```text
  lim_{k -> infinity} reg Q^M(kp, kn) / k  =  phi_M(p,n).
  ```
  Moreover `phi_M = l^M_j` on the cone `C_j`.

*Proof.*

* (a) The lower bound `>= C(M,2)` comes from Proposition 4.1 with `r = M`, and monotonicity
  in `A`. For the upper bound, Theorem 3.1 with `l_j` (if `j <= M-1`) or `l_{j-1}` (if
  `j >= 2`) gives `<= C(M,2)`, by Lemma 4.3. Then
  `C(M,2) <= reg Q^M(u_j) <= phi_M(u_j) <= C(M,2)`, so `phi_M(u_j) = C(M,2)`.
* (b) Write `(p,n)` in `C_j` as `lambda u_j + lambda' u_{j+1}` with `lambda, lambda' >= 0`.
  Let `f_j` be the alternating Laurent polynomial of (a). Then
  `f_j^{floor(k lambda)} f_{j+1}^{floor(k lambda')}`, times a monomial, is supported in
  `Q^M(kp, kn)` by Lemma 4.2. Its order is
  `>= (k(lambda + lambda') - 2) C(M,2) = k l^M_j(p,n) - 2 C(M,2)`.
  So `liminf reg Q^M(kp,kn)/k >= l^M_j(p,n) >= phi_M(p,n)`. Theorem 3.1 gives
  `reg Q^M(kp,kn) <= k phi_M(p,n)`. `[]`

## 5. The value of `tau(M)` and exact formulas for small `M`

**Theorem 5.1 (main theorem).** For every `M >= 2`:

```text
tau(M) = tau*(M) = max_{p+n>=1} reg Q^M(p,n)/(p+n) = max_{p+n>=1} reg S^M(p,n)/(p+n) = 2 - 1/ceil(M/2).
```

The maximum is attained at the following points, and at no point does the ratio exceed it:

* `M = 2k+1`: `(p,n) = (T_k, T_k)`, `t = k(k+1)`, order `k(2k+1)`. These are the Vandermonde
  measures `v = (-k..k)`.
* `M = 2k`: `(p,n) = (T_k, T_{k-1})` or `(T_{k-1}, T_k)`, `t = k^2`, order `k(2k-1)`.

*Proof.* The upper bound is Corollary 3.2. Theorem 4.4(a) at `j = k+1` (for `M = 2k+1`) or
`j = k` (for `M = 2k`) gives the ratios `(2k+1)/(k+1)` and `(2k-1)/k`. `[]`

**Corollary 5.2 (barrier, and the Ru–Vojta constant).**

* *Hypotheses.* `M >= 2` and `k = ceil(M/2)`.
* *Conclusion.*
  1. Every form `F` of degree `D` not vanishing on `X_M` has
     `ord_Y F <= (2 - 1/k) D`.
  2. In the Liouville estimate of `tau_problem.md` one has
     `|F(z)| <= K C^m R^{D - m/2}` with `D - m/2 >= D/(2k) > 0`. So the estimate never forces
     `F(z) = 0` at exponent `1/2`. To do so it would need an extra saving of
     `R^{D/(2k)}` from elsewhere.
  3. On `Bl_Y X_M`, `H - tE` is big iff `t < 2 - 1/k`. For rational `t` it is `Q`-effective
     iff `t <= 2 - 1/k`.
  4. `beta(H,E) <= 2 - 1/k`, and also `b + beta(H - bE, E) <= 2 - 1/k` for `0 <= b < tau*`.
  5. On arcs of length `C R^alpha`, the single-place Liouville method on `X_M` needs
     `m(1 - alpha) > D`. By item 1 it works exactly for `alpha < 1 - 1/tau* = (k-1)/(2k-1)`;
     the Vandermonde `F` attains this. Proposition 6.1 of `conditional.md` (which needs
     `(1-alpha) beta > 1`) works at most in the same range. Both ranges are strictly below `1/2`
     for every `M`.

*Proof.*

* Item 1 is Lemma 1.2 with Theorem 5.1.
* Item 3: bigness for `t < tau*` holds because `H` is big and `H - tau* E` is `Q`-effective
  (Theorem 4.4(a)). For `t > tau*` the class is not even pseudo-effective, by item 1 and
  superadditivity.
* Item 4: in `beta(L,D) = liminf sum_{m>=1} h^0(NL - mD)/(N h^0(NL))` every summand is `<= 1`,
  and it vanishes for `m > N tau*`. The same holds for `L = H - bE`, with `m` ranging up to
  `N(tau* - b)`.
* Item 5: `1 - 1/(2 - 1/k) = (k-1)/(2k-1)`. `[]`

So the maximal Liouville exponent from `X_M` is `(k-1)/(2k-1)`: `1/3` (`M = 3, 4`), `2/5`
(`M = 5, 6`), `3/7`, and so on. As far as I recall, these are exactly the Cilleruelo–Córdoba
exponents `1/2 - 1/(4k-2)`, which come from a Vandermonde argument. I did not re-check the
original paper. Theorem 5.1 says that on these varieties nothing better than the alternant
exists.

**Theorem 5.3 (closed forms for `M = 3, 4`).** For all `p, n >= 0`:

```text
reg Q^3(p,n) = p + n + min(p,n),
reg Q^4(p,n) = floor( min( p + 3n, 3(p+n)/2, 3p + n ) ).
```

*Proof.* The upper bounds are Theorem 3.1. Note `phi_3 = min(2p+n, p+2n)` and
`phi_4 = min(p+3n, (3/2)(p+n), 3p+n)`. For the lower bounds, use Lemma 4.2 products of these
building blocks:

* the hexagon measure `(a,b) = (1,1)`: cost `(1,1)`, order 3;
* for `M = 4`, `(a,b) = (1,2)` and `(2,1)`: costs `(3,1)` and `(1,3)`, order 6;
* the linear factors `w_1 - w_2` (cost `(1,0)`) and `w_1^{-1} - w_2^{-1}` (cost `(0,1)`),
  order 1 each.

By symmetry take `p >= n`.

* `M = 3`: `n` hexagons and `p - n` linear factors give order `p + 2n`.
* `M = 4`, `p >= 3n`: `n` blocks of cost `(3,1)` and `p - 3n` linear factors give
  `p + 3n = phi_4`.
* `M = 4`, `n <= p < 3n`: write `p - n = 2x + eps` with `eps` in `{0,1}`. Take `x` blocks of
  cost `(3,1)`, `n - x` hexagons and `eps` linear factors. The cost is `(p,n)` and the order is
  `3n + 3x + eps = floor(3(p+n)/2)`.

Both identities were also checked by exhaustive arithmetic for `p, n <= 40` (check (E)). `[]`

For `M >= 5` the finite values can lie strictly below `floor(phi_M)`. Examples:
`reg Q^5(3,2) = 7 < 8 = phi_5(3,2)` and `reg Q^6(4,2) = 9 < 10 = phi_6(4,2)`. These are exact:

* the upper bound comes from a mod-`q` rank, which is rigorous in that direction;
* the lower bound comes from a product construction.

Only the asymptotic statement, Theorem 4.4(b), is exact for every `(p,n)`.

**Remark 5.4 (products of alternants are not always optimal; generalised alternants).** Let
`v_1, .., v_q` represent the `S_M`-orbits of points of `S^M(p,n)` with pairwise distinct
coordinates. Let `d >= 0`.

* *Claim.* If `q` exceeds the dimension of the space of symmetric polynomials of degree `<= d`
  restricted to `{sum k = p - n}`, then some nonzero `c = sum_i gamma_i Alt(delta_{v_i})` has
  `ord(c) >= C(M,2) + d + 1`.
* *Proof.* For `deg P <= C(M,2) + d`, write `Alt P = Vand * S` with `S` symmetric of degree
  `<= d`. Then `sum_k c_k P(k) = sum_i gamma_i Vand(v_i) S(v_i)`. That gives
  `dim Sym_{<=d}|_H` linear conditions on `gamma`, and `q` exceeds this number. `c != 0`
  because the orbits are disjoint. `[]`
* *Example.* For `M = 5` and `(p,n) = (6,3)` there are `q = 4` orbits, and `d = 3` gives the 3
  conditions `1, e_2, e_3`. So `ord >= 14 = phi_5(6,3)`. The products of Lemma 4.2 give only 13
  here.
* *Consequence.* `reg S^5(6,3) = reg Q^5(6,3) = 14 = phi_5(6,3)`. This is certified exactly by
  `gen_alternant.py`, with `gamma = (1, 1, -2, 3)` and a 480-point support on which all moments
  of degree `< 14` vanish. The upper bound is Theorem 3.1.

So the lattice of extremal measures is richer than products of alternants. This does not affect
any threshold.

## 6. Interpretations

### 6.1 The dual polytope

`phi_M` is the support function of the region
`Lambda_M = conv{ ((M-1)/(M-j), (M-1)/j) : 1 <= j <= M-1 } + R^2_{>=0}`.

* By Theorem 4.4, `Lambda_M` is exactly the set of `(alpha, beta)` with
  `alpha T_b + beta T_a >= C(a+b+1, 2)` for all `a + b <= M - 1`. These are the linear
  inequalities satisfied by all alternating-orbit generators.
* Proof of the last statement:
  - Each vertex `l_j` satisfies all these inequalities, by Theorem 3.1 and Proposition 4.1.
  - Conversely, if `(alpha, beta) . u_j >= C(M,2)` for all `j`, then by Lemma 4.3
    `alpha p + beta n >= l_j(p,n) = phi_M(p,n)` on each cone `C_j`. That is the support-function
    description of `Lambda_M`.
* `vertices.py` lists the vertices for `M <= 9`. The vertex `j` is tight on the generators
  `(a,b) = (j-1, M-j)` and `(j, M-1-j)`, which are `u_j` and `u_{j+1}`, possibly among others.
* So the whole asymptotic Hilbert-threshold structure is generated by the full alternants on
  all `M` coordinates.

### 6.2 Toric picture (task item (ii))

`Q^M(p,n)` is the set of lattice points of the Minkowski sum `p Delta + n(-Delta)`. Its Laurent
polynomials are the restrictions of bidegree-`(p,n)` forms to the closure `Gamma_M` of the graph
of the Cremona map `w -> w^{-1}` in `P^{M-1} x P^{M-1}`.

* For `M = 3`, `Gamma_3` is the hexagonal del Pezzo surface of degree 6, which is `P^2` blown
  up at the three coordinate points. Here `O(p,n) = (p+2n)H - n(E_1+E_2+E_3)`, and `e` is a
  fourth general point. The bound `phi_3 = min(p+2n, 2p+n)` is then the classical statement
  that, on `P^2` blown up at four general points, the conic class `2H - E_1 - .. - E_4` and
  the class `H - E_4` are nef.
* For general `M`, Theorem 4.4(b) computes the maximal vanishing order at a general point,
  that is, the pseudo-effective threshold at the identity of the torus, for every bundle
  `O(p,n)|Gamma_M` given by restricted forms: it is `phi_M(p,n)`.
* `X_M` packages these bundles over `s = p - n`. Its threshold is
  `max_{p+n=1} phi_M(p,n) = 2 - 1/ceil(M/2)`.

## 7. Numerics (consistency checks; all agree with the theorems)

`upper_checks/hilb.py` computes Hilbert functions of point sets by column-ordered elimination
modulo the two primes `1000003` and `999983`. The two primes always agreed.

* A mod-`q` regularity is a rigorous upper bound for the rational one.
* The constructive bound `L_M` (Lemma 4.2 products of the generators of Prop. 4.1) is a
  rigorous lower bound.
* The column `B` is the sorted-peeling recursion of §3 evaluated exactly. It is a rigorous
  upper bound, and `B <= phi_M` always.

Selected rows from `log_M3.txt` .. `log_M6.txt`. In every computed case
`L <= reg <= floor(phi)`.

| M | (p,n) | slice `Q` reg | shell `S` reg | `L` | `B` | `phi_M` |
|---|---|---|---|---|---|---|
| 3 | (1,1) | 3 | 3 | 3 | 3 | 3 |
| 3 | (5,5) | 15 | 8 | 15 | 15 | 15 |
| 4 | (3,1) | 6 | 6 | 6 | 6 | 6 |
| 4 | (4,4) | 12 | 9 | 12 | 12 | 12 |
| 4 | (5,4) | 13 | 10 | 13 | 13 | 13.5 |
| 5 | (3,2) | 7 | 7 | 7 | 8 | 8 |
| 5 | (3,3) | 10 | 10 | 10 | 10 | 10 |
| 5 | (6,1) | 10 | 10 | 10 | 10 | 10 |
| 5 | (6,3) | 14 | 14 | 13 (14 by Rem. 5.4) | 14 | 14 |
| 5 | (5,4) | 14 | 12 | 14 | 14 | 14.67 |
| 6 | (2,1) | 4 | 4 | 4 | 4 | 5 |
| 6 | (3,3) | 10 | 10 | 10 | 10 | 10 |
| 6 | (4,2) | 9 | 9 | 9 | 9 | 10 |
| 6 | (5,2) | 10 | 10 | 10 | 11 | 11.25 |
| 6 | (4,3) | 11 | 11 | 11 | 11 | 11.67 |

Logs, all in `upper_checks/`: `log_M3.txt` (`t <= 12`), `log_M4.txt` (`t <= 10`),
`log_M5.txt` (`t <= 9`) and `log_M6.txt` (`t <= 7`). Every row passes (`L <= reg_q <= floor(phi_M)`):

* M = 3: 96 rows;
* M = 4: 70 rows;
* M = 5: 58 rows;
* M = 6: 38 rows.

For the slices `Q`, the rigorous lower bound (`L`, or Remark 5.4 at `M = 5, (6,3)`) equals the
mod-`q` upper bound in every computed case. So all slice values in these logs are exact rational
values, not just evidence: 48/48 slices for `M = 3`, 35/35 for `M = 4`, 29/29 for `M = 5`, and
19/19 for `M = 6`. The shell values are mod-`q` upper bounds only.

Further consistency points:

* The lead's values of `tau` (`1, 1.5, 1.5, 1.667, 1.667, 1.75` for `M = 2..7`) are exactly
  `2 - 1/ceil(M/2)`.
* The referee's full-slice values (`referee_checks/log_tau.txt`: `M = 4, N = 8 -> 12`;
  `M = 5, N = 6 -> 10`) agree with Theorem 4.4.
* A concurrent independent audit, `scratchpad/tau/verify.md` (not by me), reports two things:
  - it independently found the shell-decoupling error, and proved `tau <= M/2`;
  - its full-slice data satisfy `reg(A_{s,D}) = floor(phi_M)` in 202 of 203 slices
    (`M = 2..5`). The exception is `(M,p,n) = (5,3,2)`: 7 < 8, as here.

Exact checks in `certify.py` (log `log_certify.txt`), all PASS:

* (A) The peeling recursion `B_M <= phi_M`, on 5128 cases with `M <= 14`. Also the identity
  (3.1) and the slope identities of the inductive step.
* (B) The alternating-orbit certificates.
* (C) The cone tiling, `phi_M(u_j) = C(M,2)` and `max ratio = 2 - 1/ceil(M/2)`, for `M <= 60`.
* (D) The `X_2` decoupling counterexample.
* (E) The closed forms for `M = 3, 4`.

## 8. What the barrier does not exclude (question (C) of the problem)

1. **Several shells, or `N`-dependent coefficients that are polynomial in `N`.** This is
   exactly `tau*`, so there is no gain.
2. **Coefficients depending on the circle or the arc position.**
   * Here `F` must be integral with controlled height.
   * A height `R^{eta D}` costs `eta D` in the exponent.
   * So the requirement becomes `m > (2 + 2 eta) D` (on the same normalisation), which is
     still impossible by Theorem 5.1.
   * A genuinely different input would be vanishing at a special point of `Y` rather than
     along `Y`. But clusters sit at arbitrary positions, and the torus acts transitively on
     `Y minus {uv=0}`, so this needs the arithmetic of the cluster.
3. **Forced Gaussian divisibility.** It would have to supply a factor of modulus
   `> R^{D/(2 ceil(M/2))}` at every cluster. `auxiliary.md` Lemma A5 (built-in divisibility
   never exceeds the archimedean size, and is zero for squarefree `N`) already rules out the
   natural version of this.
4. **More places or divisors in Ru–Vojta/Subspace** (e.g. `S` containing primes of `N`). These
   are not covered by Corollary 5.2. They are the `S`-unit route recorded as failing on
   uniformity in `claimed_complete_proof_salvage_audit.md`.

## 9. Files

All files are in `scratchpad/tau/upper_checks/`. The parent directory `scratchpad/tau/` is
shared with a concurrent agent (`verify.md`, `taulib.py`, `e*.py` are theirs). Mine are
`upper.md` and this subfolder. Run the scripts from inside `upper_checks/`.

* `certify.py`: exact checks (A)–(E), about 6 s; log `log_certify.txt`.
* `gen_alternant.py`: the exact generalised-alternant certificate of Remark 5.4; log
  `log_gen_alternant.txt`.
* `hilb.py`: modular Hilbert functions and regularity of `Q^M(p,n)` and `S^M(p,n)`.
* `peel.py`: the peeling recursion `B_M`, the semigroup lower bound `L_M`, and the LP value.
* `data_check.py`: produces `log_M3.txt` .. `log_M6.txt` (`python3 data_check.py M tmax ball,shell`).
* `vertices.py`: the vertices of `Lambda_M` for `M <= 9`; log `log_vertices.txt`.
* `check_B_phi.py`: `B_M <= phi_M` for `M <= 10`, and `max B/t = 2 - 1/ceil(M/2)`; log `log_check_B_phi.txt`.
