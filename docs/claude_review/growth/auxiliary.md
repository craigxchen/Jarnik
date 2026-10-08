# Auxiliary polynomials, gap principles, determinant methods, larger sieves and the polynomial method at exponent 1/2

Route: methods beyond pairwise Vandermonde and beyond the project's character certificates.
Tools examined: (a) Thue–Siegel/Dyson auxiliary polynomials, (b) gap principles of Thue-counting type,
(c) Heath-Brown's p-adic and Salberger's global determinant methods, (d) larger-sieve variants,
(e) the polynomial (Stepanov) method.

Exact checkers (all pass): `aux_gauss.py`, `aux_checks.py` (output in `aux_checks_output.txt`),
`aux_paley_design.py` (output in `aux_paley_design_output.txt`), `aux_saturation.py`, all in this directory.

## 0. Verdict

**No new growth theorem.** All five tool families were pushed until they gave an explicit inequality.
Each inequality falls into one of three cases:

1. It is provably dominated by the pairwise valuation data `v_p(z_i-z_j)`, that is, by
   Vandermonde/Cilleruelo–Córdoba plus collision counting. This covers (a), (c) and (e) (Section 2)
   and the Liouville form of (b) (Section 3).
2. It is a sum of pair inequalities. The best such sum is the uniformly weighted one
   (Theorem D1, a symmetrization argument), and an explicit abstract model satisfies all
   per-pair data at `W = O(M log M)` (Proposition D3). This covers every larger-sieve variant (d)
   and Salberger's method.
3. It gains only on special arithmetic structure that generic extremal profiles need not have:
   * the Baker-improved perfect-power certificates (Proposition B2), which need a non-saturated
     difference lattice;
   * the Legendre-twisted class statistics (Proposition D2), which need a biased Walsh coefficient.

The comparison bound is the current general one: `(1+eps) log R/loglog R` (notes), or
`(2/3+eps)` in the lead's chain-lemma sketch. No method here beats it. The best general inequality
any of them produces is the sieve inequality `W >= 4 M log M - O_C(M)` (when no prime `<= M` divides
`N`). Determinant/Thue–Siegel/Stepanov/gap methods without the global sieve give only
`M = O_C(log R)`.

The following are new rigorous statements. None of them improves the point count.

| item | statement |
|---|---|
| Prop. A | Every polynomial interpolation determinant on the circle equals `2^{-mD} prod z_j^{-D} V(z) Phi(z)`, with `V` Vandermonde and `Phi` a symmetric Gaussian integer whose value at a cluster is its confluent Wronskian up to `O(delta)`. Factorwise determinant arguments are therefore implied by the Vandermonde inequality on the same points. |
| Lemma A5 | Built-in Gaussian divisibility of auxiliary values never exceeds their archimedean size. It is zero for squarefree `N`. |
| Prop. B1 / Lemma B | The Thue-type gap principle for anchor-relative heights is equivalent to the pair inequality. Heights above `N^{3/4}C` are isolated (at most one such point). |
| Prop. B2 | Perfect-power certificates: `W/4 < 2^40 V log(e(4D+1)) + log(C‖y‖_1)`. Liouville gives `D V/2` in place of the `2^40 V log` term. |
| Thm. D1 | Uniform pair weights are optimal among all weightings of pair inequalities under cut and class-count constraints. |
| Prop. D2 | Legendre-twisted sieve: `W >= (4/M) sum_l log l * max(0,(M^2+S_l^2)/nu_l - M) - 4(M-1) log C`, where `S_l` is a Walsh coefficient of the exponent code at the Legendre set `{k : (p_k/l) = -1}`. |
| Prop. D3 | An abstract Paley prime box with an explicit refining residue design satisfies every per-pair larger-sieve inequality, at all prime-power moduli and with arc order, at `W ≈ 20 M log M`. |
| Lemma D4 | Torsion regime: if all `rho_k` have order `<= T` modulo every `l ∈ 𝓛`, then `M <= 4T(1 + W/(4 theta(𝓛)))(1+o(1))`. But every such prime satisfies `log p_k >= (4 theta(𝓛) - 2T log 4)/(T(T+1))`, so [BS] already gives the same order. |

## 1. Setting and baseline (notation of the notes)

* **Normalization.** Divide by the global Gaussian gcd. This keeps the hypothesis with the same `C`,
  since `L/|g| <= C |g|^{-1/2} (N/|g|^2)^{1/4}`.
* **Points.** `N = prod_k p_k^{e_k}` with `p_k ≡ 1 (mod 4)`, and
  `z_j = eps_j prod_k pi_k^{a_jk} conj(pi_k)^{e_k - a_jk}` for `j = 1..M`.
  All points lie on an arc of length `L <= C N^{1/4}`.
* **Notation.** `W = log N`, `g_ij = gcd(z_i, z_j)`, `d_ij = sum_k |a_ik - a_jk| log p_k`,
  `u_i = z_i/g_ij`, `gamma_ij = log Norm(u_i - u_j) >= 0`.

**[P] Pair inequality.**

    log Norm(z_i - z_j) = W - d_ij + gamma_ij <= 2 log L <= W/2 + 2 log C.

**[Cut] Cut bound.**

    sum_{i<j} d_ij = sum_{k,tau} log p_k * S_{k,tau} (M - S_{k,tau}) <= (M^2/4) W.

**[Sigma] Summed inequality.** Summing [P] and using [Cut]:

    sum_{i<j} gamma_ij + sum_{k,tau} log p_k (S_{k,tau} - M/2)^2 <= W M/4 + 2 log C * binom(M,2).

**[Sieve] Collision lower bound.** Let `l^a` be a prime power coprime to `2N`. Then
`z_i ≡ z_j (mod l^a)` if and only if `l^a | u_i - u_j`. Each such collision adds `2 log l` to
`gamma_ij`. The residues of the points modulo `l^a` lie in `nu(l^a) = (l ∓ 1) l^{a-1}` classes
(split / inert). Pigeonhole and Mertens in progressions then give:

* `W >= 4 M log M - O_C(M)` if no prime `<= M` divides `N`, so `M <= (1/2 + o(1)) log R/loglog R`;
* with inert primes only, `M <= (1 + eps) log R/loglog R` (inert_prime_cofactor_bound.md);
* `W >= 3 M log M - O_C(M)` in general, by the lead's chain lemma (sketch, not re-audited here).

**[BS] Bounded support** (bounded_split_support_reduction.md, re-derived). For `C <= 1`:

    M <= 4 ceil(W/beta) + 1,      beta = max_k e_k log p_k.

The proof halves the family at the largest block and applies the Cilleruelo–Córdoba (Ramana)
determinant to `2s+1` points with `s = ceil(W/beta)`.

**Remark 1.1 (the hard regime).** Suppose `C <= 1` and `M >= eps W/log W`. Then [BS] forces every
prime-power block to satisfy `e_k log p_k < 4W/(M-5) = (4/eps + o(1)) log W`. So `N` is a product of
at least `(eps/4 - o(1)) W/log W` blocks, each of size at most `(log N)^{O(1/eps)}`, and
`M ≍ omega(N)`. Any growth improvement must handle these `(log N)^{O(1)}`-smooth `N`, with
`≍ log N/loglog N` prime factors mostly in `[M, M^{O(1)}]`.

## 2. Interpolation methods on a rational curve: (a), (c), (e)

### 2.1 Factorization (new, exact)

**Proposition A.** Let `f_1, ..., f_m ∈ Z[i][X,Y]` have total degree `<= D`. Put

    G_l(Z) := 2^D Z^D f_l( (Z^2+N)/(2Z), (Z^2-N)/(2iZ) ) = sum_{k=0}^{2D} c_{lk} Z^k ∈ Z[i][Z].

For points `z_1, ..., z_m` of `x^2 + y^2 = N`:

1. `2^{mD} prod_j z_j^D * det[f_l(z_j)] = det[G_l(z_j)]`.
2. `det[G_l(z_j)] = V(z) Phi(z)`, where `V = prod_{i<j}(z_j - z_i)` and
   `Phi = ± sum_K det(c_{l,k})_{k∈K} s_{lambda(K)}(z_1, ..., z_m)` is symmetric with `Z[i]` coefficients.
   Here `K = {k_1 < ... < k_m} ⊂ {0..2D}` and `lambda(K)_t = k_{m+1-t} - (m-t)`.
3. `Phi(w, ..., w) = Wr(w) := det[ G_l^{(t)}(w)/t! ]_{l, 0<=t<m}`.
4. If `z_j = w(1 + h_j)` with `|h_j| <= eta`, then `|Phi(z) - Wr(w)| <= ((1+eta)^Lambda - 1) Phi^+(|w|)`.
   Here `Lambda = max |lambda(K)|` and `Phi^+(rho) = sum_K |det c_K| s_{lambda(K)}(1^m) rho^{|lambda(K)|}`.

*Proof.*
1. On the circle `x = (z + N/z)/2` and `y = (z - N/z)/(2i)`. Since `a + b <= D`, the substitution of a
   monomial `X^a Y^b` multiplied by `2^D Z^D` is `2^{D-a-b} Z^{D-a-b} (Z^2+N)^a (-i)^b (Z^2-N)^b`.
2. Cauchy–Binet expands `det[sum_k c_lk z_j^k]`, and the bialternant formula gives
   `det[z_j^{k_t}] = ± V s_lambda` with integral Schur polynomials.
3. This is the confluent limit of divided differences.
4. `s_lambda` is a sum of `s_lambda(1^m)` monomials of degree `|lambda|`. ∎

**Corollary A (what the determinant knows).** Let `p` be a Gaussian prime, `p ∤ 2`. Then

    v_p(det f) = v_p(V) + v_p(Phi) - D sum_j v_p(z_j),
    v_p(V)     = sum_{i<j} [ v_p(g_ij) + v_p(u_i - u_j) ].

So the conductor and collision content of `V` is exactly the pair data, and
`log Norm(det f) + 2mD log 2 + mDW = log Norm(V) + log Norm(Phi)`.

Bombieri–Pila, Heath-Brown and Salberger arguments are **factorwise**:

* the archimedean bound is a Taylor expansion, i.e. a product of chords times derivative sup-norms;
* the arithmetic bound is `p`-adic valuations of differences at points of one residue class.

Such an argument shows `G_V G_Phi <= B_V^2 B_Phi^2`, where `G_Phi` is a valid lower bound and
`B_Phi` a valid upper bound for the nonzero Gaussian integer `Phi`. Then
`G_Phi <= Norm(Phi) = |Phi|^2 <= B_Phi^2`. Hence a determinant contradiction (`G > B^2`) forces
`G_V > B_V^2`, which is the Vandermonde (Cilleruelo–Córdoba) contradiction on the same points with
the same pair data.

The degenerate case `Wr(w) ≈ 0` gives no arithmetic gain either. Item 4 then makes `|Phi| < 1`, so
`Phi = 0` and `det = 0`: some combination of the `G_l`, of degree `<= 2D`, vanishes at all `m` points,
and Bezout only gives `m <= 2D`.

The same holds with derivative rows (multiplicities): confluent Vandermonde times a symmetric factor.
It also holds for tensor-product bases on `(circle)^k` (Kronecker structure gives powers of
one-variable Vandermondes).

*Exact check (A in `aux_checks.py`).* 24 random determinants on actual best clusters of 4 circles
(`m = 3, 4, 5`, random `f_l` with `Z[i]` coefficients):

* identity 1 holds;
* `V | det G` in `Z[i]`;
* `|Phi - Wr(z_1)|/|Wr(z_1)| <= 0.0114`.

### 2.2 Best inequalities produced

| variant | inequality | growth |
|---|---|---|
| pure Vandermonde on subsets (Cilleruelo–Córdoba mechanism) | `C < 1`: `M <= 1 + W/(4 log(1/C))`. Covering by `ceil(2C)` arcs of length `N^{1/4}/2`: `M <= ceil(2C) (1 + W/(4 log 2))` | `O_C(log R)` |
| Heath-Brown, one auxiliary prime `l ∤ 2N`, `l > C`, points in one class | the class has `m >= M/nu_l` points and every pair gains `gamma >= 2 log l`, so [Sigma] gives `(m-1) log(l/C) <= W/4`; hence `M <= nu_l (1 + W/(4 log(l/C)))` | `O_C(log R)` |
| Heath-Brown at a conductor prime (bad reduction: two lines) | a residue class is a level class, i.e. descent; gives [BS], and with the sieve on the sub-cluster `W >= 2 M log M + (M/2) log p - O_C(M)` for `p ‖ N` (Remark 5.2) | linear in `W` |
| Salberger's global method | the class statistics `sum_P binom(n_P, 2)` at all `l` are exactly [Sieve] | `W >= 4M log M - O_C(M)`, i.e. `(1/2) log R/loglog R`; in general `(2/3)` (lead) or `(1+eps)` (notes) |

For a fixed number `k+1` of points, all of these reduce to the Cilleruelo–Córdoba exponent
`1/2 - 1/(4 floor(k/2) + 2)` (plus constants). Nothing new for fixed `k`.

### 2.3 Thue–Siegel/Dyson specifics (a)

**Lemma A5 (built-in divisibility never beats size).** Let `F = sum_{t=0}^D c_t Z^t conj(Z)^{D-t}` and
`z = eps prod_k pi_k^{a_k} conj(pi_k)^{e_k - a_k}`. The gcd of all monomial values `z^t conj(z)^{D-t}` is
`prod_k p_k^{D min(a_k, e_k - a_k)}`. Its norm is at most `N^D = |z^t conj(z)^{D-t}|^2`, with equality
only if every `a_k = e_k/2`. For squarefree `N` it is `1`.

*Proof.* The `pi_k`-exponent of `z^t conj(z)^{D-t}` is `t a_k + (D-t)(e_k - a_k)`, and the
`conj(pi_k)`-exponent is `t(e_k - a_k) + (D-t) a_k`. The minimum over `t ∈ [0, D]` of each is
`D min(a_k, e_k - a_k)`. ∎

So an auxiliary polynomial's value at a cluster point never carries forced conductor content beyond
its own archimedean size. Any gain must come from vanishing at other points, which is the
determinant/Vandermonde mechanism of Section 2.1.

**Structural reason (Dirichlet quality).** Roth/Dyson-type methods, and Thue's gap principle, control
approximations to one target that are better than Dirichlet's exponent, with heights in geometric
progression.

* Relative to an anchor, the targets are rational (the direction `1`), and every point has the same
  height `N`.
* The anchor ratios `P_i/conj(P_i)` (Section 3) have heights `Q ∈ [N^{1/2}/C^2, N]` and quality
  `delta ≈ C N^{-1/4} ∈ [Q^{-1/2}, Q^{-1/4}]`.
* Dirichlet's exponent on the circle is `≈ 1/Q`, since Gaussian integers of norm `<= Q` give `≈ Q`
  directions.

The cluster points are *worse* than generic approximations. Counting them is a divisor-structure
problem, not a Diophantine-approximation problem. Liouville is attained at rational targets: the pair
`X ± i` on `N = X^2 + 1` has normalized constant `2/N^{1/4} → 0`.

### 2.4 Stepanov (e)

Stepanov's method needs the set, modulo a prime `l`, to lie in a sparse algebraic set such as a coset
of a small subgroup, where high-order vanishing can be imposed cheaply. Modulo `l ∤ 2N` the cluster
lies in `z_* <rho_k>` with `z_* = prod conj(pi_k)^{e_k}` and `rho_k = pi_k/conj(pi_k)`, inside the
cyclic group `G_l` of order `nu_l`. This is all of `G_l` unless every `p_k` satisfies a power-residue
condition at `l`. In that case the gain is exactly the class-statistics gain of Section 4.3, a factor
`r` at the prime `l` only. Lemma D4 shows that uniform congruence gains force large primes, which
descent already exploits. Over `Z`, the polynomial method is the determinant method of Section 2.1.

## 3. Gap principles (b)

**Proposition B1 (anchor-relative heights; exact).** For an anchor `z_0` put `P_i := z_i/gcd(z_i, z_0)`.
Then:

* `z_i/z_0 = (unit) P_i/conj(P_i)`;
* `P_i` is conjugate-primitive, and `Q_i := Norm(P_i) = e^{d_0i}`;
* for `i ≠ j`, `conj(P_i) P_j = G_ij Gamma_ij` with `G_ij ∈ Z_{>0}`, `Gamma_ij` conjugate-primitive,
  `G_ij^2 e^{d_ij} = Q_i Q_j` and `Norm(Gamma_ij) = e^{d_ij}`;
* `Im(conj(P_i) P_j) = G_ij t_ij` with `t_ij = Im Gamma_ij ≠ 0`.

The Thue gap argument for two approximations `P_i/conj(P_i)` and `P_j/conj(P_j)` of the same rational
target runs as follows.

* Since `z_i/z_j = (unit) conj(Gamma_ij)^2/Norm(Gamma_ij)` and the chord is at most `L`,
  `|Gamma_ij - u conj(Gamma_ij)| <= |Gamma_ij| * L/R` for a unit `u`.
* The left side is a nonzero Gaussian integer, so it is at least 1. It is nonzero because `Gamma_ij`
  is conjugate-primitive and nonunit (`d_ij > 0`).
* Hence `e^{d_ij/2} >= R/L >= N^{1/4}/C`.

Because `G_ij^2 e^{d_ij} = Q_i Q_j`, this "gap" statement `G_ij <= (L/R) sqrt(Q_i Q_j)` is *literally*
the pair inequality [P] with `gamma` dropped:

    G_ij <= C N^{-1/4} sqrt(Q_i Q_j)   <=>   d_ij >= W/2 - 2 log C.

*Proof of the identities.* Write `u_i = z_i/g`, `u_0 = z_0/g`. Then `u_0 = (unit) conj(u_i)`: at every
prime the relative exponents of `z_i` and `z_0` are opposite. The exponent bookkeeping gives the
rational content `G_ij = prod_k p_k^{m_k}`, where `m_k = min(|a_ik - a_0k|, |a_jk - a_0k|)` when the two
relative exponents have the same sign and `m_k = 0` otherwise. The primitive part has exponent
`|a_jk - a_ik|`. ∎

Checked exactly on 100 anchor pairs of actual clusters (`aux_checks.py`, B).

**Lemma B (height isolation; squarefree `N`).** At most one cluster point has
`d_0i > 3W/4 + log C`.

*Proof.* With `abar_0 = 1 - a_0`, `d_ij <= d(a_i, abar_0) + d(abar_0, a_j) = 2W - d_0i - d_0j`.
Combined with `d_ij >= W/2 - 2 log C`, this gives `d_0i + d_0j <= 3W/2 + 2 log C`. ∎

This is the only "gap" present. The anchor heights of a balanced cluster all lie near `W/2`
(Hadamard profiles), so there is no geometric growth of heights and a Thue-style count gives nothing.

**Block form ("simultaneous approximation by block products").** Write `z = A·B` with `A` the part on
a set `T` of thresholds. If `W_T > W/2 + 2 log C`, then [P] makes `z ↦ (levels on T)` injective on the
cluster. This gives `M <= prod_{T}(e_k + 1)`, which is useless when blocks are small.

**Proposition B2 (Baker-improved perfect-power certificates).** Let `y ∈ Z^M` with `sum_i y_i = 0`, and
suppose `sum_i y_i a_i = D beta` with `D >= 1` and `beta ∈ Z^s \ {0}`. Put
`Gamma = prod_k pi_k^{beta_k}` (negative powers read as powers of `conj(pi_k)`) and
`V = log Norm(Gamma) = sum_k |beta_k| log p_k`. Then:

1. *(exact)* `Pi_y := prod_{y_i>0} z_i^{y_i} prod_{y_i<0} conj(z_i)^{|y_i|} = n·eps·Gamma^{2D}` with
   `n ∈ Z_{>0}` and `eps` a unit. Checked on 511 certificates (`aux_checks.py`, C2).
2. *(Liouville / project character inequality)* `D V/2 >= W/4 - log(C ‖y‖_1/2)`.
3. *(Matveev)* `W/4 < 2^40 · V · log(e(4D+1)) + log(C ‖y‖_1)`.

*Proof.*
1. The `pi_k`-exponent minus the `conj(pi_k)`-exponent of `Pi_y` is `2 sum_i y_i a_ik`.
2. `sum y_i = 0` cancels the arc centre, so `|arg Pi_y| <= ‖y‖_1 delta/2` with `delta = L/R <= C e^{-W/4}`.
   Hence `|4D arg Gamma - k pi| <= ‖y‖_1 delta` for some integer `k` with `|k| <= 4D + 1`.
   * Liouville: `|Gamma^{2D} - i^k Norm(Gamma)^D| = |Gamma|^D |Gamma^D - i^k conj(Gamma)^D| >= Norm(Gamma)^{D/2}`.
     This uses `Gamma`, `conj(Gamma)` coprime and nonunit.
   * Matveev: use the degree-two specialization of Matveev 2000, Cor. 2.3, as audited in the notes
     (`log |u log alpha - v log(-1)| > -2^40 V log(e max(|u|,|v|))` for `alpha = Gamma/conj(Gamma)`,
     `Gamma` conjugate-primitive; items 256 and 276), with `u = 2D`, `v = k`. ∎

*Use.* Item 3 beats item 2 exactly when `D ≫ 2^41 log D`.

* Example: a normalized Hadamard core of order `s`, with `y` a non-constant column. Then
  `D = s/2`, `beta = e_k`, `‖y‖_1 = s`, so `W/4 < 2^40 log p_k log(e(2s+1)) + log(Cs)`. Choose `k`
  with `p_k = p_min` among the non-constant columns. Since `W >= (s-1) log p_min`, this is impossible
  once `s - 1 > 2^42 log(e(2s+1)) + 4 log^+(Cs)` (this recovers item 276).
* *Obstruction.* If the difference lattice `Lambda_B = {sum y_i (a_i - a_0)}` is saturated in `Z^s`,
  then `D beta ∈ Lambda_B` implies `beta ∈ Lambda_B`. Item 2 for that certificate then forces
  `V >= W/2 - O(log(C ‖y'‖_1))`, and item 3 is useless.
* In clusters with `M - 1 < s` varying primes, `e_k` is generically not even in the rational span.

*Exact data (finite evidence, not proof).*

* On 120 best clusters (`m = 3..8`) of 5 circles, the largest invariant factor of `B` is 1 in 116
  cases and 3 in 4 cases (`aux_saturation.py`).
* No coordinate certificate exists in 484 tested cases (`aux_checks.py`, C: rank deficit).

So the general Baker loss is the missing small exponent of `Z^s/Lambda_B`. Cramer bounds the exponent
by `|det| <= s^{s/2} H^s`, which costs `log D ≈ (s/2) log s` against a budget `W/(2^40 log p) ≈ s/2^40`.
This matches the obstruction recorded under item 276.

## 4. Larger-sieve variants (d) and Salberger class statistics

### 4.1 Prime powers and composite moduli

Composite moduli are redundant by CRT: the information is the set of prime powers dividing
`u_i - u_j`. Prime powers `l^a`, `a >= 2`, add at most

    sum_l sum_{a>=2} 2 log l * M^2/(2 nu(l) l^{a-1}) = M^2 sum_l log l/(nu(l)(l-1)) = O(M^2),

that is, they change only the `O_C(M)` term in `W >= 4 M log M - O_C(M)`.

### 4.2 Weights cannot help (new)

**Theorem D1 (symmetrization).** Fix any weights `rho_ij = rho_ji >= 0` (not all zero) and sum [P]
with them. Suppose an argument uses only:

* (i) the cut structure, `sum rho_ij d_ij <= W·maxcut(rho)`;
* (ii) for each modulus in any family of prime powers coprime to `2N`, a lower bound for
  `sum_{i ~ j} rho_ij` valid for **every** refining system of partitions with the given class counts.

Then the resulting bound `W >= W_rho` satisfies `W_rho <= W_unif`, the bound from uniform weights.
The same holds for conductor level structures and for any permutation-invariant constraint family.

*Proof.* Let `S*` be a maximal cut and `P*` a collision-minimizing partition system for uniform
weights, and let `sigma` be a uniformly random relabelling of the points. Then:

* `maxcut(rho) >= E cut_rho(sigma S*) = |rho| maxcut_u/binom(M,2)`;
* `min_P sum_{i ~_P j} rho_ij <= E sum_{i ~_{sigma P*} j} rho_ij = |rho| coll_u/binom(M,2)`,
  level by level.

With `W_rho = [sum_l 2 log l · coll_l(rho) - 2 log C |rho|]/(maxcut(rho) - |rho|/2)`, the numerator
scales down and the denominator scales up by the same factor `|rho|/binom(M,2)` relative to uniform
weights. ∎

Both symmetrization inequalities are checked exactly by brute force on 299 random integer
weightings (`aux_checks.py`, F2). The eigenvalue form `maxcut(rho) - |rho|/2 <= (M/4) mu(rho)`,
`|rho| <= binom(M,2) mu(rho)` (`mu = -lambda_min`) is checked on 2000 random `rho` (F).

Configuration-dependent weights amount to the full per-pair system. Its abstract feasibility at
`W = O(M log M)` is Proposition D3.

### 4.3 Class statistics (new, exact)

Fix an odd prime `l ∤ N`. Let `G_l = {x ∈ (Z[i]/l)^* : x conj(x) = 1}` (cyclic, order
`nu_l = l - (-1/l)`), `z_* = prod_k conj(pi_k)^{e_k}`, `x_i = z_i/z_* ∈ G_l` and `rho_k = pi_k/conj(pi_k)`.
Then `x_i = eps_i prod_k rho_k^{a_ik}`, and

    #{ordered colliding pairs incl. diagonal} = sum_g n_g^2 = (1/nu_l) sum_{chi ∈ hat(G_l)} |sum_i chi(x_i)|^2.

**Proposition D2.** For the quadratic character `chi_2` of `G_l`:

    chi_2(x_i) = (2/l)^{u_i} prod_k (p_k/l)^{a_ik}      (eps_i = i^{u_i}),

because `chi_2(rho_k) = (p_k/l)` and `chi_2(i) = (2/l)`. Consequently

    coll_l >= (M^2 + S_l^2)/(2 nu_l) - M/2,      S_l := sum_i (2/l)^{u_i} prod_k (p_k/l)^{a_ik},

and [Sigma] gives the **Legendre-twisted sieve**

    W >= (4/M) sum_{l odd, l ∤ N} log l * max(0, (M^2 + S_l^2)/nu_l - M) - 4(M-1) log C.

*Proof of the character values.*

* Inert `l`: `rho = pi^{1-l}`, and `rho` is a square in `G_l` iff `pi ∈ F_l^* (F_{l^2}^*)^2` iff
  `Norm(pi) = p` is a square mod `l`.
* Split `l`: `rho ≡ pi^2/p (mod l)`.
* `i` is a square in `G_l` iff `l ≡ ±1 (mod 8)`, i.e. iff `(2/l) = 1`. ∎

Checked:

* the formula on 6180 (point, `l`) pairs over 5 circles, including units and exponents `> 1`;
* the collision bound on 1560 actual point sets (`aux_checks.py`, D and D2).

`S_l` is the Walsh coefficient of the exponent code at the "Legendre set" `{k : (p_k/l) = -1}`
(by reciprocity `(p_k/l) = (l/p_k)`). In the worst case `S_l = 0` for all `l <= M`; a code balanced on
these `≈ M/log M` parity functions is consistent with everything above. The term therefore improves
nothing in general. If every `p_k` is a QR modulo every `l <= M`, it only doubles the main term
(`W >= 8 M log M`), a constant factor.

### 4.4 Congruence gains force large primes

**Lemma D4 (torsion regime).** Let `𝓛` be a set of odd primes `l ∤ N`, let
`theta(𝓛) = sum_{l ∈ 𝓛} log l`, and let `T >= 1`. Suppose that for every `l ∈ 𝓛` the subgroup of `G_l`
generated by the `rho_k` (`k` varying) has order `<= T`. Then:

1. *(sieve)* For `M'` points in one unit class (`M' >= M/4`),
   `theta(𝓛)(M'/T - 1) <= W/4 + (M'-1) log C`. So `M <= 4T(1 + W/(4 theta(𝓛)))(1 + o(1))` when
   `theta(𝓛)/T ≫ log C`.
2. *(size)* Every varying prime satisfies
   `(T(T+1)/2) log p_k + T log 4 >= 2 theta(𝓛)`.

*Proof.*

1. `G_l` is cyclic, so the `x_i` of one unit class lie in one coset of a subgroup of order `<= T`.
   There are `>= (M'/2)(M'/T - 1)` collisions, each worth `2 log l`; apply [Sigma] to the sub-cluster.
2. If `rho_k` has order `t_l <= T` at `l`, then `l | pi_k^{t_l} - conj(pi_k)^{t_l} =: X_t ≠ 0`.
   Distinct `l` with the same `t` give `prod l^2 | Norm(X_t) <= 4 p_k^t`. Sum over `t <= T`. ∎

So genuine super-pigeonhole collisions from class statistics occur only when the primes are huge:
`log p_k ≳ 4 theta(𝓛)/T^2`. There [BS] gives `M <= 4W/log p_max + 1 ≲ T^2 W/theta(𝓛)`, the same
order up to a factor `T`. The special case `T = 1` with `q` a Gaussian ideal (`rho_k ≡ unit (mod q)`,
hence `Norm q <= 4 p_k`) is the simplest instance. Index-`r` conditions with `r` bounded (e.g. all
`(p_k/l) = 1`) multiply the collisions at `l` by `r` only, a constant factor. No such condition is
forced in the hard regime of Remark 1.1.

### 4.5 Abstract feasibility of all per-pair data (rigorous obstruction)

**Proposition D3.** Fix `C > 0`. For every prime `q ≡ 3 (mod 4)` with `M = q+1` large, there exist:

1. an exponent code: the rows of the normalized Paley–Hadamard matrix restricted to its `q`
   non-constant columns, each column placed on `5` distinct actual split primes. These are
   `r = 5q` split primes `p_j` with `log p_j ∈ [tau, tau + 1/r)` and `tau ∈ [4 log M, 4 log M + log 2]`,
   so `W = 20 q log M (1 + o(1))`. They exist by the notes' pigeonhole: partition `[M^4, 2M^4]` into
   `r` equal intervals and use the fixed-modulus PNT, so no short-interval prime theorem is needed;
2. a refining residue design. For each prime power `l^a` it uses at most `nu(l^a)` classes, giving
   point `x ∈ {1..M}` the class:
   * `x mod l^a` for `l = 2` or `l ≡ 3 (mod 4)`;
   * `x mod l^a` if `l ∤ x`, and `(x - 1/2) mod l^a` if `l | x`, for `l ≡ 1 (mod 4)`.
   Conductor primes `> M^4` get injective classes.

These satisfy, for every pair, `demand_xy <= excess_xy`, where

    demand_xy := sum_{colliding levels} 2 log l <= 2 log n + 2 log(2n+1) + 2 log(2n-1) <= 6 log(2M+1),
    excess_xy := d_xy - W/2 + 2 log C >= (5/2) 4 log M - o(1) + 2 log C,      n = |x - y|.

They also satisfy the *ordered-arc* version: points at normalized positions `xC/M` with
`gamma_xy := excess_xy - 2 log C + 2 log(C n/M) ∈ [demand_xy, excess_xy]`.

*Proof.*

* Paley–Hadamard orthogonality gives `sum_j s_xj s_yj = -5` over the `5q` physical columns. With
  `log p_j = tau + eps_j`, `0 <= eps_j < 1/r`, we get
  `d_xy - W/2 = -(1/2) sum_j (tau + eps_j) s_xj s_yj >= (5 tau - 1)/2 >= 10 log M - 1/2`.
* Deleting rows preserves every pair inequality: columns that become constant join the gcd, which
  lowers `W` and raises every excess. So every large `M` (with `q <= 2M`) is covered.
* Collisions of the design at `l^a` occur exactly when `l^a` divides `n`, or (split `l`, mixed type)
  `2n ± 1`. Hence `demand <= 2 log n + 2 log(2n+1) + 2 log(2n-1)`.
* The class maps land in the residues prime to `l` (split case) and refine across `a`.
* The ordered version needs `10 log M + 2 log(Cn/M) >= 6 log(2n+1)`, true for `M` large. ∎

*Exact check* (`aux_paley_design.py`, `C = 1`, `q = 11, 19, 23, 31, 43`). It uses the `5q`
consecutive split primes above `M^4`, whose logarithms differ by less than `2·10^{-3}`:

* all class counts and refinements hold;
* min over pairs of `excess - demand` is `17.5 .. 24.1`;
* min of `gamma_line - demand` is `16.7 .. 21.4`;
* `W/(M log M) ≈ 18.4 .. 19.6`.

*Meaning.* Every inequality produced by (c) (any auxiliary primes, Heath-Brown or Salberger) and by
(d) (any moduli, prime powers, weights, cut-profile dependence) uses only the per-pair collision data,
the cut structure, the class-count constraints and the arc order. All of these hold simultaneously in
this abstract model at `W ≈ 20 M log M`. So none of these methods can prove `W ≥ ω(M log M)`, i.e.
improve the growth `log R/loglog R`.

The model is not a lattice-point configuration. It violates the homomorphic structure of the
residues: actual residues are `x ↦ eps ∏ rho_k^{a_k}`, linked to the exponents. That structure is the
class statistics of Section 4.3. The exact Ptolemy relations are not in the model either (checked
separately: `n_A a + n_B b = n_C c` on 39 actual quadruples, `aux_checks.py`, E).

## 5. Natural restricted regimes

**5.1 Bounded `omega(N)`.** Already uniform by [BS] (`M <= 4K+1`, notes). Nothing new.

**5.2 Dyadic or large-prime regimes (descent + sieve; combination of known tools).** Let `p ‖ N` and
assume no prime `<= M` divides `N`. Take the larger level class of `p` (`m >= M/2` points). It is a
cluster on `N/p` with constant `<= C p^{-1/4}`, and the sieve applies inside it. So

    (M/2)(log(M/2) + (1/4) log p - log C - O(1)) <= W/4,
    W >= 2 M log M + (M/2) log p - O_C(M).

If all prime factors of `N` lie in `[P, 2P]` with `P = M^A`, then `M <= (2/(4+A) + o(1)) W/log M`. This
improves the sieve's `1/4` iff `A > 4`, and two-step descent is no better. It is `Θ(W/log W)` unless
`A → ∞`, where [BS] already gives `o(W/log W)`. No growth change.

**5.3 Fixed `k`.** Every method here reduces to Cilleruelo–Córdoba's `k`-point exponent plus constants.

## 6. Missing lemma and why the methods stall

**Missing Lemma (super-pigeonhole reduced-chord energy).** There is `psi(M)` with
`psi(M)/log M → ∞` such that for each `C` and every primitive cluster of `M >= M_0(C)` points on an arc
of length `<= C N^{1/4}`:

    sum_{i<j} log Norm( (z_i - z_j)/gcd(z_i, z_j) ) >= psi(M) M^2.

By [Sigma] this gives `W >= 4 M psi(M) - 4M log C`, hence `M = o(log R/loglog R)`. Growth
`psi(M) >= c M^kappa` would give `M = O_C((log R)^{1/(1+kappa)})`. The uniform theorem needs the
`W`-dependent form `sum gamma > WM/4 + 2 log C binom(M,2)` for `M >= M(C)`. That is the handoff's
`H_sum >= c w` (item 532) in another normalization.

**Why each method stalls.**

* **(a), (c), (e):** On the rational curve `x^2 + y^2 = N` every interpolation determinant is
  `V × Phi`. `Phi` has no forced content beyond its size (Prop. A, Lemma A5). The method sees only
  `v_p(z_i - z_j)`, and their sum is Plotkin-critical.
* **(b):** Anchoring makes the target rational, and Liouville is attained. Baker gains only on
  perfect-power certificates, which need a non-saturated difference lattice. Generic and actual small
  clusters have saturated lattices (data).
* **(d):** All per-pair information is feasible at `W = O(M log M)` (Prop. D3), and no reweighting
  helps (Thm D1). The only non-pairwise remnant inside (c)/(d) is the homomorphic class structure. Its
  computable consequence, the Legendre-twisted term (Prop. D2), vanishes on Legendre-balanced codes.

The missing input must tie the finite class maps (discrete logarithms of `rho_k = pi_k/conj(pi_k)`
modulo `l`) to the archimedean arguments `arg pi_k`. Within (a)–(e) the only such tie is an exact
multiplicative relation, i.e. a character certificate, which is outside this route and handled
elsewhere in the project.

## 7. Exact computations performed

| check | content | size |
|---|---|---|
| A | Prop. A: identity 1, `V \| det G` in `Z[i]`, Wronskian proximity | 24 determinants on actual clusters |
| B | Prop. B1 identities (conjugate-primitivity, `G^2 e^{d_ij} = Q_i Q_j`, `G \| Im`, `t ≠ 0`) | 100 anchor pairs |
| C, C2 | Prop. B2(1): `Pi_y = n·unit·pi_k^{±2D}` (`D` up to 11); no coordinate certificate in 484 cluster cases | 511 certificates |
| saturation | largest invariant factor of difference lattices of best clusters: 1 in 116 of 120 | 120 clusters |
| D, D2 | Prop. D2 character formula; Legendre-twisted collision bound | 6180 point–prime pairs; 1560 point sets |
| E | Ptolemy: `X + Y = Z` collinear with 0, `n_A \| x`, `n_B \| y`, `n_C \| z` | 39 quadruples |
| F, F2 | Thm D1: eigenvalue form; exact brute-force symmetrization inequalities (maxcut and partition minima, `nu = 2, 3`) | 2000 random matrices; 299 exact weightings |
| G | Lemma D4(2) torsion divisibility; Lemma B isolation | 768 `(p, l)` pairs; 225 anchors |
| design | Prop. D3 at `q = 11..43` (pair and ordered-arc margins) | 5 values of `q` |

Finite checks support the stated identities. They are not used as proofs of any asymptotic claim.
