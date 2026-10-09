# Outside the barred classes? GAP/Fourier, energy, incidences, equidistribution: reductions and a fake-circle barrier

Round 4, target: an argument outside the barred classes of `round4_context.md` along four directions:
(a) the global GAP/torus-orbit structure of `E_R` (Fourier/Riesz products, Bohr sets, inverse
Littlewood-Offord); (b) sum-product and multiplicative energy; (c) incidence geometry; (d)
equidistribution of `E_R` along subsequences. Checks are in `round4/outside_checks/`. Every
script runs from that directory with `python3 <script>`, and its output is stored next to it as
`*_out.txt`.

Conventions. *Theorem*, *Proposition* and *Lemma* mean proved here with all quantifiers, or
cited with the exact source. *Exact check* means integer or `Fraction` arithmetic. *Illustration*
means floating point; it is never used as a proof. `W = log(R^2) = log N`. An arc of length
`<= C sqrt(R)` has angular width `Delta <= C N^(-1/4) = C e^(-W/4)`.

---------------------------------------------------------------------------------------------

## 0. Verdict

**No growth improvement and no uniform theorem.** Each of the four directions reduces to a
barred class. For (a) and (b) the reduction is now a theorem, not a heuristic. The main new
results are:

1. **Theorem A.1 (completion).** Assume all varying primes are large in the explicit sense
   `prod_j (1 + 2/(sqrt(p_j) - 1)) < 1.75`, for example `p_j >= 21 r^2`. Then the Liouville
   ("axis/diagonal gap") inequalities for **all** monomials `prod pi_j^(v_j)`, `v in Z^r`, carry
   exactly the same information about a cluster as the rational-character inequalities of
   `growth/walsh-extraction.md`, imposed simultaneously. Monomials outside the rational span of the
   cluster's exponent differences impose nothing: a Haar-positive proportion of all real angle
   vectors realising the cluster satisfies every one of their gaps.
2. **Theorem A.3 (fake circles, rigorous barrier).** Take any `C > 0` and any `M >= M_0(C)` that
   is the order of a Hadamard matrix. Then there is an explicit **fake circle**: genuine primes
   `p_j = 1 mod 4`, a genuine exponent profile, and real angles. It has an `M`-point cluster on an
   arc of length `C sqrt(R)` with `log R <= (15 + o(1)) M^2 log M`. It satisfies:
   * every monomial gap, for every nonzero integer vector;
   * every character gap strengthened by **any** residue-collision information modulo prime powers
     `<= 2M`, the range in which pigeonhole forces collisions;
   * every exact exponent identity (gcd norms, cut identity, levels).

   **Corollary A.4.** No argument that uses the lattice points only through these inputs can prove
   `M = o((log R / loglog R)^(1/2))`. Such arguments may also use any property of the angles that
   holds off a Haar-null set. The class includes:
   * every Fourier/Riesz-product, Bohr-set, inverse Littlewood-Offord and Delsarte-type argument on
     the angle GAP;
   * every multiplicative-energy argument;
   * all pair, character, residue-pigeonhole, cut and chain certificates of the notes, and any
     combination of them.

   A variant (Theorem A.3') adds residue strengthening for *all* monomials. The ceiling then becomes
   `(log R)^(1/3)`.
3. Directions (c) and (d) are orthogonal or local:
   * **(c), Prop. C.1.** The triangle-area (Jarnik) information is an exact product of pair
     identities. As a formula,
     `4(2 Area)^2 = Norm(gcd(z_a,z_b,z_c))^2 nu_ab^2 nu_bc^2 nu_ac^2`.
   * **(c), Prop. C.2.** Point-line incidences give a uniform bound only for arcs near a
     low-height lattice direction `v`: `M <= C(K + C)|v| + 2` when the arc is within `K R^(-1/2)`
     of `v`. In general they give `2 C^(3/2) R^(1/4) + 2`.
   * **(d), Prop. D.1.** Clusters live in a density-zero set: `#{N <= X : M_N >= 2} << C X^(3/4) log X`.
   * **(d), Prop. D.2.** A cluster carries at most `M 2^(-M-1)` of the angular mass. So
     equidistribution and weak-limit statements are blind to it.

| direction | strongest rigorous inequality obtained here | compared with the known `(2/3+eps) log R/loglog R` | class |
|---|---|---|---|
| (a) GAP / Bohr / Fourier / inverse LO | `M <= ceil(C/sqrt2) (1 + rho(W))`, where `rho(W) log(4 rho/e) = W`, so `M <= (2 ceil(C/sqrt2)+o(1)) log R/loglog R` (Prop A.4). Inverse-LO hypotheses force `r = O(A log A)` (Prop A.5) | weaker by the factor `3 ceil(C/sqrt 2)`; same growth | = rational characters (Thm A.1); capped at `sqrt(log R/loglog R)` (Thm A.3) |
| (b) energy / sum-product | the cluster is additively and multiplicatively Sidon; energies are 4-term characters; these are subcritical on average by `W/4` (Prop B.2) | no bound | inside the character class |
| (c) incidences | triangle identity = product of pair identities (C.1); lattice-line bound (C.2) | uniform only on a thin class of arcs | pair/local, or auxiliary polynomials |
| (d) equidistribution | `#{N<=X: M_N>=2} <= 16 C X^(3/4)(1+log(C X^(1/4))) + 14 C^2 X^(1/2)` (D.1); cluster mass `<= M 2^(-M-1)` (D.2) | says nothing for a single `N` | orthogonal to the worst case |

What lies outside every class covered by Theorem A.3 is listed in Section 6. It consists of:

* exact algebraic identities among the points (barred: `threshold/upper.md`, Ptolemy);
* the algebraicity of the prime arguments beyond integrality (Baker; useless here by `round4/flipped.md`, Prop B);
* residue information for *all* monomials at *all* moduli. This holds with equality at the actual
  angles, and in effect it is the actual primes themselves (§2.4);
* multi-place Diophantine approximation (barred).

---------------------------------------------------------------------------------------------

## 1. Setting, the information in `E_R`, and the gap lemma

Normalise as in `two_thirds/proof.md`, Lemma 1.1. Remove the Gaussian gcd of the cluster. This does
not increase `W` or the arc constant, and leaves `N` odd and free of primes `3 mod 4`. Then

```text
z_x = eps_x prod_{j=1}^r pi_j^(a_xj) conj(pi_j)^(e_j - a_xj),   x = 1..M,  0 <= a_xj <= e_j,
N = prod p_j^(e_j),  W = sum e_j log p_j,  eps_x in {1, i, -1, -i},
```

where `pi_j` is a fixed Gaussian prime above the split prime `p_j`. Fix real lifts `phi_j` of
`arg pi_j`. The full set `E_R` consists of the points `R i^u exp(i <2a - e, phi>)` with `a` in the
box `B = prod [0, e_j]` and `u in Z/4`. Its angle set is the generalised arithmetic progression with
generators `2 phi_j` and lengths `e_j + 1`. For `v in Z^r` put

```text
A_v = prod_j pi_j^(v_j^+) conj(pi_j)^(v_j^-),   Norm A_v = e^(w(v)),  w(v) = sum |v_j| log p_j,
arg A_v = <v, phi>  (mod 2 pi).
```

**Lemma 0.1 (gaps; walsh-extraction Lemma 1.1).** For every `v in Z^r \ {0}`:

```text
dist(<v,phi>, (pi/2) Z) >= arcsin(e^(-w(v)/2)),    dist(<v,phi>, pi/4 + (pi/2) Z) >= arcsin(e^(-w(v)/2)/sqrt 2).
```

In particular `dist(<v,phi>, (pi/4) Z) >= arcsin(e^(-w(v)/2)/sqrt 2)`.

*Proof.* Write `A_v = X + iY`. It is a conjugate-primitive nonunit of odd norm.
* If `Y = 0` or `X = 0`, then `A_v` is associate to its conjugate. That is impossible, because each
  `p_j` contributes only one of `pi_j` and `conj pi_j`.
* `X + Y` is odd, since `X^2 + Y^2` is odd. Hence `|X| >= 1`, `|Y| >= 1` and `|X +- Y| >= 1`.
* The sine of the distance to the axes is `min(|X|,|Y|)/|A_v|`. The sine of the distance to the
  diagonals is `min|X +- Y|/(sqrt 2 |A_v|)`. QED

*Realisation.* The cluster lies in an arc `theta_0 + [0, Delta]`. So there are `delta_x in [0, Delta]`
and integers `k_x` (absorbing units and windings) with

```text
<2 a_x - e, phi> = theta_0 + delta_x + (pi/2) k_x     (x = 1..M).                          (1.1)
```

**Definitions (angle models).** Fix distinct primes `p_j = 1 mod 4`, exponents `e`, and rows
`a_x in B`.

* A *profile* is `P = (p, e, (a_x), (delta_x), (k_x))`.
* `phi in R^r` *realises* `P` if (1.1) holds for some `theta_0`.
* `phi` is *gap-admissible* if the two inequalities of Lemma 0.1 hold for every `v != 0`.
* An *angle model* is a profile together with a gap-admissible realising `phi`. Its points are
  `R i^u exp(i <2a-e, phi>)`, `a in B`, with `R = N^(1/2)`.

By Lemma 0.1 every actual normalised cluster is an angle model. Methods (a) and (b) see only this
data: the GAP of angles, the cluster, and the integrality of monomials. Possibly they also see
residues, as in Section 2.4.

**Characters.** Let `L = span_R{a_x - a_y} cap Z^r`. Each `v in L` is `sum lambda_x a_x` with
`sum lambda = 0` and `lambda` rational. By (1.1),

```text
<v, phi> = (1/2) sum_x lambda_x (delta_x + (pi/2) k_x) =: chi_P(v)                        (1.2)
```

for every realising `phi`, because `theta_0` and `<e, phi>` cancel. These are exactly the
integral/rational characters of `growth/walsh-extraction.md` (1.2). Call `P`
*character-admissible* if `chi_P(v)` satisfies the two gap inequalities with `w = w(v)` for
**every** `v in L \ {0}` simultaneously.

---------------------------------------------------------------------------------------------

## 2. Direction (a): GAP, Bohr sets, Fourier/Riesz products, inverse Littlewood-Offord

### 2.1 The completion theorem

**Theorem A.1.** Let `P` be a profile that has a realisation. Let `A~` be the `M x r` matrix with
rows `2a_x - e`, and let `U = {u : A~ u in R 1}`.

(i) *Determinacy.* The realisations form the coset `phi_0 + U`. For `v in Z^r`, the value
`<v, phi>` is the same for all realisations iff `v in L`, and then it equals `chi_P(v)`.

(ii) *Completion.* Suppose `Pi := prod_j (1 + 2/(sqrt(p_j) - 1)) < 1.75` and `P` is
character-admissible. Then the gap-admissible realisations have Haar proportion at least
`1 - 1.33 (Pi - 1) > 0` in the realisation torus `(phi_0 + U)/2 pi Z^r`.

*Proof.* (i) Let `phi, phi'` realise `P` with `theta_0, theta_0'`. Then
`A~(phi - phi') = (theta_0 - theta_0') 1`, so `phi - phi' in U`. The converse is the same
computation read backwards. `<v, .>` is constant on `phi_0 + U` iff `v` is orthogonal to `U`.

For any linear map `T` and subspace `L'`, `(T^(-1) L')^perp = T^*(L'^perp)`, because
`Tu in L'` iff `u` is orthogonal to `T^* mu` for every `mu in L'^perp`. Hence

```text
U^perp = A~^T (1^perp) = { sum_x lambda_x (2a_x - e) : sum lambda = 0 } = span_R{a_x - a_y}.
```

The value on `L` is (1.2).

(ii) `U` is cut out by integer equations, so `U cap Z^r` spans `U`. Hence `T_U := U/(U cap 2 pi Z^r)`
is a compact torus, and the realisations modulo `2 pi Z^r` form the coset `phi_0 + T_U`; give it
Haar measure.

* For `v in L \ {0}`, the gap inequalities hold by character-admissibility and (i).
* For `v` not in `L`, `u -> <v, u> mod 2 pi` is a nontrivial continuous homomorphism `T_U -> R/2 pi Z`.
  It is nontrivial because `v` is not orthogonal to the rational subspace `U`. Its image is a
  connected nontrivial closed subgroup, hence the whole circle. So `<v, phi> mod 2 pi` is
  Haar-uniform on the coset.
* The union of the two bad sets of Lemma 0.1 lies in `{dist(<v,phi>, (pi/4)Z) < arcsin(e^(-w(v)/2))}`.
  This set has probability `(8/pi) arcsin(e^(-w(v)/2)) <= 2.64 e^(-w(v)/2)`. Here
  `arcsin y <= 1.0367 y` for `y <= 5^(-1/2)`, and `w(v) >= log 5`.
* `v` and `-v` give the same event. Summing over `v != 0` and using
  `sum_(v in Z^r) e^(-w(v)/2) = prod_j sum_(k in Z) p_j^(-|k|/2) = Pi` bounds the bad proportion by
  `1.32 (Pi - 1)`. QED

The condition `p_j >= 21 r^2` for all `j` implies `1.33(Pi - 1) < 1`
(`check_completion.py`, Part 3).

**Corollary A.2 (reduction of (a) and (b)).** Suppose an argument proves "`M <= f(W, C)`" for every
angle model, possibly using in addition any property of `phi` that fails only on a Haar-null subset
of each realisation torus. Then `M <= f(W, C)` holds for every character-admissible profile with
`Pi < 1.75`.

*Proof.* Theorem A.1(ii) gives a Haar-positive set of gap-admissible realisations, and it meets the
conull set. QED

So, in the large-prime regime, every Fourier-analytic, Riesz-product, Bohr-set, Delsarte/SDP or
inverse-Littlewood-Offord argument on the angle GAP lies in the rational-character class of the
notes. This holds whether the argument is magnitude-only or phase-sensitive, and whatever
generic Diophantine or equidistribution property of the angles it invokes. Such arguments see
exactly the cluster's angle data and the integrality of all monomials. The notes' `Theorem F`
(residue pigeonhole on winding numbers) is one of these consequences: it is Lemma 0.1 applied to
`v in L`.

### 2.2 The fake-circle barrier

The notes' Prop P (`growth/walsh-extraction.md` §7.1) checks that fully flipped Paley profiles pass
each character inequality **individually**. Those are necessary conditions for a realisation.
Theorem A.1 shows that the right object is **simultaneous** admissibility. The next theorem
constructs simultaneously admissible profiles, and therefore genuine angle models with large
clusters. Its rate is worse than Prop P's `W = Theta(M log M)`, but nothing in it is heuristic.

**Construction.** Let `H` be a Hadamard matrix of order `M`, normalised so its first column is all
ones. Labels are `a = 1..M-1`, and `b = 6M`.

* *Columns.* Each label `a` gets `b` columns, partitioned as `J_a = F_a cup S_a cup O_a` with
  `|F_a| = |S_a| = 2` and `|O_a| = b - 4`. There are `r = b(M-1)` columns in all.
* *Flips.* Each column in `F = cup F_a` is *flipped* at one row `x(j)`. Every row receives one or two
  flips; this is possible since `M <= 2(M-1) <= 2M`. Each row `x` gets one *designated* flip `f(x)`,
  injectively.
* *Siblings.* Inside each label, a bijection `F_a -> S_a` assigns siblings. Write `s(x)` for the
  sibling of `f(x)`.
* *Sign matrix.* `S_xj = H(x, a(j)) (-1)^[j in F, x = x(j)]`. Rows are `a_x = (1 + S_x)/2` in
  `{0,1}^r`, `e = 1`, and all `k_x = 0`.
* *Primes and balancing.* The `p_j` are distinct primes `= 1 mod 4` in `[P, 2P]`, and
  `tau := log P`. The sets `O_a` are chosen so that their log-sums `Omega'_a` all lie in a window
  `[Omega'_0, Omega'_0 + 1/M]` (Lemma A.6).

**Lemma A.6 (balancing, elementary).** Let `Q` be a set of `n` distinct primes in `[P, 2P]`, and
suppose `n > k^2 m (k log 2/eta + 1)`. Then `Q` contains `m` pairwise disjoint `k`-subsets whose
log-sums lie in a common interval of length `eta`.

*Proof.*
1. The `k`-subset log-sums lie in an interval of length `k log 2`. Cut it into
   `J <= k log 2/eta + 1` intervals of length `eta`. One of them contains at least
   `binom(n,k)/J` subsets.
2. Pick subsets greedily, discarding after each pick the at most `k binom(n-1,k-1)` subsets that
   meet it.
3. The greedy step succeeds `m` times because `binom(n,k)/J > (m-1) k binom(n-1,k-1)`, that is,
   `n > k^2 (m-1) J`. QED

Take `k = 6M - 4` and `eta = 1/M`. The lemma then needs about `150 M^5` primes in `[P, 2P]`, plus
`4(M-1)` more for `F` and `S`. By the prime number theorem for the progression `1 mod 4`, we may
take `P ~ 400 M^5 log M`. Then:

* `tau = 5 log M + loglog M + O(1)`;
* `W = sum log p_j <= r(tau + log 2) = (30 + o(1)) M^2 log M`;
* `p_j >= 21 r^2`, so Theorem A.1(ii) applies.

**Lemma A.7 (character margins).** Assume `tau >= 10`. Let `c in Z^M`, `c != 0`, `sum c = 0`,
`n = ||c||_1`. Put `T_a = sum_x c_x H(x,a)`, `F(c) = sum_(a >= 1) |T_a|`,
`v(c) = (1/2) c^T S in Z^r` and `V_c = w(v(c)) = (1/2) sum_j w_j |c^T S_j|`. Then:

```text
(i)   F(c) >= M;
(ii)  L = { v(c) : c in Z^M, sum c = 0 }   (saturation: every rational character is integral);
(iii) V_c - W/2 >= tau (n + 0.86 M).
```

*Proof.*

(i) `T_0 = sum c = 0` and `sum_(a=0)^(M-1) T_a^2 = M ||c||_2^2` (Hadamard). With `|T_a| <= n`
this gives `F >= M ||c||_2^2 / n >= M`, since `c_x^2 >= |c_x|`.

(ii) Columns `f(x)` and `s(x)` differ only in row `x`, by `-2H(x,a)`. Hence for
`v = sum lambda_x a_x` we get `v_f(x) - v_s(x) = -lambda_x H(x,a)`. So integrality of `v` forces
`lambda in Z^M`, and then `v = (1/2) c^T S`. Integrality of `v(c)` holds because `c^T S_j` is
even.

(iii) Use the designated twins and the columns of `O_a`, and drop the remaining nonnegative terms.
Sibling and `O` columns are unflipped, so `c^T S_s(x) = T_a` and `c^T S_f(x) = T_a - 2c_x H(x,a)`.
Since `|T - 2c_x h| + |T| >= 2|c_x|` and all weights are `>= tau`,

```text
V_c >= tau n + (1/2) sum_a Omega'_a |T_a|,
V_c - W/2 >= tau n + (1/2) sum_a Omega'_a (|T_a| - 1) - (1/2) sum_a (Omega_a - Omega'_a).
```

Every `T_a` is even, so `|T_a| - 1` is `>= 1` or equals `-1`. The window condition then gives

```text
sum_a Omega'_a (|T_a| - 1) >= Omega'_0 (F - M + 1) - (M-1)/M.
```

Next, `Omega'_0 >= (b-4) tau`. Also `Omega_a - Omega'_a <= 4(tau + log 2)`, since that difference is
the weight of the four columns in `F_a cup S_a`. Finally `F >= M`. Together,

```text
V_c - W/2 >= tau n + (3M - 2) tau - 1/2 - 2(M-1)(tau + log 2) >= tau (n + M) - 1.4 M >= tau(n + 0.86M)
```

for `tau >= 10`. QED

The equal-weight forms of (i)-(iii) are checked exactly in `check_fake_barrier.py`. This covers
`F >= M`, the two bounds behind (iii), the rank/saturation statement (ii), and the twin structure.
The profiles tested are Sylvester `M = 8, 16, 32` and Paley `M = 12, 20, 24`, with
`b = 6` and `b = 6M`, on about 389,600 characters for each value of `b`: all characters with `||c||_1 <= 4`, all columns,
all two-column half-sums, and random larger ones. All pass.

**Theorem A.3 (fake circles).** Take `C > 0`, and let `Q_M = lcm(1, ..., 2M) <= e^(2.08 M)`. There
is `M_0(C)` such that for every `M >= M_0(C)` that is the order of a Hadamard matrix (for example
`M = 2^m`, or `M = q + 1` with `q = 3 mod 4` prime) the construction above yields an angle model
with the following properties:

1. it has `M` distinct points on an arc of length `<= C sqrt(R)`, where `log R = W/2 <= (15 + o(1)) M^2 log M`;
2. every `v in Z^r \ {0}` satisfies `dist(<v, phi>, (pi/4) Z) >= arcsin(e^(-w(v)/2))`;
3. every character `v in L \ {0}` satisfies the residue-strengthened gap
   `dist(<v, phi>, (pi/4) Z) >= arcsin(Q_M e^(-w(v)/2))`.

In particular `M >= (7.5^(-1/2) - o(1)) (log R / loglog R)^(1/2)`.

*Proof.*

1. *Choice of `delta`.* Draw `delta` uniformly from `[0, Delta]^M`, `Delta = C e^(-W/4)`. For a
   character `c`, the value `X_c = chi(v(c)) = (1/2) c^T delta` has density at most
   `2/(Delta ||c||_inf)` and range of length `Delta n/2`.
2. *One character.* With `eta_c = arcsin(Q_M e^(-V_c/2))`, which is at most `1.04 Q_M e^(-V_c/2)`
   by Lemma A.7(iii),

   ```text
   P(dist(X_c, (pi/4)Z) < eta_c) <= (2 Delta n/pi + 2) * 4 eta_c / Delta
                                   <= 1.04 Q_M e^(-(V_c - W/2)/2) (2.55 n e^(-W/4) + 8/C)
                                   <= 1.04 Q_M P^(-(n + 0.86M)/2) (2.55 n + 8/C).
   ```

3. *Union bound.* There are at most `3^M 2^(n-1)` vectors `c` with `||c||_1 = n`. Summing over
   `c != 0` gives at most

   ```text
   1.04 (2.55 + 8/C) (3 e^(2.08) P^(-0.43))^M sum_(n>=2) n (2/sqrt P)^n / 2.
   ```

   This is `< 1` for `M >= M_0(C)`. Note `P^(0.43) >= 13 M^(2.15)`; the sum over `n` is
   `O(1/P)`; and `M_0(C) = O(1 + log^+(1/C))`.

   So some `delta` makes the profile character-admissible in the strengthened form 3. This uses
   Lemma A.7(ii): `L` consists exactly of the `v(c)`.
4. *Completion and distinctness.*
   * `S` has rank `M`, because the twin differences span every `e_x`. So the profile is realisable
     for every `delta`.
   * Theorem A.1(ii) supplies a gap-admissible realisation. For `v` outside `L` it gives 2 directly,
     since the union bound there was run for the `(pi/4) Z` form.
   * The `M` points are distinct by the pair gaps (`c = e_x - e_y`).
5. *Size.* The bound on `log R` was computed above. With `loglog R = (2 + o(1)) log M` it gives the
   last statement. QED

`check_fake_instance.py` is an illustration of a concrete instance. With Sylvester `M = 8`, 336
consecutive primes `= 1 mod 4` above `10^9` and `C = 1/2`, all 3,868 tested characters meet
property 3 with log-margin at least 196. The plain gaps hold at 200,000 random monomials, with
minimum ratio 45.

**Corollary A.4 (barrier).** Consider an argument that proves a cluster bound `M <= f(R)` for all
large `R`, using about the configuration only:

* (I) its exponent profile and the exact identities between exponent data: gcd norms
  `Norm g_xy = N/P_xy`, the cut identity, levels and chains;
* (II) the gap inequalities of Lemma 0.1 for arbitrary monomials;
* (III) residue-collision strengthenings, modulo prime powers `<= 2M`, of the gap inequalities for
  character monomials (for any residue assignment);
* (IV) any property of the angles valid off a Haar-null subset of each realisation torus.

Then `f(R) >= (7.5^(-1/2) - o(1)) (log R/loglog R)^(1/2)` for infinitely many `R`. In particular no
such argument proves the uniform theorem.

*Proof.* The fake circles of Theorem A.3 satisfy (I)-(III), and (IV) by Theorem A.1(ii). QED

Classes (I)-(III) contain every input of `two_thirds/proof.md`:
* pair chords are the pair gaps, strengthened by the collision product of each pair, and that
  product divides `Q_M` because the collisions used there are modulo `q^a` with
  `(q - 1) q^(a-1) < M`;
* the cut and chain identities are exponent identities;
* the class contains pair, integral and rational character certificates and Theorem F;
* it contains the aggregate collision and capacity inequalities of Prop P.

So the combined local + character + GAP-Fourier toolkit has a proved ceiling between
`c (log R/loglog R)^(1/2)` (this theorem) and `(2/3) log R/loglog R` (`two_thirds`).

*Consistency checks against refereed theorems.*
* Theorem C of `walsh-extraction.md` (flipped Walsh, a character argument) forces
  `W/(M log M) -> infinity` at the slow item-335 rate. Our Sylvester fakes have
  `W/(M log M) ~ 30 M`.
* Theorems B and B'' need pure or almost pure cores; our fakes have every row flipped.
* `two_thirds` requires `W >= 3 M log M - O(M)`.
* `linear_allocation_affine_rigidity.md` requires `M <= r + 1 = 6M(M-1) + 1`.

**Theorem A.3' (all monomials strengthened).** Run the same construction with
`P >= 35 r^2 Q_M^2`, so `tau = (4.16 + o(1)) M`. Then also every `v` outside `L` satisfies
`dist(<v,phi>, (pi/4)Z) >= arcsin(Q_M e^(-w(v)/2))`, with `log R <= (12.5 + o(1)) M^3`. So adding
residue strengthening modulo prime powers `<= 2M` for **all** monomials still leaves a ceiling
`>= c (log R)^(1/3)`.

*Proof.* In Theorem A.1(ii) the bad proportion becomes `1.32 Q_M (Pi - 1) <= 1/2`, because
`Pi - 1 <= 2.2 r/sqrt P`. The character part only improves. QED

**Open sharpening.** Fully flipped Paley profiles with `b = O(1)` copies would give the sharp
ceiling `W = Theta(M log M)` claimed heuristically by Prop P, if they are simultaneously
character-admissible. With `b = O(1)` the argument of Lemma A.7 still handles:
* pairs (margin `b - 8`);
* characters with `||c||_1 >= 2.4M`, via the twins.

A fourth-moment Weil argument for `H^T c` should handle `||c||_1 <= M^(1/2)/7` for Paley `H`.
This is a sketch only, not written out and not claimed. Nothing handles `||c||_1` in
`[M^(1/2), 2M]`. What is missing there is an `l^1`-stability
(inverse-uncertainty) theorem for the Paley transform: if `||H^T c||_1 <= (1 + gamma) M`, then `c`
is close to a pair, a column `H_a`, or a two-column half-sum `(H_a +- H_b)/2`. Evidence from the
exact census `check_flat_vectors.py`:
* Paley `M = 20, 24`: no flat character (`F(c) = M`) of support `<= 6` other than pairs, and no flat
  combination of `<= 4` columns other than columns and half-sums.
* Paley `M = 12`: the 220 extra hits of support `6 = M/2` are exactly the `+-(H_a +- H_b)/2`.
* Sylvester: many flat characters of support 4, namely the aligned flats, as Theorem C requires.

In the profile of §2.2 the flip terms rescue columns and half-sums. Each flipped column whose label
`a(j)` has `T_a(j) = 0` contributes `+2|c_x(j)|`. Hence `G(H_a) >= b + 4M - 12`, where
`G = sum_j |c^T S_j| - b(M-1)`.

### 2.3 The strongest inequality from the GAP itself; inverse Littlewood-Offord

**Proposition A.4 (GAP-rank inequality).** Let `rho(W) = max{r >= 0 : r log(4r/e) <= W}`, which is
`(1 + o(1)) W/log W`. Every cluster on an arc of length `<= C sqrt(R)` satisfies
`M <= ceil(C/sqrt 2)(1 + rho(W))`. Hence `M <= (2 ceil(C/sqrt 2) + o(1)) log R/loglog R`.

*Proof.*
1. Split the arc into `ceil(C/sqrt 2)` arcs of length `<= sqrt 2 sqrt R`, and normalise each
   sub-cluster. Its `W_i <= W`, and its arc constant stays `<= sqrt 2`.
2. By `linear_allocation_affine_rigidity.md` (1), each sub-cluster's exponent rows are affinely
   independent, so `M_i <= r_i + 1`, where `r_i` counts its varying split primes. That proof uses
   only the pair gaps: the vectors `f_i` have pairwise negative inner products.
3. The `j`-th prime `= 1 mod 4` is at least `4j + 1`, and `r! >= (r/e)^r`. So
   `W >= W_i >= sum_(j <= r_i) log(4j+1) >= r_i log(4 r_i/e)`.

This is checked for `r <= 10^5` in `check_energy_rank.py`. QED

This is the best inequality obtained from the GAP structure alone. It recovers the growth rate
`log R/loglog R` with constant `2 ceil(C/sqrt 2)`, against `2/3` in `two_thirds`. It is a pair-gap
argument, so it lies in the class of Corollary A.4.

**Proposition A.5 (inverse Littlewood-Offord is vacuous).** Take `C <= sqrt 2`. Let `rho = M_1/|B|`
be the concentration probability of the cluster's unit class under the uniform random walk on the
box `B`, where `M_1 <= M`. If `rho >= r^(-A)` for some `A >= 1`, then `r <= r_0(A) = O(A log A)` and
`M <= r_0(A) + 1`.

*Proof.* `|B| = prod(e_j + 1) >= 2^r` and `M <= r + 1`, so `r^(-A) <= (r+1) 2^(-r)`. QED

So the polynomial-concentration hypothesis of every inverse Littlewood-Offord theorem (Tao-Vu,
Nguyen-Vu) holds only for bounded clusters. In the unbounded regime `rho` is super-polynomially
small. There the small-ball bounds of Rudelson-Vershynin carry an additive `e^(-c r)` floor, which
gives at best `M <= |B| e^(-c r)`. This is exponentially large unless `c > log 2`. The "weighted
rare-fiber inverse theorem" asked for in `deep_global_methods_second_pass.md` §7 would be a statement
about angle models. By Corollary A.4 it cannot be proved from (I)-(IV) below the square-root
ceiling.

### 2.4 A remark on residues outside the character span

Theorem A.3 adds residue strengthening for all monomials only at moduli `<= 2M`. With residue data
at **every** modulus, an angle model would have to satisfy, for all `v` (taking the real-axis
version), `|sin <v,phi>| e^(w(v)/2) >= Q_v`. Here `Q_v` is the part of `Im A_v` detected by the
residues. If the residues are those of actual Gaussian primes, then `Q_v = |Im A_v^actual|`, and
the condition reads `|sin <v,phi>| >= |sin <v, phi_actual>|` for every `v`. It holds with equality
at the actual angles. I do not know whether it determines `phi`, and claim nothing beyond the
following heuristic. Using it amounts to using the actual primes, which is the original problem.

For a heuristic picture with random residues at all primes `q <= Q`: the expected number of
violated monomials is about `(Pi - 1) 2^(pi(Q))`. Random completions therefore survive only for
`pi(Q) <~ log_2(M)`. This is the one place where "monomial integrality" carries information beyond
characters. But it constrains only angles **not** determined by the cluster. It couples to the
cluster only through monomials `v_out` that are themselves near-resonant, and nothing forces those
to exist. I found no mechanism by which it could exclude a cluster, and none is claimed.

---------------------------------------------------------------------------------------------

## 3. Direction (b): sum-product and multiplicative energy

**Proposition B.1.**

(i) For any `C`, the cluster is additively Sidon in `Z^2`. If `z_a + z_b = z_c + z_d = s` on a
circle centred at `0`, then `{z_a, z_b}` are the two intersections of `|z| = R` with `|s - z| = R`.

(ii) For `C < 2 sqrt 2` it is multiplicatively Sidon (`multiplicative_rectangle_separation.md`).
For `C <= sqrt 2` it is multiplicatively independent modulo units
(`linear_allocation_affine_rigidity.md` §3).

(iii) For every quadruple,

```text
z_a z_b / (z_c z_d) = (unit) A_v / conj(A_v),   with v = a_a + a_b - a_c - a_d in L.
```

Likewise `z_i conj(z_j) = (unit) Norm(g_ij) A_ij^2`, which places the products near the real axis.
So every approximate multiplicative-energy statement, at any scale, is a statement about the
characters `e_a + e_b - e_c - e_d` or `e_i - e_j`.

Sum-product inequalities bound `max(|A+A|, |A.A|)` from below, or energies from above. By (i) and
(ii) both are already extremal, so these inequalities are vacuous. By (iii) every approximate
version is a character argument, covered by Corollary A.4.

The products `z_i conj z_j` are lattice points of the circle of radius `N` in the axis strip
`|Im| <= C N^(3/4)`. Chan's near-axis theorems need the narrower strip `N^(1/2)(log N)^(1/7)`.
The notes' `chan_almost_square_axis_audit.md` already explains why they do not apply.

**Proposition B.2 (four-term characters are subcritical on average).** Take a 0/1 column with
`S` ones among `M` rows, and sum over all ordered quadruples in `[M]^4`. Then exactly

```text
sum |x_a + x_b - x_c - x_d| = 4 S (M - S)(M^2 - S M + S^2),     sum_(pairs) |x_a - x_c| * M^2 = 2 S (M-S) M^2.
```

So, column by column, the quadruple average of `w(a+b-c-d)` is `2(1 - s + s^2) >= 3/2` times the
pair average (`s = S/M`). The pair average is at most `W/2` by the cut identity, which is exactly
critical. Hence the 4-term characters have average weight `>= (3/4) W` on balanced binary profiles:
an average slack of `W/4` over their requirement `W/2 - O(log C)`.

For multi-level columns one still has `E|D + D'| >= E|D|` (Jensen, `D'` symmetric). The exact search
in `check_energy_rank.py` found the minimal ratio `127/90`.

*Proof.* Let `u = x_a + x_b` and `u' = x_c + x_d`, each in `{0,1,2}` with counts
`(M-S)^2, 2S(M-S), S^2`. The sum of `|u - u'|` over the 9 cases equals the stated polynomial.
This is checked exactly for all `M < 40` and `0 <= S <= M`. QED

So energy-type arguments gain nothing on average. Only aligned structures can be critical, namely
those whose `L^1-L^2` defect `Delta_c` is small (walsh-extraction Theorem A). That is precisely
the Walsh / aligned-flat mechanism already in the notes. Direction (b) is a subclass of the
character class.

---------------------------------------------------------------------------------------------

## 4. Direction (c): incidence geometry

**Proposition C.1 (triangle identity).** Let `z_a, z_b, z_c` be distinct with `|z|^2 = N`, and put
`g_ab = gcd(z_a, z_b)`, `g_abc = gcd(z_a, z_b, z_c)`, `nu_ab^2 = |z_a - z_b|^2 / Norm(g_ab)`, which
is an integer, and `2 Area = |det(z_b - z_a, z_c - z_a)|`. Then

```text
Norm(g_ab) Norm(g_bc) Norm(g_ac) = N Norm(g_abc)^2,      4 (2 Area)^2 = Norm(g_abc)^2 nu_ab^2 nu_bc^2 nu_ac^2.
```

*Proof.*
1. At a split `p` with allocations `alpha, beta, gamma`, the `p`-exponent of the left-hand product
   is `3e - (|alpha-beta| + |beta-gamma| + |alpha-gamma|) = 3e - 2(max - min)`. That equals the
   exponent of `N Norm(g_abc)^2`. Inert and ramified primes have the same exponent in all three
   points.
2. Then `2 Area = |z_a-z_b||z_b-z_c||z_a-z_c|/(2R)` (circumradius formula). Substitute
   `|z_a - z_b| = |g_ab| nu_ab`. QED

`check_triangle_identity.py` verifies this exactly on 16,594 triples, from 401 circles and 1,172
clusters with `C <= 8`, including the 10-point cluster on `n = 1176852625`.

So the Jarnik integrality `2 Area >= 1` and its refinements are products of three pair identities.
Summing the logarithm over all triples reproduces `(M-2)` times the summed pair identity (2.3) of
`two_thirds`, because of the triple cut identity
`sum_triples log Norm g_abc = binom(M,3) W - (M-2)/2 sum_pairs d_ab`. Triangle-area and
distance-integrality arguments are therefore pair/local.

**Proposition C.2 (lattice-line incidences).** Let `v in Z[i]` be primitive, and let `Gamma` be an
arc of `|z| = R` of angular width `D < pi/2` and length `L = RD`. Let `eta` be the angular distance
from `Gamma` to `arg v + pi Z`. Then

```text
#(Gamma cap Z^2) <= |v| L^2 / R + 2            if eta = 0,
#(Gamma cap Z^2) <= |v| L (eta + L/(2R)) + 1   if eta > 0.
```

*Proof.*
1. Every lattice point lies on exactly one line `Re(z conj v) = k`, `k in Z`. Each line meets the
   circle in the two points at angles `arg v +- alpha`.
2. If `eta > 0`, the arc contains at most one of the two points. Otherwise it would contain an arc
   joining them, hence `arg v` or `arg v + pi`.
3. On `Gamma`, `Re(z conj v) = R|v| cos(psi - arg v)` ranges over an interval of length at most
   `R|v|(cos eta - cos(eta + D)) <= |v| L (eta + D/2)`. If `eta = 0` the length is at most
   `R|v|(1 - cos D) <= |v| L^2/(2R)`, with two points per line. QED

*Consequences.*
* If the arc is within `K R^(-1/2)` of a lattice direction `v`, then `M <= C(K + C)|v| + 2`. This is
  a uniform bound for that class of arcs. The axes and diagonals give `M <= C^2 + 2` and
  `M <= sqrt 2 C^2 + 2`. The 10-point cluster sits on the diagonal and gives the bound 98.
* For a general arc, Dirichlet's theorem on lattice directions gives `M <= 2 C^(3/2) R^(1/4) + 2`.
  This is worse than the divisor bound.

`check_lattice_lines.py` checks C.2 on 3,484 clusters against all primitive `|v| <= 40`. Its
docstring describes the floating-point use; this is a sanity check, not a proof.

**Remark C.3 (Szemeredi-Trotter, Pach-Sharir, unit distances).** These theorems bound the number of
**rich** curves, or of incidences, in large families. A counterexample needs a single rich circle
per radius.

* The incidence structures a cluster generates have incidence counts inside the
  Szemeredi-Trotter / Pach-Sharir ranges. Example: the `M(M-1)` difference points `z_a - z_b` lie on
  the `M` translated circles `|z + z_b| = R`. Then `I = M(M-1) <= c((M^2)^(2/3) M^(2/3) + M^2)`.
* Within a cluster all chord lengths are distinct, by the Sidon property, so unit-distance counts
  are trivial.
* Real-plane incidence arguments also hold verbatim for the fake circles of Theorem A.3, which have
  the same real geometry.
* What lattice integrality adds is distance and area integrality, which is pair-local by C.1. Exact
  low-height lattice conics through many points are auxiliary polynomials, barred by
  `threshold/upper.md`.

So direction (c) yields only C.2, which is uniform on a thin class of arcs.

---------------------------------------------------------------------------------------------

## 5. Direction (d): equidistribution along subsequences

**Proposition D.1.** For `X >= 1` and `T := C X^(1/4) >= 1`,

```text
#{ N <= X : two distinct lattice points of x^2+y^2=N at distance <= C N^(1/4) }
     <= 16 C X^(3/4) (1 + log(C X^(1/4))) + 14 C^2 X^(1/2).
```

So for every `M >= 2` and every `C`, the `N` admitting a cluster at scale `C` have density zero.

*Proof.*
1. Map such an `N` to `(z', w)` with `w = z - z' != 0`, `|w| <= T` and `|z'| <= sqrt X`. Then `z'`
   lies on the line `2 Re(u conj w) = -|w|^2`.
2. Write `w = g w_0` with `w_0` primitive. The lattice points of that line are spaced `|w_0|` apart,
   so at most `2 sqrt X/|w_0| + 1` of them have `|u| <= sqrt X`.
3. Sum over `g <= T` and `0 < |w_0| <= T/g`, using `sum_(0<|w|<=S) 1/|w| <= 8S` (max-norm shells) and
   `#{0 < |w| <= S} <= 8 S^2`. QED

Exact counts in `check_density.py`, chord version (which is `>=` the arc version and is also bounded
by the proof):

* `X = 10^7`, `C = 1`: 179,581 values of `N` have a close pair (`= 1.01 X^(3/4)`), 204 have 3 points,
  none has 4.
* `X = 10^7`, `C = 2`: 360,572; 4,585; 32.

The data suggest `E_2 ~ C X^(3/4)` and `E_3 ~ X^(1/2)`. This is evidence, not proof.

**Proposition D.2 (clusters are invisible to weak limits).** Take `C <= sqrt 2` and a cluster of `M`
points on `x^2 + y^2 = N`. The proportion of `E_R` it occupies is at most `M 2^(-M-1)`.

*Proof.*
1. After normalisation, `M <= r + 1` by Prop A.4's input, so `r_2(N') >= 4 * 2^(M-1)`.
2. Multiplication by the removed gcd injects `E_(R')` into `E_R`. QED

So fixed Fourier coefficients and weak limits of `mu_N` change by `O(M 2^(-M))` under insertion or
deletion of a cluster. The classification of attainable limits (Kurlberg-Wigman; audited in
`kurlberg_wigman_shrinking_arc_audit.md`) is blind to clusters. Concentration results along
density-zero subsequences act on fixed scales, while D.1 already confines clusters to density zero.
Any "arithmetic rigidity" added on top must hold for each single `N`. It is then a statement about
angle models and falls under Corollary A.4, unless it uses one of the inputs of Section 6.

---------------------------------------------------------------------------------------------

## 6. What is not covered, and why I did not find an argument there

Corollary A.4 covers every argument whose arithmetic input is integrality of monomials, characters
with small-modulus residues, exponent identities and generic properties of the angles. An argument
outside it must use one of the following.

1. **Exact algebraic identities among the points** (polynomials in the `z_x`, Ptolemy relations,
   additive relations such as `z_a + z_b = 2 Re(A_ab) g_ab`). Single-place auxiliary polynomials are
   capped at `tau*(M) < 2` (`threshold/upper.md`). Forced Gaussian divisibility is only critical
   (`threshold/construct.md`). The Ptolemy route needs `(L_8)` (`threshold/L8.md`; `round4/L7.md`).
2. **The algebraicity of `e^(2 i phi_j) = pi_j/conj pi_j`** beyond Lemma 0.1, that is,
   transcendence measures. Baker/Matveev bounds miss by a factor `~ c log M`
   (`round4/flipped.md`, Prop B); Lang-Waldschmidt would suffice for part of the Paley class (ibid.,
   Thm L). This is conditional.
3. **Residues at all moduli for monomials outside `L`** (Section 2.4). This holds with equality at
   the actual angles, so in effect it is the actual primes. A heuristic count shows random
   completions survive only with `pi(Q) <~ log_2 M`. A contradiction would need a **covering**
   statement for the strengthened slabs that genuinely uses the cluster. The same statement without
   the cluster is false, because the actual angles of any actual circle avoid all of its slabs.
4. **Multi-place Diophantine approximation** (Subspace / Ru-Vojta with `S` = primes of `N`). This is
   barred: the exceptional sets depend on `S`.

Two concrete open problems come out of this round:

* **(P1)** Prove fake-realisability of fully flipped Paley profiles with `b = O(1)`. This needs the
  `l^1`-stability theorem of §2.2. It would make the ceiling in Corollary A.4 sharp at
  `Theta(log R/loglog R)`, so the whole local + character + GAP toolkit would be proved unable to
  improve the growth rate.
* **(P2)** Conversely, show that some character family is **not** simultaneously admissible on
  fully flipped Paley profiles at `W = O(M log M)`, for every choice of positions
  `delta in [0, Delta]^M` and winding data `k`. That would prove the fully flipped Hadamard
  phase lemma by a purely real (Diophantine-free) argument about vectors `delta in [0, Delta]^M`.
  Theorem A.1 makes this an honest reformulation: the phase lemma's archimedean content **is**
  simultaneous admissibility.

---------------------------------------------------------------------------------------------

## 7. Files (`round4/outside_checks/`)

| script | what it checks | type |
|---|---|---|
| `gx.py` | Gaussian-integer toolkit (exact representations by factorisation, gcds, cluster finder) | helper |
| `check_completion.py` | Thm A.1(i): `U^perp = span{a_x-a_y}` on 200 random multi-level profiles; constants `1.32`, `p >= 21 r^2` | exact; Part 2 is an illustration |
| `check_fake_barrier.py` | Lemma A.7 (equal-weight forms), saturation, rank, twin structure; Sylvester 8/16/32, Paley 12/20/24, `b = 6, 6M` | exact |
| `check_fake_instance.py` | a concrete fake circle (`M = 8`, 336 primes near `10^9`) | illustration |
| `check_flat_vectors.py` | census of flat characters `F(c) = M` (Sylvester, Paley) for the open sharpening (P1) | exact, evidence |
| `check_energy_rank.py` | Prop B.2 identity; multilevel `E|D+D'| >= E|D|`; the prime bound of Prop A.4 | exact |
| `check_triangle_identity.py` | Prop C.1 on 16,594 triples | exact |
| `check_lattice_lines.py` | Prop C.2 on 3,484 clusters | float sanity check |
| `check_density.py` | Prop D.1 counts for `X = 10^6, 10^7` | exact counts |

Notes searched before claiming novelty (`research/docs/`):
* `fourier_riesz_product_route.md` and `fourier_moment_stress_test.md`: heuristic magnitude barrier only;
* `deep_global_methods_second_pass.md` and `global_methods_broad_audit.md`: heuristic surveys;
* `kurlberg_wigman_shrinking_arc_audit.md`, `vector_box_route.md`, `research_reset_inverse_concentration.md`,
  `linear_allocation_affine_rigidity.md` (used in Prop A.4), `chan_almost_square_axis_audit.md`;
* `growth/walsh-extraction.md` (Prop P: individual character passes; no simultaneous realisation)
  and `round4/flipped.md`.

None contains the completion theorem, a simultaneous realisation (fake circle), or the triangle and
lattice-line identities in the stated form. The density bound D.1 and the near-axis case of C.2 are
elementary and are probably folklore; no novelty is claimed for them.
