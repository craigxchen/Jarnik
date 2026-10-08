# Walsh extraction: from structured-profile growth theorems towards general profiles

Route: items 311--335 (repeated-Walsh growth), the Paley obstruction (items 334, 336, 391--394),
and the question of an extraction lemma transferring item 333 to arbitrary circles.

**Status (honest).** No general growth improvement and no uniform theorem is proved here.
The general bound `M <= (1+eps) log R/loglog R` (inert_prime_cofactor_bound.md) is unchanged.
What is proved (prose proofs below, every identity and every numerical constant checked by exact
scripts in `walsh_extraction_checks/`):

| result | class | conclusion | relation to the notes |
|---|---|---|---|
| Thm A | every two-valued (or threshold-layer) profile | exact "L1-L2 defect" form of every balanced character phase inequality | generalises item 334 (1) from 4 rows to any support; recovers item 322 (2) exactly |
| Thm F | every profile | *residue pigeonhole*: the character inequalities extended to rational (fractional) characters, using the integer unit/winding data `k_x` | new method; vacuous exactly when the rational saturation is half-integral (Lemma F.1) |
| Thm B | pure Hadamard cores of any order and any type (Sylvester, Paley, ...) | for `C<=1`: `M<=24`; for every `C`: `M < max(28, M_1(C))`, `M_1(C)=C^(2+o(1))`, effective, elementary | items 276--277 (Matveev, `C<=sqrt 2` unit regimes, astronomically large thresholds) and the Roth notes; here: every C, small explicit constants, no Baker/Roth |
| Thm B'' | Hadamard cores with up to `M/(40 log M)` deleted rows (equivalently: arbitrary modifications confined to that many rows) | uniform bound, every `C` | item 277 allows only a bounded number `d` of deleted rows |
| Thm C | repeated-Walsh model with flips on an **arbitrary subset** of rows (one flip per flipped row, distinct columns) | `W0/(M log M) -> infinity`, with the item-335 rate; uniformly, unflipped rows are `O(M^(1-1/S_0))` | strictly larger class than items 333/335, which need every row flipped; includes pure Walsh |
| Prop P | fully flipped capacity-two Paley profiles on actual nearby primes | pass pair, **all rational** characters (every Theorem-F consequence, every m), aggregate small-prime collisions, copy capacity, at `W=Theta(M log M)` | sharpens the barrier of items 334/393; pure Paley and few-flip Paley *fail* (Thms B, B'') |
| Lemmas 7.2-7.3 | fully flipped Hadamard cores | exact "flip-shifted pinning" form of phase localisation; rigidity (unflipped blocks determined by the flipped primes and `M-1` residues) | precise form of the remaining Diophantine problem |

Part (2) of the task: the extraction lemma is formulated precisely (Section 8), its conditional
consequence is proved, and it is tested on actual clusters (Section 9). The data show **no
enrichment of approximate (flipped) Walsh structure** in tight point sets beyond the
binomial chance rate, and an enrichment of *exact* multiplicative rectangles that is explained by
a Bohr-set mechanism whose absolute frequency tends to zero on large circles. Together with Prop P
this locates the missing step precisely: a Diophantine phase lemma for fully flipped,
quasi-random (Paley-like) profiles, which no rational-character, collision or capacity argument
can supply.

---------------------------------------------------------------------------------------------

## 1. Setting and notation

After removing the complete Gaussian gcd `g` (content `D=log Norm g >= 0`) the `M` points are

```text
z_x = g eps_x prod_j pi_j^(s_xj),     s_xj in {+1,-1},   pi^(+1)=pi, pi^(-1)=conj(pi),
```

with `eps_x` units and `pi_j` Gaussian primes over distinct split primes `p_j` (two-valued model;
Theorems B, B'', F only use *blocks* `G_a` = products of such factors with pairwise disjoint rational
prime supports, each conjugate-primitive, so nested threshold layers grouped into whole classes are
covered exactly as in item 257). `w_j=log p_j`, `W=sum_j w_j`, `log R^2 = W + D`.
The points lie in an arc of angular width `Delta <= C R^(-1/2) = C exp(-(W+D)/4)`.

Pair data, as in the three-family and pair-slack notes:

```text
G_xy = sum_j w_j s_xj s_yj,   kappa = 4 log C - log 4 - D,   u_xy = kappa - G_xy >= 0.
```

**Winding integers.** Fix real lifts `phi_j` of `arg pi_j`. Since `arg z_x` lies in an interval
`[theta*-Delta/2, theta*+Delta/2]` modulo `2 pi`, there are integers `k_x` (absorbing the units) with

```text
(S phi)_x := sum_j s_xj phi_j = theta' + delta_x + (pi/2) k_x,    |delta_x| <= Delta/2.     (1.1)
```

For rational `lambda` with `sum lambda_x = 0` and integral column image `v = lambda^t S`, put
`beta(lambda) = prod_j pi_j^(v_j)` (conjugate factor for negative exponent), `V(lambda)=sum_j |v_j| w_j
= log Norm beta`. Then exactly

```text
arg beta(lambda) = (pi/2) lambda^t k + lambda^t delta  (mod 2 pi),   |lambda^t delta| <= ||lambda||_1 Delta/2.   (1.2)
```

**Lemma 1.1 (axis gap).** A conjugate-primitive nonunit Gaussian integer `beta` of odd norm is not a
unit times a rational integer nor a unit times a multiple of `1+i`. Hence if `arg beta` is within `e`
of `(pi/2)Z` then `|beta| sin e >= 1`; within `e` of `(pi/4)Z`, `|beta| sin e >= 1/sqrt 2`.

*Proof.* Rotating by a unit, `beta = X+iY` with `|Y| = |beta| |sin(arg)|`; `Y=0` would make `beta` and
`conj beta` associates, contradicting conjugate-primitivity. Diagonal case: `X=Y` gives `(1+i) | beta`,
impossible for odd norm; `|X-Y| >= 1`. QED

For integral zero-sum `c` (with `v = c^t S/2`) this is exactly the standard half-height inequality
`V_c >= (W+D)/2 + log 8 - 2 log(||c||_1 C)` of the notes.

## 2. Theorem A: the defect form of every balanced character inequality

**Theorem A.** Let `c in {0,+1,-1}^M` be supported on `h` rows `P` with `sum c = 0`, and
`T_j = sum_x c_x s_xj`, `V_c = sum_j (w_j/2)|T_j|`. Then exactly

```text
V_c = W/2 - kappa/2 - (1/(2h)) sum_{x!=y in P} c_x c_y u_xy + Delta_c,
Delta_c = (1/2) sum_j w_j |T_j| (h - |T_j|)/h >= 0.                                    (A.1)
```

More generally for any integral zero-sum `c` and any `Lambda >= max_j |T_j|`,
`V_c = (||c||_2^2/(2 Lambda))(W-kappa) - (1/(2 Lambda)) sum_{x!=y} c_x c_y u_xy
+ (1/2) sum_j w_j |T_j| (Lambda-|T_j|)/Lambda`.

*Proof.* `sum_j w_j T_j^2 = sum_{x,y} c_x c_y G_xy = hW + sum_{x!=y} c_x c_y (kappa - u_xy)
= hW - kappa h - sum_{x!=y} c_x c_y u_xy`, using `(sum c)^2 = 0`. Then
`V_c = (1/2h) sum w_j T_j^2 + (1/2) sum w_j (|T_j| - T_j^2/h)`. QED

**Corollary A.2.** For such a nonunit character the endpoint phase inequality
`V_c >= (W+D)/2 + log 8 - 2 log(hC)` is *equivalent* to

```text
Delta_c - (1/(2h)) sum_{x!=y in P} c_x c_y u_xy  >=  -2 log(h/2).                     (A.2)
```

(`-kappa/2 = -2 log C + log 2 + D/2` cancels every C and D.)

Consequences:

* `h=4`, `c=(1,1,-1,-1)`: `Delta_c = (W - T_1234)/4`, which is item 334 (1).
* Aligned Walsh flat (coset `P=x0+V`, character `(-1)^ell`, all assigned flips in `P` restricting to
  `ell`): exactly `c_x c_y u_xy = (-1)^ell(d) u~_d + 2(f_x+f_y)` with `u~_d = kappa - What(d)` the
  unflipped core slack, `Delta_c = (1-2/h) F_p`, and (A.2) becomes `Psi~ >= 2F_p - 4 log(h/2)`.
  With `u_d = u~_d + (4/M) Fhat(d)` this is *exactly* item 322 (2) for aligned cosets. So the whole
  repeated-Walsh mechanism is: flips make `Delta_c` small while forcing a large *signed* slack sum.
* Any extraction argument must therefore produce characters whose L1-L2 defect `Delta_c` is
  small compared with their signed slack sum. Paley profiles have `Delta_c >= c_0 W` for every
  character in the ranges of items 336/391/393 -- this is the precise content of the Paley barrier.

Checked exactly on 2,999 random weighted profiles (Fractions), 500 four-row specialisations, and 180
literal Gaussian characters on flipped Paley-12/Walsh-8/Paley-20 rows with content `2+i` and random
units (`check_defect_identity.py`).

## 3. Theorem F: residue pigeonhole for rational characters

**Theorem F.** Let `U` be the row set of a sub-cluster (possibly all rows), `S_U` its sign matrix,
`m >= 1`, and `B` a finite subset of `(1/m) Z^U` such that for all distinct `lambda, lambda'` in `B`,
`mu = lambda - lambda'` satisfies: `sum mu = 0`; `mu^t S_U` is integral and nonzero on some
positive-weight column; and

```text
V(mu) < (W+D)/2 - 2 log(||mu||_1 C / 2).                                              (F.1)
```

Then the residue map `rho(lambda) = m lambda^t k_U  mod m` is injective on `B`; hence `|B| <= m`.
If moreover every nonzero `lambda_1+lambda_2-lambda_3-lambda_4` (elements of `B`) satisfies the same
three conditions, then `rho(B)` is a Sidon set in `Z/m` and `|B| <= (1+sqrt(4m-3))/2`.

*Proof.* If `rho(lambda)=rho(lambda')` then `mu^t k` is an integer, so by (1.2) `arg beta(mu)` lies
within `||mu||_1 Delta/2` of `(pi/2)Z`. Lemma 1.1 gives `V(mu) >= 2 log(2/(||mu||_1 Delta))`, i.e.
the reverse of (F.1). The Sidon statement is the same argument for 4-term combinations; in a Sidon
set of `Z/m` all ordered differences of distinct elements are distinct and nonzero. QED

For `m in {1,2}` Theorem F is the usual character inequality. It is strictly stronger exactly when
the profile has *fractional saturation*: rational `lambda` with integral column images and
denominators larger than two.

**Lemma F.1 (saturation; items 321, 393).** If every row of `U` carries a flip at a physical column
that has an unflipped sibling (same label, same core values), then every rational `lambda` with
`lambda^t S_U` integral has `2 lambda` integral. *Proof:* the flipped column minus its sibling is
`+-2 e_x`. QED

So for fully flipped profiles every Theorem-F consequence, for every `m`, reduces to a half-height
inequality for an integral character. This is the structural reason why flips are essential both to
item 333 (they create the aligned-flat contradiction) and to the Paley barrier (they kill fractional
characters).

## 4. Theorem B: pure Hadamard cores (any type, any arc constant)

Let `H` be a normalised Hadamard matrix of order `M` (first column all ones; the other columns have
zero sum and `H^t H = M I`). A *pure Hadamard-core cluster* is

```text
z_x = g eps_x prod_a G_a^(H(x,a)),   x = 1..M,   a over the nonconstant columns,
```

with `G_a` conjugate-primitive (possibly `1`) with pairwise disjoint rational prime supports; this is
the merged-block form of the pure (unflipped) repeated-Walsh (`H` Sylvester) and pure Paley profiles
with any number of copies per label. `W_a = log Norm G_a`, `A = {a : W_a > 0}`.

**Lemma B.1 (label count).** Let `U` be any set of `u >= 2` rows and suppose
`G_xy = sum_a W_a H(x,a) H(y,a) <= gamma` for all distinct `x,y in U`. Then

```text
|A| >= u W / (2W + 2 u gamma_+).
```

*Proof.* `sum_{x,y in U} G_xy = sum_a W_a (sum_{x in U} H(x,a))^2 >= 0`, so `sum_{x!=y} G >= -uW`.
Orthogonality-free lower bound: `sum_{x,y in U} G_xy^2 = sum_{a,b} W_a W_b (sum_{x in U} H(x,a)H(x,b))^2
>= u^2 sum_a W_a^2`, so `sum_{x!=y} G^2 >= u^2 sum W_a^2 - u W^2`. Since `G <= gamma_+`,
`|G| <= (gamma_+ - G) + gamma_+`, whence `sum_{x!=y}|G| <= 2 gamma_+ u^2 + uW`, and `|G| <= W`
gives `sum_{x!=y} G^2 <= W(2 gamma_+ u^2 + uW)`. Hence `sum W_a^2 <= 2W^2/u + 2 gamma_+ W`; finish by
Cauchy-Schwarz `W^2 <= |A| sum W_a^2`. QED

The pair condition `|z_x - z_y| >= |gcd(z_x,z_y)|` gives `G_xy <= gamma := 4 log C - D` for all pairs.

**Theorem B.** Let `M >= 28` be a Hadamard order, `s_M = floor((1+sqrt(4M-3))/2) + 8`. A pure
Hadamard-core cluster of order `M` on an arc of length `<= C sqrt R` must satisfy

```text
W <= max( 8 M s_M log^+ C / (M - 2 s_M),  36 log^+(2C) ).                              (B.1)
```

*Proof.* Pinning. By (1.1) and `H^t H = M I`, for every nonconstant `a`
`phi_a := arg G_a = (1/M) sum_x H(x,a)(S phi)_x = (pi/(2M)) m_a + eta_a`, with
`m_a = sum_x H(x,a) k_x` and `|eta_a| <= Delta/2` (the `theta'` term cancels since columns are
balanced). Equivalently, in Theorem F take `B = {lambda_a = H(.,a)/M}`: `lambda_a - lambda_b` has zero
sum, column image `e_a - e_b`, `||.||_1 = 1`, and 4-term combinations have `||.||_1 <= 4`.

Light labels: `W_a < theta_L := [(W+D)/2 - 2 log(2C)]/4`. Any nonzero 2- or 4-term combination
`mu` of light labels has `V(mu) < 4 theta_L <= (W+D)/2 - 2 log(||mu||_1 C/2)` (as `||mu||_1 <= 4`),
i.e. it satisfies (F.1). Theorem F (Sidon form, `m = M`) gives `#light <= floor((1+sqrt(4M-3))/2)`.
Heavy labels: if `W > 36 log^+(2C)` then `theta_L > W/9` and there are at most 8. So
`|A| <= s_M`. Lemma B.1 with `u = M`, `gamma_+ <= 4 log^+ C` gives `|A| >= MW/(2W + 8M log^+C)`.
If `W` exceeded both terms of (B.1), these two bounds would contradict each other. QED

**Corollary B.2 (`C <= 1`).** Every pure Hadamard-core cluster on an arc of length `<= sqrt R` has
`M <= 24`. *Proof.* `gamma_+ = 0`, so `|A| >= M/2`; for `M >= 28` at least 14 distinct split primes
occur, so `W >= log(5*13*...*113) > 52 > 36 log 2`, and `M/2 > floor((1+sqrt(4M-3))/2) + 8` for every
Hadamard order `M >= 28` (table in the checker output: 28 and 32 are the first excluded orders). QED

**Corollary B.3 (every `C`).** By the audited inert-prime inequality
`F(M) <= (M/4) log R' + binom(M,2) log C` applied to the primitive normalisation (constant `<= C`),
`W >= 8F(M)/M - 4(M-1) log^+ C = (2+o(1)) M log M - 4M log^+ C`, which exceeds the right side of
(B.1) once `log M > (2+o(1)) log C`. Hence pure Hadamard-core clusters have `M < max(28, M_1(C))` with
`M_1(C) = C^(2+o(1))`. Numerically (exact `F(M)`, all powers of two up to `2^17`) the hypothesis of
Theorem B holds automatically from `M = 64` (C=1), 128 (C=2), 256 (C=4), 1024 (C=8), 4096 (C=16).

*Scope and comparison.* The pinning `phi_a in (pi/2M)Z + O(Delta)` is the Walsh/Hadamard isolation of
affine_cube_walsh_rigidity.md, fixed_hadamard_roth_phase_obstruction.md and items 257, 276--278,
which conclude via Roth (ineffective, large N) or Matveev (effective but with `2^40`-size thresholds,
unit regimes `C <= sqrt 2`). The new step is to compare *different* labels through their residues
`m_a mod M` (pigeonhole/Sidon in `Z/M`) instead of approximating a single block to a fixed ray. This
needs no transcendence input, works for every `C`, applies to every Hadamard type including Paley,
and gives `M <= 24` at `C <= 1`. Pure Paley profiles pass every *integral* character inequality
(`B(c) >= bM/2 > r/2` by Hadamard inversion), so Theorem B is a genuinely non-integral obstruction.

## 5. Theorem B'': row-deleted Hadamard cores

**Theorem B''.** Let `M >= 4096` be a Hadamard order, `F` a set of `f <= M/(40 log M)` rows, `U` the
others, and suppose the sub-cluster `{z_x : x in U}` is a pure Hadamard-core cluster restricted to the
rows `U` (blocks as above). Put `L = floor(M/60)`. If
`W >= max(8M log^+ C, 32 log^+(LC))` (automatic for `M >= M_4(C) = C^(O(1))` by the inert-prime
bound), the configuration does not exist. In particular, for any Hadamard core and *any*
modifications (any number of flips, any columns) confined to at most `M/(40 log M)` rows, the
unmodified rows already violate the endpoint hypothesis.

*Proof.* Lemma B.1 with `u = M-f`, `gamma_+ <= 4 log^+ C`: `|A| >= (M-f)/3`. Light: `W_a < theta`,
`2L theta = (W+D)/2 - 2 log(LC)`; the hypothesis gives `theta >= 7W/(32 L)`, so at most `32L/7` heavy
labels and `n >= (M-f)/3 - 32L/7 >= 0.2558 M` light ones. For `v in {0,1}^light` with `|v| = L` record
`((Hv)_x)_{x in F} in [-L,L]^f` and `sum_{x in U} (Hv)_x k_x mod M`. There are `binom(n,L) >= 15.3^L`
such `v` and at most `(2L+1)^f M` records; `15.3^L > (2L+1)^f M` because
`0.0455M - 2.73 > M/40 + log M` for `M >= 4096` (exact integer comparisons at several orders in the
checker). Two `v != v'` share a record. Then `mu = H(v-v')/M` vanishes on `F`, has zero sum, column
image `v - v'` on `U` (since `H^t H = M I`), `||mu||_1 <= 2L`, and `mu^t k in Z`. Its height is
`V(mu) < 2L theta <= (W+D)/2 - 2 log(||mu||_1 C/2)`, so `mu` satisfies (F.1) while
`rho(v) = rho(v')`, contradicting Theorem F. QED

Item 277 (deleted_hadamard_uniform_family_bound.md) proves the bounded-deletion case with Matveev in
unit regimes; the residue pigeonhole allows `M/(40 log M)` deletions and every `C`. Exactness of the
certificates (support on `U`, zero sum, column image) is verified on 90 random collisions for
Walsh-32/64 and Paley-44 (`check_theorem_C_ingredients.py`).

## 6. Theorem C: repeated-Walsh profiles with an arbitrary set of flipped rows

**Model WF(t,b,F).** Rows `x in F_2^t`, `M = 2^t`; `b >= 5` physical columns (distinct split primes)
for each nonzero label; core signs `sigma_j (-1)^(a_j . x)`; an *arbitrary* set `F` of rows, each
`x in F` flipped at one physical column `j(x)` of nonzero label `a(x)`, with `j` injective (copy
capacity `b`); rows outside `F` unflipped. Arbitrary prime sizes, orientations, row units and common
content. Item 333's model is `F = F_2^t`; `F = empty` is the pure Walsh profile.

**Theorem C.** Fix `C`. For WF(t,b,F) configurations on arcs of length `<= C sqrt R`,
`W0/(M log M) -> infinity` as `M -> infinity`, uniformly in `b, F`, prime sizes and units; more
precisely `W0 >= c M log M sqrt(T)/log_2 T` and `M = O_C(log R log_2 T_R /(loglog R sqrt(T_R)))`,
`T = log_2^* M` (the rate of item 335). Moreover, uniformly in `R`, the unflipped set satisfies
`|F_2^t \ F| < 2 S_0^(1/S_0) M^(1 - 1/S_0)`, where `S_0 = S_0(C)` is the least power of two
`>= max(32, M_1(C))` (`S_0 = 32` for `C <= 1`).

*Proof.* (i) Density lemma (item 326, re-proved and machine-checked): a subset of `F_2^t` of density
`>= 2(K/M)^(1/K)` contains an affine flat of size `K`. (Doubling: among `n >= 2` disjoint cosets of
`V`, some direction occurs for `>= n^2 h/(2M)` ordered pairs; the pairs `{P, P+e}` are disjoint and
give cosets of `V + <e>` with density `>= rho^2/2`.)
(ii) If the unflipped set `U` has density `>= 2(S_0/M)^(1/S_0)`, it contains an affine `S_0`-flat `Y`.
On `Y` no column is flipped, and every column is an affine character of `Y` (constant ones are
content). So `{z_x : x in Y}` is a pure Sylvester-core cluster of order `S_0` with arc constant
`<= C` after primitive normalisation, contradicting Theorem B / Corollaries B.2-B.3.
(iii) Otherwise `|U| < 2 S_0^(1/S_0) M^(1-1/S_0)`. Run the proof of item 333 (and its item-335
parameters) verbatim, with `F` the total flipped mass and the rows of `U` added to the bad rows. All
identities used there (pair-slack Fourier identity, (4)-(6), the four-row bound
`u_d <= 6W0/M + 2a_C + 4 log 2`, the Ramsey colouring) only involve flipped masses `F_a` and hold for
any `F` (unflipped rows contribute zero). In the bad-row average over the `M/N` cosets of the Ramsey
subspace `H` (`N = 2^s` fixed, or `N <= 2^(log^* M)` in the quantitative version) the new term is
`N|U|/M <= 2 N S_0^(1/S_0) M^(-1/S_0) -> 0`, so a coset with at most `b+1` bad rows still exists;
patching, extraction and discarding are unchanged, and the surviving aligned flat consists of flipped
rows with flipped prime logs `>= (log M)/2`. The final contradiction is item 333 (11). QED

The uniform part (ii) is the `k = 0` case of item 278 (dense Walsh restriction), here with the
elementary threshold `S_0 = 32` at `C <= 1`; the new content of Theorem C is that the dichotomy
"dense unflipped part / almost all rows flipped" closes the gap between pure Walsh profiles and item
333, so the growth improvement holds for every flip pattern. Exact checks: density-lemma flats in 44
random dense sets; for a partially flipped Walsh-64 with an unflipped affine 16-row cube, all 60
fractional characters `chi_a|_Y/16` have integral images (`check_theorem_C_ingredients.py`).

## 7. Question (1): which constraints do Paley-type profiles fail?

### 7.1 The sharpened barrier

**Proposition P.** Let `q = 3 mod 4` be a large prime, `M = q+1`, `b >= 5`, and take the fully
flipped capacity-two Paley profile (one distinct flipped copy per row, each label flipped at most
twice) on distinct split primes with `tau <= log p_j < tau + 1/r`, `tau = 4 log M + O(1)`, arbitrary
orientations and units (`W = (4b+o(1)) M log M`). Then:

1. every pair inequality holds (item 334);
2. every integral character inequality holds: every nonzero half-integral character has
   `V - W/2 > (tau-1)/2` (items 336 (12), 391, 393; for the canonical assignment the margin is
   `> tau(b/2-2) - 1/2`, item 392), far above the required `log 8 - 2 log(||c||_1 C)`;
3. every **rational** character inequality, i.e. every consequence of Theorem F for every `m`, holds:
   by Lemma F.1 all rational characters are half-integral, so this reduces to 2;
4. the aggregate small-prime collision inequality used by the general bounds holds: every column has
   `S_j(M-S_j) >= M^2/4 - 1`, so `sum_{x<y}(d_xy - W/2) >= W(M/4 - 1) = (b+o(1)) M^2 log M`, which
   exceeds the forced collision mass `(3/4+o(1)) M^2 log M` (two-thirds sketch; inert part
   `(1/4+o(1))M^2 log M` audited) for every `b >= 1`;
5. copy capacity holds by construction.

In contrast:

* **pure Paley** (no flips) passes 1, 2, 4, 5 but **fails 3**: Theorem B excludes it for
  `M >= max(28, M_1(C))`;
* Paley with flips on at most `M/(40 log M)` rows **fails 3** (Theorem B'');
* flipped Walsh fails the aligned-flat characters (items 333/335; Theorem C for partial flips).

Hence, within Hadamard-core profiles, every argument built from pair, character (integral or
rational), aggregate collision and capacity inequalities is blocked *exactly* on the class
"non-Walsh core with at least `M/(40 log M)` (in particular all) rows flipped". This answers the
"ideally near-uniform slack" question negatively for these methods: the profiles with the most
uniform pair slack (pure cores, slack constant up to `O(1)`) are excluded by Theorem B, while the
barrier (fully flipped Paley) has slack uniform up to `4 tau + 1 = O(W/M)` and passes everything
rational. Near-uniformity of slack is not the dividing line; fractional saturation is.

### 7.2 Small-prime residue collisions, pair by pair

The aggregate statement 4 is the only form used by the general bound. The per-pair form
`coll_xy := sum_{q <= M} v_q(Norm(u_x - u_y)) log q <= slack_xy` is necessary for an arc but is
*implied by* (and in strength close to) smallness of the reduced chord `u_x - u_y`, i.e. by the arc
condition itself. Literal Paley realisations (`paley_collision_check.py`, `C = 1`, `D = 0`, primes from
the densest window of log-width `1/r`):

| order | b | max coll (q<=M) | max(coll - slack) | violating pairs |
|---|---|---|---|---|
| 12 | 5 | 18.66 | +4.15 | 6/66 |
| 20 | 5 | 24.48 | +16.00 | 21/190 |
| 44 | 5 | 37.37 | +29.79 | 95/946 |
| 12 | 9 | 17.24 | -20.43 | 0/66 |
| 20 | 9 | 19.98 | -17.60 | 0/190 |
| 44 | 9 | 34.43 | -3.75 | 0/946 |

So at these orders the per-pair inequalities can be met by taking more copies (still
`W = Theta(M log M)`), but the margin shrinks. Heuristic only: for a non-arc realisation the reduced
chords are integers of size `exp(Theta(W))` whose `M`-smooth parts behave like those of random
integers; the maximum over `M^2` pairs of `log` of the `M`-smooth part is then about
`2 (log M)^2/loglog M`, which beats the slack `~ 2(b-4) log M` for any fixed `b` once `M` is huge. Thus
per-pair collisions are probably *not* jointly satisfiable by a Paley realisation at `W = Theta(M log M)`
asymptotically, but a proof would need exactly the arithmetic of actual reduced chords, i.e. the arc
condition; I found no way to turn this into an independent rigorous inequality. Per-pair collisions
restricted to `q <= 0.9 (b-4) log M` are satisfied for any residues (the minimal slack is
`((b-4) tau - 1)/2 = 2(b-4) log M + O(1)` by item 334, and `2 theta(Q) < 2.04 Q`), but they carry only
a `loglog M` share of the forced collision mass.

### 7.3 Phase localisation of all rows: exact form for flipped Hadamard cores

**Lemma 7.2 (flip-shifted pinning).** In a Hadamard-core profile with flipped row set `F` (row `x`
flipped at prime `pi_{j(x)}` of label `a(x)`, angle `psi_x`), let `Phi_a` be the summed angle of all
primes of label `a` (oriented as the core) and `sigma_x = H(x,a(x))`. Then for every nonconstant `a`

```text
Phi_a - (2/M) sum_{x in F} H(x,a) sigma_x psi_x  in  (pi/(2M)) Z + [-Delta/2, Delta/2].
```

*Proof.* `(S phi)_x = sum_a H(x,a) Phi_a - 2 sigma_x psi_x 1_{x in F}`; invert with `H^t H = M I` in
(1.1). QED

For `F` empty this is the pinning of Theorem B. A relation `sum n_a Phi_a` is free of the `psi`'s iff
`(Hn)_x = 0` for all `x in F` -- exactly the certificates of Theorem B''. For `F` = all rows only
`n = 0` remains, so the residue method has nothing to pigeonhole: the targets of the `M-1` block
angles move linearly with the `M` flipped-prime angles.

**Lemma 7.3 (sector uniqueness and rigidity).** If `X alpha < 1`, a sector of angle `alpha` contains at
most one conjugate-primitive Gaussian integer of norm `<= X` (two non-collinear lattice points and
`0` would span a triangle of area `< 1/2`; collinear points on a ray are integer multiples of one
primitive point). Consequently, in a fully flipped Hadamard-core endpoint configuration whose
unflipped block of label `a` has `log Norm < (W+D)/4 - log C` (true for nearby primes, where it is
`~ (b-1)tau << W/4`), the `M` flipped Gaussian primes together with the residues `m_a mod 4M`
determine every unflipped block uniquely: at most `(4M)^(M-1)` completions per choice of flipped
primes.

Heuristically each of the `M-1` targets is hit with probability `~ e^(W_a^u) Delta`, so the expected
number of completions is `exp(-(1-o(1)) M W/4)`: phase localisation is the condition that fails, but
it fails "probabilistically", through `M-1` simultaneous Diophantine approximations whose individual
(one-character) consequences are all satisfied (Prop P).

### 7.4 Summary answer to (1)

| constraint | pure Paley | fully flipped Paley | flipped Walsh |
|---|---|---|---|
| pair | pass | pass | pass |
| integral characters, any support/amplitude | pass | pass (393) | fail (aligned flats, 333) |
| rational characters / residue pigeonhole (Thm F) | **fail (Thm B)** | pass (Lemma F.1) | vacuous |
| aggregate small-prime collisions | pass | pass | pass |
| per-pair collisions, `q <= M` | open | open (finite data: pass for b=9 up to order 44) | open |
| copy capacity | pass | pass | pass |
| joint phase localisation | fail (Thm B) | **only remaining candidate**; Lemmas 7.2-7.3 | fail (333) |

## 8. Question (2): the extraction lemma, precisely

**Definition 8.1 (partially flipped Walsh frame of capacity B).** A sub-cluster `Z` of `2^t` points with
a bijection `F_2^t -> Z` such that, after removing `Z`'s common Gaussian content, every rational prime
varying on `Z` takes two exponent levels on `Z`, and there are labels `a_p in F_2^t \ {0}` and
orientations with sign `sigma_p (-1)^(a_p . x)` at every point except a set of flips with at most one
flip per point and per prime; each label carries between 5 and `B` primes.

**Extraction hypothesis EH(alpha, B).** For each `C` there are `B = B(C)`, `M_0(C)` and a function
`alpha_C(M) in (0,1]` such that every primitive endpoint cluster with `M >= M_0(C)` points contains a
partially flipped Walsh frame of capacity `<= B` and size `>= alpha_C(M) M`.

**Theorem 8.2 (conditional).** If EH(alpha, B) holds with
`alpha_C(M) sqrt(T_M)/log_2 T_M -> infinity` (`T_M = log_2^* M`), then `M = o_C(log R/loglog R)`;
quantitatively `M = O_C(log R/(loglog R alpha(M) sqrt(T)/log_2 T))` whenever `alpha >= M^(-1/2)`.

*Proof.* The frame is itself an endpoint cluster on the same arc; its content is retained in `D`.
Theorem C gives `log R^2 >= W0(frame) >= c alpha M log(alpha M) sqrt(T)/log_2 T`; invert as in
item 335. QED

**What EH demands.** Because the item-335 saving is only `sqrt(log^* M)/log log^* M`, EH needs a
single frame containing essentially a fixed proportion of the cluster; many disjoint small frames do
not add (each gives the same lower bound for the same `R`). Literal Walsh-core extraction from a
full-fair profile already loses most of the mass (item 279). In additive language: a frame is a set of
sign vectors in `F_2^k` within one flip per point of an affine subspace, i.e. it has *approximate*
additive energy `~M^3` at defect weight `O(W/M)`, while having no exact additive quadruples.
A fully flipped Paley cluster contains no Walsh frame with more than 3 points (item 334: every
quadruple has quartic defect `(1-o(1))W`, whereas a 2-frame has odd quartic support of at most 4
primes), so EH in particular *implies* the missing lemma of Section 10 for Paley cores: EH is at
least as hard as excluding the approximately-Sidon (quasi-random) class of Prop P. The uniform variant
"every large cluster contains a pure (or `m/(40 log m)`-row-deleted) Hadamard sub-core of order
`>= S_0(C)`" would imply the full uniform theorem by Theorem B''; both statements are, in the vacuous
sense, consequences of the conjecture, and both contain the Paley barrier as their essential case.

## 9. Tests on actual clusters

All points below are literal Gaussian integers on `x^2 + y^2 = N`, `N` a squarefree product of split
primes; clusters are verified with exact integer chords (`cluster_search*.py`).

**9.1 Largest observable clusters are Poisson fluctuations.** On circles with `N >= 10^20` the
smallest arc constants found for `M = 4,...,16` points occur at Poisson means
`lambda = (#points) C/(2 pi sqrt R)` between `0.17` and `7.3` (`large_clusters_s11_L20.json`); on
smaller circles the record clusters are macroscopic (`N ~ 10^8`, arc ~1/9 of the circle). Fixed-`C`
clusters with large `M` cannot be observed (as the conjecture predicts), so all tests concern
structure at small `M`.

**9.2 Approximate (flipped) Walsh structure occurs exactly at the chance rate.** Census of all tight
4-sets (angular span `<= C N^(-1/4)`) on 3 circles per `k` (`tight_quadruple_census.py`): the fraction
that are Walsh 2-frames with flips (odd quartic support `<= 4` primes) versus the random-sign
prediction `P(Bin(k,1/2) <= 4)`:

| k | C=64 | C=256 | binomial |
|---|---|---|---|
| 8 | 0.640 | 0.634 | 0.637 |
| 10 | 0.394 | 0.385 | 0.377 |
| 12 | 0.190 | 0.194 | 0.194 |
| 13 | 0.117 | 0.136 | 0.133 |
| 14 | 0.067 | 0.086 | 0.090 |
| 15 | 0.076 | 0.055 | 0.059 |
| 16 | 0.013 | 0.046 | 0.038 |

There is no enrichment of approximate Walsh structure: tight quadruples are quasi-random in the sense
relevant to aligned flats.

**9.3 Exact structure: multiplicative rectangles, from a Bohr-set mechanism.** Exact parallelograms
(pure Hadamard-4 cores) are enriched over the chance rate `2^-k`, by a factor growing from about 2
(`k = 8`) to about 400 (`k = 15`), and almost all of them are *integer rectangles*
`z_1 z_2 = unit * z_3 z_4` (verified exactly in Gaussian integers for every one found): e.g. `k = 13`,
`C = 64`: 15 of 1,461 tight 4-sets; `k = 16`, `C = 256`: 6 of 3,022. In windows of `M = 16` points:
extremal windows contain on average 6.1 parallelograms (5.9 rectangles), typical windows of the same
size 2.7 (2.6), random subsets of the sign cube 0.16 (0.02); affine 3-cubes (pure Hadamard-8 cores) are
almost absent (0.02 per extremal 16-window, 0.001 per typical one) (`window_structure_out.txt`).

Mechanism (heuristic): if `x_1, x_2` lie in a window, swapping a subset `J` of their differing
coordinates gives `x_3, x_4` with `x_1 + x_2 = x_3 + x_4` and `theta_3 + theta_4 = theta_1 + theta_2`
exactly, so one angle condition is saved. This predicts that the rectangle fraction among tight
4-sets is about `2^(-k/2)/lambda`, `lambda = (2/pi) 2^k C N^(-1/4)`:

| (k, C) | (8,64) | (10,64) | (12,64) | (13,64) | (14,64) | (15,64) | (10,256) | (12,256) | (13,256) | (14,256) | (16,256) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| predicted | .0063 | .0089 | .0083 | .0105 | .0139 | .0156 | .0042 | .0020 | .0021 | .0035 | .0069 |
| observed | .0074 | .0065 | .0056 | .0103 | .0111 | .0120 | .0058 | .0029 | .0026 | .0031 | .0020 |

(agreement within a factor 1.5 except the sparsest case, 6 rectangles). The absolute number of tight
rectangles per circle, of order `lambda^2 2^(k/2)`, tends to zero at fixed `C` because `N^(1/4)` grows
faster than any `2^(ck)` (`N >= exp(c k log k)` for products of `k` distinct primes). Rectangles are
2-dimensional pure cores (one matching absent); pure cores of order `>= 28` are excluded by Theorem B.

**9.4 Known extremal families.** The eight-point translated Pell family (item 511) is an affine
3-cube in sign space (odd-weight words of `F_2^4`) with four heavy coordinate labels plus one light
non-affine class (`F = 1+2i`, majority pattern): a pure Sylvester core of order 8 with one light
augmentation. Checked exactly at `n = 1` and `n = 61` (`check_pell_eight_structure.py`: eight distinct
equal-norm primitive points, `log10 N = 24` and `622`, largest exact chord `10^4.34` and `10^7.65`
times `sqrt R`, block signs an affine 3-cube, the `F` class not affine). It is consistent with
Theorem B (order `< 28`) and illustrates that the known structured extremal configurations are
*pure-core-like*, not flipped-frame-like.

**Verdict on plausibility.** The observable data give no support to EH: approximate Walsh structure
appears only at the chance rate, and the only enriched exact structure (rectangles, small pure cores)
is of the kind already excluded at larger orders by Theorem B. If large clusters existed, the data
and the angle-condition heuristic suggest they would look quasi-random, i.e. like the Paley barrier.
An extraction into Walsh frames is therefore not a plausible intermediate step; a successful route
must treat quasi-random fully flipped profiles directly.

## 10. Missing lemma and why the method stalls

**Missing lemma (precise).** *Fully flipped Hadamard phase lemma.* For every `C` and `K` there is
`M_0` such that no endpoint cluster with `M >= M_0` points and `W <= K M log M` has a fully flipped
Hadamard-core profile (any normalised Hadamard matrix of order `M`, `b >= 5` copies, one flip per row
at distinct primes, capacity `<= B`). Equivalently, by Lemma 7.2, the system

```text
Phi_a^unflipped  in  (pi/(2M)) m_a + (2/M) sum_x H(x,a) sigma_x psi_x - sum_{a(x)=a} psi_x + [-Delta/2, Delta/2]
```

has no solution in actual Gaussian blocks. Together with Theorems B'' and C this would settle every
Hadamard-core profile with at most one flip per row; a general growth improvement would further need
an extraction of Hadamard-core structure from arbitrary clusters, which Section 9 does not support.

**Why the method stalls.**
* Every rational-character consequence (Theorem F for all `m`, including all residue pigeonholes)
  is satisfied by the barrier, because flips make the rational saturation half-integral (Lemma F.1)
  and item 393 controls all half-integral characters.
* Aligned-flat extraction needs small L1-L2 defect `Delta_c` (Theorem A); quasi-random cores have
  `Delta_c >= c_0 W` for all characters (items 334, 336, 391, 393).
* Small-prime collisions contribute only the aggregate `M^2 log M`; per-pair forms are equivalent in
  strength to smallness of the reduced chords, i.e. to the arc condition itself.
* What remains is `M-1` simultaneous Diophantine approximations of Gaussian blocks to targets moving
  with `M` free prime angles; single-character (half-height) consequences are all satisfied, and no
  available tool (Liouville, Roth, Baker/Matveev, Subspace with uniform exceptional sets) controls the
  joint system uniformly in `M`.

## 11. Files (all under `growth/walsh_extraction_checks/`; run from that directory)

* `check_defect_identity.py` -- Theorem A, exact (Fractions) and literal Gaussian heights. PASS.
* `check_hadamard_residue_pigeonhole.py` -- Lemma B.1 on 280 random weightings (Walsh 16/32/64,
  Paley 12/20/44/68, full and row-deleted); literal pure cores versus (B.1); Theorem B/B'' constants
  in exact integers; thresholds `M_1(C)` from the exact inert-prime function `F(M)`. PASS.
* `check_theorem_C_ingredients.py` -- density lemma, Theorem B'' certificates, saturation. PASS.
* `paley_collision_check.py`, `paley_collision_out.txt` -- per-pair collision data of Section 7.2.
* `cluster_search.py`, `cluster_search_large.py`, `analyze_clusters.py`, `tight_quadruple_census.py`,
  `window_structure.py`, outputs `census_C64.txt`, `census_C256.txt`, `window_structure_out.txt`,
  `*clusters*.json` -- Section 9.
* `check_pell_eight_structure.py` -- Section 9.4.
* `gauss_common.py` -- exact Gaussian integer helpers.

The repository `/home/user/Jarnik` was not modified.
