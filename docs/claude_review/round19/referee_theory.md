# Round 19 referee report on theory.md (abc lane, unit pencils)

Status: complete. Referee of `round19/theory.md`. Outcome: **survives** (no fatal or major defect; 7 minor/cosmetic items).
Labels: *proved* = re-derived here; *exact* = independent exact computation; *evidence* = finite sample, used in no proof.

## 0. Verdict

**The result survives.** Theorem 2.4, Corollary 2.5, Lemmas 3.1–3.2, Theorem 3.3 and Corollary 3.4 are correct as
proved. Positivity (F3) is the only deep input, and it was refereed earlier (`L8.md` Thm A(a), Keel–McKernan). The
conclusion is that **no unit pencil `lambda A/C` on `Mbar_{0,8}`, of any degree and for any `lambda ∈ C^*`, has
positive abc-gain**, and neither does any composite `h(A/C)`. So round-18 Prop 2.6 never applies.

* *Object.* The gain is computed for exactly round-18's object. The `X`-level fibres reproduce round-18's `F_red`
  term by term, for every `lambda` (`ord_{D_T} F_1 = nu_T(B) - mu_T`). Exceptional divisors of the base-locus blow-up
  have `beta~`-degree 0 under round-18's convention and are not counted by Prop 2.6's proof (§1).
* *Proof.* I re-derived every step (§2). The log-cotangent filtration, its quotients `pi_k^* psi_k`, the `df`-divisor
  formula, the line-subsheaf bound and the Boolean-cube count are all right. The margin really is `32 - 31 = 1`.
* *Computation.* My code is independent and exact:
  * every class number, plus the class identity `Lambda_k = pi_k^* psi_k` checked against all 1701 F-curves;
  * `N(u) >= 32` on 1.43 million vectors;
  * every claimed extreme case, and 22 Belyi-type special `lambda`, each attaining `-32` with equality in Cor 2.5;
  * the author's `D <= 3`, 14-`lambda` evidence family in full, matching histogram entry by histogram entry;
  * random `D <= 5` unit pencils at all candidate special `lambda`;
  * a stress test of Theorem 2.4 on non-unit pencils.

  All agree, and no inequality is violated (§3).
* *Defects.* None is fatal or major. There are seven minor or cosmetic items (§5). Two of them change how the
  headline should read:
  * **m1**: `k` must be `k(f)`, the least number of points through which `f` factors, including the anchor when needed;
  * **m2/m3**: the closed lane is the **unit-pencil** lane of Prop 2.6. Round 18's wider gap ("abc on auxiliary
    pencils", which includes non-unit pencils with boundary base locus) stays open, as theory.md §4 itself says. The
    stop rule's literal threshold `-32` is proved only when `k(f) = 4`.
* Nothing here touches `(L_8)`, the uniform bound or the growth rate, and the write-up does not claim otherwise.

## 1. Object check: is the gain the right object?

Read in full: `round19/theory.md`; `round18/README.md`; `round18/theory.md` §2.4 (Prop 2.6, identity, partial
closure); `round18/referee_theory.md` §2.6 and §3.4 (items B, F); `threshold/L8.md` Thm A(a) and the referee's item 3.

**The gain of theory.md §1 is exactly round-18's gain, for every rational (indeed complex) `lambda`.** Re-derived:

* *Fibre over 1, boundary part.* `f - 1 = (lambda A - C)/C = -B/C`, so on `X`
  `ord_{D_T} F_1 = ord_{D_T}(f - 1) + ord_{D_T} F_inf = (nu_T(B) - nu_T(C)) + max(0, nu_T(C) - nu_T(A)) = nu_T(B) - mu_T`.
  Likewise `ord_{D_T} F_0 = nu_T(A) - mu_T` and `ord_{D_T} F_inf = nu_T(C) - mu_T`. So `D_T` lies in `F_0 + F_1 + F_inf`
  iff `max(nu_T(A), nu_T(B), nu_T(C)) > mu_T`, which is round-18's boundary count (including pairs, `mu = 0`). The three
  fibres share no component, so no `D_T` is counted twice.
* *Interior part.* `A`, `C` are units, so the only interior components are the strict transforms `Y~` of the
  non-difference irreducible factors `Y` of `B`, with `Y~.beta = 31 deg Y - sum_{|T|>=3} nu_T(Y)` (`H.beta = 31`,
  `E_T.beta = 1`). A `Q`-irreducible squarefree factor splits over `C` into distinct conjugates whose `beta`-degrees add
  up to the same number, so working over `C` changes nothing.
* *Fibre class.* `L.beta = F_0.beta = sum_{e(T)>0} e(T) = 31 D - sum_T mu_T` (the total transform of a product of `D`
  differences has `beta`-degree `31 D`, since each `x_i - x_j` vanishes on the 31 divisors `D_T`, `T ⊇ {i,j}`).
* *Exceptional divisors of the base-locus blow-up.* The base locus of `lambda u` on `X` is `F_0 ∩ F_inf`, a union of
  codimension-2 boundary strata `D_T ∩ D_T'` (`e(T) > 0 > e(T')`). Round 18 uses `beta~ = pi^* beta` (functional
  `E ↦ pi_*E . beta`), which gives every exceptional divisor `beta~`-degree 0, and Prop 2.6's proof never counts them
  (on `Bal^pc` no prime divides two `d_T`). theory.md works on `X` and drops them; that is the same number. ✓
* *All `lambda`.* Round-18 item F: abc applies to `qC - pA = qB` for `lambda = p/q`, so Prop 2.6 covers every
  `lambda ∈ Q^*`. theory.md proves the bound for `lambda ∈ C^*`, a superset. ✓
* *Units.* `O^*(M_{0,8})/C^*` = degree-0 products of the 21 differences (complement of a hyperplane arrangement), rank 20;
  `ord_{D_T}(u) = e(T)` because the order of a pulled-back function along `E_T` is the multiplicity at the generic point
  of `L_T`, where the iterated Kapranov blow-up is a single smooth blow-up. ✓ So "unit pencil" in theory.md is exactly
  round-18's `lambda A/C`, of any degree `D` (the degree is invisible: only `e = e(A) - e(C)` and `lambda` matter).

**Verdict on the object:** correct. One caveat on the *scope* of the word "lane" (issue m2, §5).

## 2. Step-by-step re-derivation

**Lemma 2.1 (filtration) ✓.** The forgetful map is the universal curve, hence log smooth for the divisorial log
structures, so `0 -> pi^* Omega_Y(log) -> Omega_{Y'}(log) -> omega_pi(sum sigma_i) -> 0` is exact and locally split;
the three local charts in the proof are the standard ones, and the boundary of `Mbar_{0,m+1}` is indeed `pi^{-1}∂`
plus the `m` sections (a `D_S` with `S^c = {i, m+1}` is a section, all others are components of pull-backs).
Pulling back preserves exactness and local splitting, so each `F_{k-1} ⊂ F_k` is a subbundle with invertible
quotient. Degrees: I recomputed both routes independently and, more strongly, checked the **class** identity
`pi_k^*(K_k + ∂_k) - pi_{k-1}^*(K_{k-1} + ∂_{k-1}) = pi_k^* psi_k` against all 1701 F-curves of `Mbar_{0,8}`
(which span `N_1`), for `k = 4..8`, with 0 mismatches (`r1_classes.py`). `Lambda_k.beta = 16, 24, 28, 30, 31`. ✓
Also `min_F Lambda_k.F = 0` (each `Lambda_k` is nef, as a pull-back of `psi`; not needed).

**Lemma 2.2 ✓.** Minimality of `k` makes `O(Delta) -> Lambda_k` non-zero (the kernel of `F_k -> Lambda_k` is
`F_{k-1}`), so `Lambda_k - Delta ~ E >= 0`. Only `E.beta >= 0` is needed. (F3) is `L8.md` Thm A(a), which I re-read:
symmetrize a non-zero effective `E` over `S_8`, use Keel–McKernan Thm 1.3 (effective cone of `Mbar_{0,n}/S_n`
simplicial on `B_2..B_{floor(n/2)}`), and `S_8`-invariance of `beta`. It was refereed (`threshold/referee_L8.md`
item 3). ✓

**Lemma 2.3 ✓.** Re-derived each of the four local cases with the log frame (`ds/s` along boundary components):
orders `m - 1` (interior, finite fibre), `m` (boundary, finite fibre), `-m` (boundary, polar), `-(m+1)` (interior,
polar), `0` otherwise (the `dy`-part of `df` restricts to `d(f|_Y) ≠ 0`). Summing gives (2.1). The map
`g ↦ g·theta`, `O(div theta) -> Omega^1(log ∂)`, is regular in codimension one, hence regular (`X` normal, target
locally free), and injective. Sign sanity check on `Mbar_{0,4} = P^1`: `div(dz) = [0] + [1] - [inf]`, degree
`1 = deg Omega^1_{P^1}(log 3 pts)`. ✓

**Theorem 2.4 ✓.** `f = g ∘ pi_S` gives `theta = pi_S^* dg`, a rational section of `F_{|S|}` (with `S` ordered first);
since `F_{|S|}` is saturated, `O(div theta)` lands in it, and Lemma 2.2 applies with `k' <= k(f)`. All fibres
`F_p ~ L` (`div(f - c) = F_c - F_inf`), so `(div theta).beta = sum_{c ∈ C}(L.beta - R_c^int) - L.beta - R_inf^int`,
and the rearrangement is correct. Termwise `L.beta - R_p^int >= n_p` uses `Y.beta >= 0` for every prime `Y` (F3). ✓
Equality cases re-checked by hand: cross-ratio `3·16 = 2·16 + 16`; `4z(1-z)`: `32 + 16 + 32 = 2·32 + 16`
(fibre over 1 is `2·chi^*(1/2)`, interior, `R^int = 16`; all other fibres reduced interior). ✓

**Corollary 2.5 ✓.** `gain = sum_{0,1,inf}(L.beta - R_p^int) - 2 L.beta - n(f)` and Theorem 2.4. A `D_T` on which `f` is
constant `c` lies in `F_c` (and conversely), so `n'(f) = sum_{p ∈ P^1} n_p`, giving `gain <= Lambda_{k(f)}.beta - n'(f)`. ✓

**Lemma 3.1 (Schwartz–Zippel on the cube) ✓.** All four cases of the induction are right; the bound may be fractional
for `d > n - 1`, which is harmless (a non-zero multilinear polynomial is determined by its values on the cube, so it
has `>= 1` non-zero there). This is the classical Reed–Muller minimum-distance bound, valid over any field.

**Lemma 3.2 ✓.** `q(x) = sum e_ij x_i x_j` is non-zero iff `e ≠ 0` iff `u` non-constant (distinct irreducible
differences). `q(1_T) = 0` for `|T| <= 1` and for `T = [7]` (`sum e = 0`), so the `>= 32` non-zeros are among the 119
divisors. Cross-ratio: `q = x_1(x_2 - x_3)`, 32 non-zeros. ✓ (Independent check `r2_support.py`, §3.)

**Theorem 3.3 ✓, with a definitional caveat (m1).** `F_0, F_inf ⊂ ∂` with `n_0 + n_inf = N(u)`; Corollary 2.5 and
Lemma 3.2 give `gain <= Lambda_{k(f)}.beta - 32 - n_1 = -2^{8-k(f)} - n_1`. The proof is valid **with `k := k(f)` of
§2.3** (least `|S|` with `f = g ∘ pi_S`). The statement's "number of the 8 points the unit involves" must be read that
way, *including the anchor 8 when it is needed*: for `u = prod (x_i - x_j)^{e_ij}` in the anchor-8 chart,
`k(f) = |I| + [some vertex sum sum_j e_ij ≠ 0]`, `I` the indices touched by `e` (`u` is independent of `x_8 = inf`
iff it is `PGL_2`-invariant iff all vertex sums vanish, by `x ↦ 1/x`). Counting only chart indices gives `k = 3` for
`(x_1 - x_2)/(x_1 - x_3)`, where `Lambda_3` is undefined and `-2^{8-3} = -32` would be an unproved claim. Since
`Lambda_k.beta` increases with `k`, any `k >= k(f)` gives a valid bound, so the sign conclusion `gain <= -1` is
unaffected.

**Corollary 3.4 ✓.** `f^*(p) = u^*(h^*p)` in codimension one; `R_p(f) = sum_{t ∈ h^{-1}(p)} R_t(u)` (distinct `t` give
disjoint supports); `|S| >= deg h + 2` by Riemann–Hurwitz for `h`; Theorem 2.4 for `u` with the `t = 0, inf` terms equal to
`L_u.beta`; and `L_u.beta = sum_{e(T)>0} e(T) >= N_0`, `= sum_{e(T)<0} |e(T)| >= N_inf`. The cases `j = 0, 1, 2` all
give `gain(f) <= Lambda_k.beta - N(u)`. ✓

**§4.2 necessary conditions ✓** (1: `gain >= 1` and Cor 2.5; 2: if two fibres lie in `∂`, a Möbius change of the
target makes `f` a composite `h ∘ u`; 3: `k(f) = 4` means `f = h ∘ chi`; 4: `gain = E_0 + E_1 + E_inf - 2 L.beta` and
`E_q < L.beta` because `R_q > 0`; 5: Cor 2.5). The Campana remark is not used; it is imprecise (m3).

## 3. Independent exact recomputation

All scripts are in `referee_theory_checks/`, written from scratch; none imports the author's code. Arithmetic is exact
(Python integers, `Fraction`, python-flint `fmpz_mpoly` factorization and composition over `Z`). Every job ran under
`prlimit --as=3500000000`, at most two at a time, with jsonl output per orbit/chunk and resumable drivers.

### 3.1 Classes — `r1_classes.py` -> `r1_classes.out`

* 119 boundary divisors; `beta` satisfies the Keel relations (all 70 quadruples x 3 pairings give 16). ✓
* `K.beta = 10`; `(K + ∂).beta_k = 1, 5, 17, 49, 129` = `(k-4)2^{k-3} + 1`; pulled back: `16, 40, 68, 98, 129`;
  `Lambda_k.beta = 16, 24, 28, 30, 31`. ✓
* **Stronger than the author's check:** `Lambda_k = pi_k^* psi_k` holds as an identity of classes (0 mismatches against
  all 1701 F-curves, which span `N_1`), and is independent of the choice `(a, b)` in the `psi` formula. Sanity of the
  F-curve pairing: `K.F = 2, 1, 0, 0, -1` on the types `2222, 1223, 1133, 1124, 1115`, as in round 18. ✓
* `psi_j.beta = 31` for every `j` and every `(a, b)`. ✓

### 3.2 Boolean support count — `r2_support.py` -> `r2_support.out`

* Exhaustive over all 1,121,610 vectors `e ∈ {-1,0,1}^21` with `sum e = 0` and support `<= 6`: `min N = 32`.
* Exhaustive over 8,820 vectors with entries in `{±1, ±2}`, support `<= 3`: `min N = 32`.
* 300,000 random vectors (entries `-3..3`, any support): `min N = 32`.
* For every vector, `q(1_T) = 0` off the 119 divisors (asserted). Cross-ratio: `N = 32`, `k(f) = 4`. ✓
* *Evidence:* smallest `N` found per `k(f)` (my definition, m1): `32, 48, 40, 48, 44` for `k = 4..8` (author: `..., 50`).

### 3.3 Gains: every claimed extreme case — `r3_lib.py`, `r3_extreme.py` -> `r3_extreme.out`

`r3_lib.gain_data` computes, for arbitrary coprime `A`, `C` (units or not) and `lambda = p/q`, the exact factorization
of `A`, `C` and `B = qC - pA` over `Z`, the exact orders `nu_T` of every factor along all 98 centres `L_T`
(substitution `x_i -> x_{i0} + t x_i`, lowest `t`-power) and all 21 pairs (exact division), and from these the fibres
`F_0, F_1, F_inf` **on `X`**. Built-in consistency checks, asserted on every pencil: (a) `F_0.beta = F_1.beta =
F_inf.beta` (three independent computations of `L.beta`); (b) the three fibres share no `D_T`; (c) the multiplicity form
`sum ((m_Y - 3)/3) Y.beta` equals `L.beta - F_red.beta`; (d) `nu_T(B) >= mu_T`.

| pencil | `lambda` | `k(f)` | `L.beta` | gain | `n(f)` | `N(u)` | `n_1` | Thm 3.3 bound `Lambda_k - N - n_1` | `gain + n <= Lambda_k` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `z^D`, `z = (x_1-x_2)/(x_1-x_3)`, `D = 1..5` | 1 | 4 | `16D` | -32 | 48 | 32 | 16 | -32 | 16 <= 16 (equality) |
| `4z(1-z)` | 4 | 4 | 32 | -32 | 48 | 48 | 0 | -32 | equality |
| round-18 referee example `(x_1-x_2)^2/((x_1-x_3)(x_2-x_3))` | -1/4 | 4 | 32 | -32 | 48 | 48 | 0 | -32 | equality |
| `lambda z^a (1-z)^b`, `a + b <= 5`, critical value 1 (10 cases, e.g. `27/4`, `256/27`, `3125/108`) | | 4 | `16(a+b)` | -32 | 48 | 48 | 0 | -32 | equality |
| `lambda z^a/(1-z)^b`, `5 >= a > b`, critical value 1 (10 cases, e.g. `-4/27`, `1/16`, `108/3125`) | | 4 | `16a` | -32 | 48 | 48 | 0 | -32 | equality |
| Ptolemy | 1 | 4 | 16 | -32 | 48 | 32 | 16 | -32 | equality |
| Ptolemy | -1 | 4 | 16 | -32 | 32 | 32 | 0 | -16 | 0 <= 16 |
| `(x_1-x_2)/(x_3-x_4)` | 1 | 5 | 24 | -48 | 48 | 48 | 0 | -24 | 0 <= 24 |
| 8-point unit `(x_1-x_2)(x_3-x_4)(x_5-x_6)/((x_1-x_3)(x_2-x_5)(x_4-x_7))` | 1 | 8 | 30 | -58 | 59 | 58 | 1 | -28 | 1 <= 31 |
| star unit (min-`N` at `k = 8`) | 1 | 8 | 30 | -44 | 45 | 44 | 1 | -14 | 1 <= 31 |

Every value agrees with theory.md (cross-ratio `-32`, `n = 48`, equality in (2.2); `4z(1-z)` `-32`) and with round 18
(Ptolemy `-32`; `(x_1-x_2)/(x_3-x_4)` `-48`; referee's `-1/4` example `-32`). Theorem 3.3's bound
`Lambda_k - N(u) - n_1` is **attained** by every `k = 4` composite of the cross-ratio, so it cannot be improved
in general for `k = 4`; for `k >= 5` it is far from sharp in the examples.

### 3.4 Gains: author's evidence family, reproduced — `r3_enum.py` -> `r3_enum_D{1,2,3}.jsonl`, `r3_compare.out`

The family is: `A` ranges over one representative per `S_7`-orbit of multisets of `D` pairs (orbits computed by brute force
over all 5040 permutations); `C` ranges over all multisets of `D` pairs sharing no pair with `A`; and `lambda` ranges
over the author's 14 values. I ran one process per orbit.

| `D` | orbits | triples | max gain | `hist(gain + n)` (referee = author) | max `gain + n - Lambda_{k(f)}.beta` |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | 280 | -32 | 16:10, 0:270 | 0 |
| 2 | 3 | 8,260 | -32 | 16:18, 8:28, 4:54, 2:24, 0:8136 | 0 |
| 3 | 8 | 138,600 | -32 | 16:10, 8:52, 4:148, 2:270, 1:193, 0:137927 | 0 |

**The histograms agree with the author's `c3_D*.jsonl` entry by entry**, and so do the triple counts and the maximum
gain. One more observation, my own and in the author's favour: Corollary 2.5 with the *refined* constant
`Lambda_{k(f)}.beta` (not just 31) holds on all 147,140 triples and is attained with equality.

### 3.5 Random unit pencils with all candidate special `lambda`, `D <= 5` — `r3_sample.py` -> `r3_sample.jsonl`, `r3_summary.out`

**Sample.** 6,000 random pairs `(A, C)`, 200 chunks with seeds `1000 + chunk`, weighted towards `D = 4, 5`, and
supported on 3 to 7 random points. For each pair I evaluated `lambda = ±1` and **every candidate special `lambda`**
(23,105 pencils in all). A candidate is either:
* (i) a ratio at which the initial forms of `A` and `C` along some `D_T` with `e(T) = 0` (pairs and centres `L_T`)
  are proportional; or
* (ii) a rational root of `disc_s(C - lambda A)` on either of two random lines (lines on which the degree drops are
  redrawn).

**Why these candidates suffice.** For `lambda` outside (i) ∪ (ii), `B` is reduced and has no extra boundary order, so
`gain = -N(u)`. I also asserted this on a generic `lambda = 1000003/7919` for every pair.

**Results.** All checks passed (§3.3 (a)–(d), Thm 3.3's bound, `gain + n <= Lambda_{k(f)}.beta`, `N >= 32`).

| `D` | pairs | pencils evaluated | max gain | max `gain + n - Lambda_{k(f)}.beta` | min `N` |
| --- | --- | --- | --- | --- | --- |
| 1 | 757 | 2,969 | -32 | 0 | 32 |
| 2 | 740 | 3,226 | -32 | 0 | 32 |
| 3 | 735 | 2,999 | -32 | 0 | 32 |
| 4 | 1,483 | 5,947 | -32 | 0 | 48 |
| 5 | 2,285 | 8,964 | -32 | 0 | 48 |

This is *evidence*, not an exhaustive `D = 5` run. Candidate set (ii) is complete only with probability 1, and the
sample is random. It is consistent with the theorem, and the `D = 5` maximum is the round-18 value `-32`. No proof
uses it, and the theorem makes it unnecessary.

### 3.6 Stress test of Theorem 2.4 on non-unit pencils — `r4_nonunit.py` -> `r4_nonunit.jsonl`

Theorem 2.4 is claimed for **every** rational function, and for non-unit pencils it is the only bound. I tested
`gain + n(f) <= Lambda_{kb}.beta` on non-unit pencils. Here `kb = min(8, #indices + 1) >= k(f)`, which gives a
weaker inequality than the theorem's.

**Sample.** 5,000 random coprime pairs `(A, C)`, `D <= 5`, built as products of powers (up to 3) of:
* differences;
* interior linear forms `a x_i + b x_j - (a+b) x_k` and `x_i + x_j - x_k - x_l`;
* interior quadrics `(x_i-x_j)(x_k-x_l) + c (x_i-x_k)(x_j-x_l)`.

Their base loci can meet the interior. Each pair was evaluated at `lambda = ±1` and at all candidate special
`lambda`: 23,454 pencils in all.

**Results.**
* All consistency checks hold. In particular `F_0.beta = F_1.beta = F_inf.beta` held on every pencil, which tests the
  `nu_T` computations on interior factors.
* The maximum of `gain + n - Lambda_{kb}.beta` is **0** (attained), so there are **no violations**.
* The maximum gain is `-32`.

This is evidence that the general inequality is right, beyond the unit case.

## 4. Completeness / evidence vs proof

* **No enumeration and no list of special `lambda` is used in any proof.** Theorem 3.3 covers every `lambda ∈ C^*`
  at once, through Corollary 2.5 (which holds for every pencil) and Lemma 3.2 (which depends only on `e`, not on
  `lambda`). So the completeness questions that dominated round 18 (item F: which `lambda` are special; which `D`)
  do not arise. I checked that `lambda` enters nowhere in the proof except through `n_1 >= 0`, which is dropped. ✓
* **The 14-value list in `c3_gain_vs_n.py` is not a complete list of special `lambda`** (it omits, e.g.,
  `27/4` for `z(1-z)^2`, `256/27`, `-4/27`, `3125/256`, all of which I evaluated in §3). It is labelled evidence and is
  used only as a consistency test, so this is not a defect of the proof; but the sentence in theory.md §6 "It extends
  the referee's `D <= 3` special-`lambda` result to the 14 listed values of `lambda`" should not be read as covering
  all special `lambda` at `D <= 3` (it is the round-18 referee's run that did that, with probability-1 completeness).
* **Evidence labels are correct** apart from one number: theory.md §6 / `c4_support.out` report "smallest `N` found
  for units involving exactly `k = 4..8` points is 32, 48, 40, 48, 50". My search (`r2_support.py`) finds `N = 44` at
  `k(f) = 8` (`u = (x_1-x_5)(x_1-x_6)(x_1-x_7)/((x_1-x_2)(x_1-x_3)(x_1-x_4))`, confirmed by the exact gain code:
  `N = 44`, gain `-44`). Either the author's `k` is defined differently (see m1) or the search missed it. It is
  evidence-only and is not used, so this is cosmetic (m5).
* **Proof vs citation.** The proof rests on three cited facts: (F1)-(F2) (Knudsen, Keel; standard), (F3) (Keel–McKernan
  via the refereed `L8.md` Thm A(a)) and (F4) (Keel–McKernan/Pandharipande formula for `K`). The identification of
  the quotient `Lambda_k` with `pi_k^* psi_k` is not needed for the proof (the determinant route suffices) and I
  verified it as a class identity anyway. Nothing depends on Bruno–Mella (round 18's "read from the abstract only"
  citation is no longer needed: Theorem 3.3 and Cor 3.4 subsume the morphism case `h ∘ chi`).

## 5. Issues

No fatal or major defect was found. Every proof step was re-derived (§2) and every number recomputed (§3).

| item | where | severity | defect | fix |
| --- | --- | --- | --- | --- |
| m1 | Thm 3.3, §0, §7(1), abstract | minor | `k` = "the number of the 8 points the unit involves" is ambiguous about the anchor. Counting chart indices gives `k = 3` for `(x_1-x_2)/(x_1-x_3)`, which is undefined, and `k = 7` instead of 8 for units that also involve `x_8 = inf`. Read that way, the bound would be twice as strong as what is proved. | Define `k := k(f)` as in §2.3. In the anchor-8 chart, `k(f) = |I| + [some vertex sum sum_j e_ij ≠ 0]`. The sign conclusion `gain <= -1` is unaffected. |
| m2 | §0 first line, abstract ("the abc lane ... is closed by theory") | minor (interpretive) | Round 18 defined the lane more widely than Prop 2.6. Its theory §2.4 (after Prop 2.6, and the closing sentence "open only through unit pencils of degree `>= 5` (or non-unit pencils with boundary base locus)") and README §4 ("abc applied to auxiliary pencils") both include non-unit pencils whose base locus lies in codimension-2 boundary strata. theory.md §4 correctly leaves these open, and its Theorem 2.4 allows them `gain <= 31`. | Say: "the **unit-pencil** abc lane (round-18 Prop 2.6, the only proved part of the lane) is closed, together with all composites `h ∘ u`. The barrier-map gap 'abc on auxiliary pencils' stays open for non-unit pencils in `P°`, and there it needs an unproved gcd extension as well." |
| m3 | §1 / §0 "stop rule fires in a stronger form" | minor | The stop rule's literal threshold (`gain <= -32` at `D = 5`, and at `D = 4` with `lambda ≠ ±1`) is **not** proved. What is proved is `gain <= -2^{8-k(f)}`, and `<= -32` only when `k(f) = 4`, where `lambda u = h ∘ chi` and Riemann–Hurwitz applies. The stop rule exists to decide whether Prop 2.6 can ever apply, and for that the sign is enough. | Say "the stop rule's purpose is met by proof (no positive gain, any `D`, any `lambda`). Its numerical threshold `-32` is proved only for `k(f) = 4`, and for `k(f) >= 5` it is evidence." |
| m4 | §4.2 Campana remark | minor (not used) | The orbifold base has to be computed on the resolved pencil `X~ -> P^1`, where exceptional components of the base-locus blow-up can lower the gcd multiplicity of a fibre. The claim "at most two fibres of gcd-multiplicity `> 1`" holds for `X~`-fibres, not for fibres on `X`. | Add "on `X~`". |
| m5 | §6 (`c4_support.out`, evidence) | cosmetic | The smallest `N` found at `k = 8` is 44, not 50 (`u = (x_1-x_5)(x_1-x_6)(x_1-x_7)/((x_1-x_2)(x_1-x_3)(x_1-x_4))`). | Correct the evidence line, or state the author's `k`. |
| m6 | §6 last sentence on `c3` | cosmetic | The 14 fixed values of `lambda` are not the special values. They miss, for example, `27/4`, `256/27` and `-4/27`, so the run does not "extend the referee's special-`lambda` result". | Reword as a consistency test at 14 fixed `lambda`. |
| m7 | §7(4) | cosmetic (not used) | "Holds verbatim for any movable class `gamma` with 31 replaced by `max_k Lambda_k.gamma`". For a class `gamma` that is not `S_8`-invariant, `Lambda_k.gamma` depends on the ordering, and the ordering must put `S` first. "Movable" must mean dual to the pseudo-effective cone (BDPP). | State the bound as the max over `k` and over orderings, for `gamma` non-negative on effective divisors. |

## 6. Corrected statement

**Theorem (round 19, refereed; proved).** Let `X = Mbar_{0,8}` and let `beta` be the balanced class (`D_T.beta = 1`).

1. *(beta-weighted log Riemann–Hurwitz.)* Let `f` be any non-constant rational function on `X`, and let
   `k(f) = min{|S| : f = g ∘ pi_S}`. Then
   `sum_{p ∈ P^1} (L.beta - R_p^int) <= 2 L.beta + Lambda_{k(f)}.beta`, where `Lambda_k.beta = 32 - 2^{8-k}`
   (`= 16, 24, 28, 30, 31`). Equality holds for the cross-ratio and for `4z(1-z)`. Hence
   `gain(f) <= Lambda_{k(f)}.beta - n'(f) <= 31 - n'(f)`.
2. *(Unit pencils.)* For every non-constant unit `u` of `M_{0,8}` and every `lambda ∈ C^*`, in particular for every unit
   pencil `lambda A/C` of round-18 Prop 2.6, of any degree `D` and any `lambda ∈ Q^*`:
   `gain(lambda u) <= Lambda_{k(f)}.beta - N(u) - n_1 <= -2^{8 - k(f)} <= -1`, with `N(u) >= 32`. For `k(f) = 4`,
   `gain <= -32`. The same bound `gain(h ∘ u) <= Lambda_{k(u)}.beta - N(u) <= -1` holds for every composite with a
   non-constant `h ∈ C(z)`.
3. *(Consequence.)* Round-18 Prop 2.6 never applies. **The unit-pencil abc lane is closed by proof**, and the round-18
   stop rule's purpose is met for all `D` and all `lambda`. Its literal threshold `-32` is proved only for
   `k(f) = 4`, and is evidence for `k(f) >= 5`. The nested-chain extension (README §3.3, J2) is moot for unit
   pencils.
4. *(Not closed.)* An abc route through **non-unit** pencils with base locus in codimension-2 boundary strata would
   need two things: (i) a gcd extension of Prop 2.6, which is unproved, and (ii) a pencil with
   `n'(f) <= 31 - 2^{8-k(f)} <= 30`, three non-reduced fibres with `sum 1/mu_p < 1`, `k(f) >= 5`, not composite with a unit.
   None is known, and Theorem 2.4 allows such pencils `gain <= 31`.
5. *Unchanged:* `(L_8)` is open. The uniform bound is conditional on Vojta. The best growth bound is
   `M <= (2/3 + eps) log R/loglog R`.

## 7. Not re-verified (not load-bearing)

* `theory_checks/c2_psi_cone.py` (`beta` is not in the cone of `psi`-monomials) and `t2.py` (`(K + ∂)^4` pairings).
  They only explain why (F3) is taken from Keel–McKernan rather than from explicit curves, and no proof uses them.
* `mbar.py` / `t_mbar.py` (author's intersection engine), which only those two scripts use.

## 8. Files (`referee_theory_checks/`)

* `r1_classes.py` -> `r1_classes.out`: Keel relations, `K.beta`, `(K + ∂).beta_k`, `Lambda_k.beta`, class identity
  `Lambda_k = pi_k^* psi_k` on 1701 F-curves, `psi_j.beta`.
* `r2_support.py` -> `r2_support.out`: `N(u) >= 32` (exhaustive small support plus random), min `N` per `k(f)`.
* `r3_lib.py`: exact fibre/gain engine for arbitrary pencils (factorization over `Z`, exact `nu_T`, three-way `L.beta`
  check, multiplicity form, candidate special `lambda`).
* `r3_extreme.py` -> `r3_extreme.out`: every claimed extreme case and 22 Belyi-type special `lambda`.
* `r3_enum.py`, `run_enum_D3.sh` -> `r3_enum_D{1,2,3}.jsonl`; `r3_compare.py` -> `r3_compare.out`: the author's
  `D <= 3`, 14-`lambda` family, compared histogram by histogram.
* `r3_sample.py`, `r4_nonunit.py`, `run_samples.sh` -> `r3_sample.jsonl`, `r4_nonunit.jsonl`; `r3_summarize.py` ->
  `r3_summary.out`: random unit pencils (`D <= 5`, all candidate special `lambda`) and non-unit stress test.
* Notes: `../referee_theory_notes.md`.
