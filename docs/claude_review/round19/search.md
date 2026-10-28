# Round 19: exact abc-gain search for unit pencils on Mbar_{0,8} (stop-rule test)

Labels: **proved** (proof given here), **exact** (finite exact computation, complete by a proved reduction),
**data/evidence** (computation that no proof uses). Nothing here has been refereed.

## 0. Summary and verdict

**Verdict: the round-18 stop rule fires.** Every unit pencil `f = lambda A/C` on `Mbar_{0,8}` with
`deg A = deg C = D <= 6` and **any** `lambda in Q^*` has `gain(f) <= -32` (exact). Every recorded pencil attaining the
maximum `-32` (all `lambda = ±1` and all special ones) is a forgetful pull-back of a genus-0 Belyi map (data). No pencil with gain `> -32`, let alone `> 0`, exists in
the searched class, so there is no input for the nested-chain extension. **No growth improvement and no uniform
bound follow**; `M <= (2/3 + eps) log R/loglog R` is unchanged.

1. **Excess identity (Lemma 1.1, proved).** `gain = -N0(A,C) + sum_{Y ⊂ F_1} (m_Y - 1) Y.beta`, where
   `N0 = #{T visible : a_T ≠ c_T}` and the sum runs over the components of the fibre over 1. Positive excess needs
   a non-reduced fibre over 1.
2. **`N0 >= 32` for all `D` (Prop 2.2, proved).** So for every `D` and every coprime unit pair, the gain is `<= -32`
   at every non-special `lambda`.
3. **Special `lambda` (Prop 2.1, Thm 2.3, proved).** For `lambda ≠ ±1` the fibre over 1 has no pair component and
   no exceptional component, so only a repeated interior factor of `B_lambda` can help. Every such `lambda` is a
   rational root of the discriminant `Delta_ℓ(lambda) = Disc_D(C|_ℓ - lambda A|_ℓ)` on **one** line `ℓ`. This turns
   the referee's probability-1 candidate list (round 18 §3.4) into a provably complete one, for every `D`.
4. **Fast exact path for `lambda = ±1` (Prop 3.1, proved).** It is cross-checked against the full referee formula
   on every pencil with `D <= 5` (21020 + 256920 + 2556680 pencils with `lambda = ±1`, plus every squarefree
   candidate `lambda`; 0 mismatches; §4, §7.3).
5. **Cross-checks with round 18.** `D <= 4`, `lambda = ±1`: histogram identical entry by entry (43 gain values,
   256920 pencils). `D = 2, 3` special `lambda`: the referee's 17 values, all `-32`.
6. **Results (exact).**

| `D` | pairs `(A,C)` | `lambda = ±1` pencils: max (count) | generic `lambda`: max `-N0` | genuinely special `lambda ≠ ±1`: count, gains | max gain | second value |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 128460 | `-32` (24) | `-32` | 19, all `-32` | **`-32`** | `-40` |
| 5 | 1278340 | `-32` (20) | `-32` | 22, all `-32` | **`-32`** | `-44` |
| 6 | 11313931 | `-32` (24) | `-32` | 31, all `-32` | **`-32`** | `-40` |

7. **Not searched:** `D >= 7` (`D = 7` has 80687282 pairs, beyond the time budget), and non-unit pencils with
   base locus in codimension-2 boundary strata (also admissible in Prop 2.6).

What the stop rule means: the judge's decision rule says to close the abc lane and stop internal rounds. That is a
decision, not a theorem; the abc-pencil question ("does any pencil on `Mbar_{0,8}` have `K^log_rel.beta < 10`?")
stays open for `D >= 7` and for non-unit pencils, and a positive answer would still give only `(L_8^pc)`.

## 1. Setup: parametrisation, normalisation, gain (recalled from theory.md §2.4 + referee corrections)

Everything in this section is copied from `round18/theory.md` §2.4 (Prop 2.6 and the definitions before it) and
`round18/referee_theory.md` §2.6, §3.3, §3.4 (item F); the code conventions are those of
`round18/referee_theory_checks/ref_gain.py` and `ref_gain_lambda.py`.

* **Chart.** Kapranov's model: `Mbar_{0,8}` is the iterated blow-up of `P^5 = P(Q^7 / Q(1,...,1))` along the linear
  spaces `L_T = {x_i = x_j : i, j in T}`, `T ⊂ [7]`, `3 <= |T| <= 6` (98 of them); the 21 hyperplanes `x_i = x_j`
  are the boundary divisors `D_{ij}`. The 119 *visible* `T` are the 21 pairs and the 98 sets `Ts`.
  `beta` is the balanced curve class: `H.beta = 31`, `E_T.beta = 1`, every boundary divisor has `beta`-degree 1,
  `K.beta = -6·31 + sum_T (|T|-2) = 10`.
* **Unit pencil.** `f = lambda A/C` with `lambda in Q^*`, `A, C in Z[x_1..x_7]` products of `D` differences
  `x_i - x_j` (multisets of edges of `K_7`), coprime, i.e. disjoint edge multisets. `B = C - lambda A`
  (for `lambda = p/q` the code uses `B = qC - pA`, which has the same divisor). The fibres over `0, 1, infinity` are
  `{A = 0}`, `{B = 0}`, `{C = 0}` (strict transforms plus exceptional divisors).
* **Orders.** `nu_T(F)` = order of `F` along `L_T` (`|T| >= 3`), or the exponent of `x_i - x_j` (`T = {i,j}`);
  for a unit polynomial `nu_T(A) = a_T := #{edges of A inside T}`; `mu_T = min(a_T, c_T)` (0 for pairs, since the
  edge multisets are disjoint).
* **Gain** (theory.md §2.4, verbatim):

```text
L.beta     = 31 D - sum_{|T|>=3} mu_T,
F_red.beta = #{T visible : max(nu_T(A), nu_T(B), nu_T(C)) > mu_T} + sum_Y (31 deg Y - sum_{|T|>=3} nu_T(Y)),
gain(f)    = L.beta - F_red.beta,
```

  the last sum over the `Q`-irreducible factors `Y` of `B` that are not differences. Prop 2.6 (proved, conditional
  on abc): `gain > 0` implies `(L_8^pc)` off `{B = 0}`. Equivalent forms: `gain = 10 - K^log_rel.beta~ =
  sum_Y ((m_Y - 3)/3) Y.beta` over all components of the three fibres.
* **Symmetry.** `S_7` (relabelling `x_1..x_7`) preserves the chart, the class `beta` and the formula, so `A` may be
  taken up to `S_7`; `C` runs over all multisets of `D` edges disjoint from `A`. (The orbit counts and pair counts
  for `D <= 4` agree with round 18: 1, 3, 8, 22 orbits; 20, 590, 9900, 128460 pairs.)
* **Item F (referee).** Round 18's theory enumeration used `lambda = ±1` only; the referee added every special
  rational `lambda` for `D <= 3` (17 values, all gain `-32`), with a candidate set that was complete only with
  probability 1 (two random lines). `D = 4` with `lambda ≠ ±1` and all `D >= 5` were untested.

### 1.1 The excess identity

**Lemma 1.1 (exact identity).** Let `N0(A,C) = #{T visible : a_T ≠ c_T}`. Let `b_e` be the multiplicity of the
pair factor `e` in `B`, `m_Y` the multiplicity of the interior `Q`-irreducible factor `Y`,
`Y.beta = 31 deg Y - sum_{|T|>=3} nu_T(Y)`, and `S = {T : |T| >= 3, a_T = c_T, nu_T(B) > mu_T}`. Then

```text
gain(lambda A/C) = -N0(A,C) + sum_{e | B} (b_e - 1) + sum_{T in S} (nu_T(B) - mu_T - 1) + sum_Y (m_Y - 1) Y.beta .   (1.1)
```

So `gain = -N0 + sum_{Y ⊂ F_1} (m_Y - 1) Y.beta`, a sum over the components of the fibre over 1 (pair divisors,
exceptional divisors `E_T`, `T in S`, with multiplicity `nu_T(B) - mu_T`, and interior components).

*Proof.* In the formula, `#{T in Ts : max(a,b,c) > mu} = #{T in Ts : a_T ≠ c_T} + #S` (if `a_T ≠ c_T` the max
exceeds the min; if `a_T = c_T = mu_T` the max exceeds `mu_T` iff `nu_T(B) > mu_T`), and the pair count
`#{e in A ∪ C ∪ pairs(B)} = #{pairs T: a_T ≠ c_T} + #supp(pairs of B)`, since `A, B, C` are pairwise coprime (a pair
dividing two of `A, B, C` divides the third). So `F_red.beta = N0 + #S + #supp + sum_Y Y.beta`. On the other side
`nu_T(B) = mu_T` for `T in Ts \ S` (if `a_T ≠ c_T`, the term of lower order survives in `C - lambda A`), so
`sum_{Ts} mu_T = sum_{Ts} nu_T(B) - sum_S (nu_T(B) - mu_T)`. With `D = sum_e b_e + sum_Y m_Y deg Y` and
`sum_{Ts} nu_T(B) = sum_e 30 b_e + sum_Y m_Y sum_{Ts} nu_T(Y)` (each pair lies in exactly 30 sets of `Ts`) this gives
`L.beta = sum_e b_e + sum_S (nu_T(B) - mu_T) + sum_Y m_Y Y.beta`. Subtracting gives (1.1). `[]`

In particular `gain = -N0` whenever the fibre over 1 is reduced. The cross-ratio `(x_1-x_2)/(x_1-x_3)` has
`N0 = 16 + 16 = 32` (the visible `T` containing `{1,2}` but not 3, or `{1,3}` but not 2).

## 2. Special lambda: the generic value, and a provably complete candidate list

Fix a coprime unit pair `(A, C)` of degree `D` (any `D >= 1`). Call `lambda in Q^*` *special* for `(A, C)` if the
fibre over 1 of `lambda A/C` is non-reduced, i.e. some term of (1.1) after `-N0` is nonzero.

**Proposition 2.1 (the generic value; proved, all `D`).** For every `lambda in Q^* \ {±1}`:
1. `B_lambda = C - lambda A` has no pair factor;
2. `S = ∅`, i.e. `nu_T(B_lambda) = mu_T` for every `T in Ts`;
3. hence `gain(lambda A/C) = -N0(A,C) + sum_Y (m_Y - 1) Y.beta` (interior factors only), and
   `gain(lambda A/C) = -N0(A,C)` whenever `B_lambda` is squarefree.

*Proof.* (1) Suppose `x_i - x_j | B_lambda`. Then `{i,j}` is not an edge of `A` or `C` (else it divides both, by
`C = B + lambda A`, contradicting coprimality), and substituting `x_j := x_i` gives
`C|_{x_j = x_i} = lambda A|_{x_j = x_i}`, an identity between two nonzero products of linear forms each of which is
`±(x_u - x_v)` with `u ≠ v` (a factor `x_k - x_j` becomes `x_k - x_i`; none becomes 0 because `{i,j}` is not an
edge). By unique factorisation in `Q[x]` the factors match up to constants, and two forms `±(x_u - x_v)`,
`±(x_{u'} - x_{v'})` are proportional only with ratio `±1`. So `lambda = ±1`, a contradiction.
(2) Let `T in Ts` with `a_T = c_T = mu_T` (if `a_T ≠ c_T`, `nu_T(B) = mu_T` automatically). In the local
coordinates `x_i = x_{i0} + eps·y_i` (`i in T \ {i0}`, `i0 = min T`, `y_{i0} := 0`), each factor `x_i - x_j` of `A` or
`C` has initial form `eps(y_i - y_j)` (both ends in `T`), `x_{i0} - x_j` (only `i` in `T`), `x_i - x_{i0}` (only
`j`), or `x_i - x_j` (neither). So `in_T(A)` and `in_T(C)` are nonzero products of forms `±(u - v)`, `u ≠ v` distinct
variables among `{y_i} ∪ {x_k}` with `y_{i0} = 0`, and `nu_T(B_lambda) > mu_T` iff `in_T(C) = lambda in_T(A)`. As in
(1) this forces `lambda = ±1`. (3) follows from (1), (2) and Lemma 1.1. `[]`

**Proposition 2.2 (`N0 >= 32`; proved, all `D`).** For every coprime unit pair `(A, C)` of equal degree,
`N0(A, C) >= 32`. Hence **for every `D` and every `(A, C)`, all non-special `lambda` give gain `<= -32`**, with
equality for the cross-ratio.

*Proof.* Put `w = 1_A - 1_C` (signed edge multiplicities) and `w(T) = sum_{e ⊂ T} w_e = a_T - c_T` for `T ⊆ [7]`.
Then `w(T) = 0` for `|T| <= 1`, and `w([7]) = D - D = 0`; so every `T` with `w(T) ≠ 0` is visible. Pick an edge
`{i,j}` with `w_{ij} ≠ 0` (one exists since `A ≠ ∅` and `A, C` are disjoint). For each of the `2^5 = 32` sets
`T' ⊆ [7] \ {i,j}`, `w(T' + i + j) - w(T' + i) - w(T' + j) + w(T') = w_{ij} ≠ 0`, so one of these four sets has
`w ≠ 0`. The 32 quadruples are disjoint. `[]`

**Theorem 2.3 (complete candidate set; proved).** Let `ℓ: (s:t) -> s p + t q` be a line (`p, q in Z^7`), and let
`a(s,t) = A(sp + tq)`, `c(s,t) = C(sp + tq)`, binary forms of degree `D`. Let
`Delta_ℓ(lambda) = Disc_D(c - lambda a)`, the discriminant of a binary form of formal degree `D`, a polynomial of
degree `<= 2D - 2` in `lambda`. If `Delta_ℓ ≢ 0`, then every special `lambda ≠ ±1` is a rational root of
`Delta_ℓ`. Consequently the set of special `lambda ≠ ±1` is contained in the finite, explicitly computable set
`Λ_ℓ(A, C) = {rational roots of Delta_ℓ} \ {0, ±1}`, and for `lambda` outside `Λ_ℓ ∪ {±1}` the gain is exactly
`-N0(A, C) <= -32`.

*Proof.* By Prop 2.1 a special `lambda ≠ ±1` has `B_lambda = Y^2 G` with `Y` a nonconstant form (an interior factor
of multiplicity `>= 2`). Then `B_lambda(sp + tq) = c - lambda a` equals `Y(sp+tq)^2 G(sp+tq)`; either
`Y(sp + tq) ≡ 0`, and then `c - lambda a ≡ 0`, or `Y(sp + tq)` is a nonzero binary form of degree `deg Y >= 1`,
which has a zero on `P^1(Qbar)`, a double zero of `c - lambda a`. In both cases `Disc_D(c - lambda a) = 0`, since the
discriminant of a binary form of formal degree `D` vanishes exactly when the form is `0` or has a repeated zero on
`P^1` (including the point at infinity, where the order of the zero is the drop in degree). The degree bound follows because
`Disc_D` is homogeneous of degree `2D - 2` in the coefficients, each linear in `lambda`. `[]`

*Remarks.* (i) The referee's construction (round 18 §3.4) intersected the roots for two random lines and called
the result complete "with probability 1"; Theorem 2.3 shows **one** line with `Delta_ℓ ≢ 0` already gives a
provably complete superset, and any number of further lines may be used to discard candidates (each line gives a
necessary condition). (ii) The candidate types (i) and (ii) of the referee (proportional initial forms along `L_T`;
proportional restrictions to `x_i = x_j`) produce only `lambda = ±1` (Prop 2.1), so for `lambda ≠ ±1` only the
discriminant type (iii) matters. (iii) `Delta_ℓ1 ≢ 0`, and `t`-degree exactly `D` of `c - lambda a` on line 1, are checked for every pair; both
counters (`DISC_ZERO_line1`, `deg_drop_line1`) stayed 0 in every run, at `D = 4, 5, 6`.

**How the search uses it.** For each pair `(A, C)`: compute `Delta_ℓ1` exactly (`fmpz_mpoly.discriminant` in
`Z[t, lambda]`; on line 1 the `t`-degree of `c - lambda a` is `D`, which is checked, so the affine discriminant equals
`Disc_D`), factor it over `Z`, take the linear factors. For each candidate `lambda = p/q`: if `q c - p a` is
squarefree as a binary form of degree `D` on line 2 or line 3, then `B_lambda` is squarefree (by the same argument)
and `gain = -N0` (Prop 2.1(3)); otherwise the exact gain is computed by the referee's formula
(`exact.exact_gain`). The squarefree test on a line: the dehomogenised `b(t)` has degree `>= D - 1` and
`gcd(b, b') = 1`.

## 3. Algorithm, symmetry reduction, chunking, exactness

**Enumeration.** `A` runs over `S_7`-orbit representatives of `D`-edge multisets of `K_7` (brute-force canonical
relabelling, as in `ref_gain.py`; for `D = 6` the representatives are produced by augmentation, `orbits_aug.py`,
which reproduces the brute-force lists exactly for `D <= 5`). `C` runs over **all** multisets of `D` edges from the
`21 - #supp(A)` edges not in `A` (itertools order), so every coprime unit pair is covered up to `S_7`. For each pair
three things are computed:

1. **Generic value** `-N0(A, C)` (numpy count of `a_T ≠ c_T` over the 119 visible `T`). By Thm 2.3 this is the
   exact gain for every rational `lambda` outside `{±1} ∪ Λ_ℓ1(A, C)`.
2. **`lambda = s = ±1`, exact gain**, via the fast path of Prop 3.1 below. The fast path never failed at `D <= 6`
   (the fallback, the referee's full formula, was never needed for `lambda = ±1`), and it was checked against the
   full formula on every pencil with `D <= 4` and on every pencil with `D = 5` (§4).
3. **`lambda ≠ ±1`, special values**: all rational roots of `Delta_ℓ1`, each decided as in §2: squarefree on line 2
   or 3 (gain `-N0`), otherwise the exact gain by the referee formula (`exact.exact_gain`, which factors
   `qC - pA` over `Q` with flint and computes every `nu_T(Y)` by exact substitution). The exact records of all
   genuinely special `lambda` (non-squarefree `B_lambda`) are written to disk.

**Proposition 3.1 (fast path for `lambda = s = ±1`; proved).** Let `S_s = {T in Ts : in_T(C) = s·in_T(A)}`.
(a) `S_s` is exactly the set `S` of Lemma 1.1, and `mu_T = 0` for `T in S_s`.
(b) If `B = C - sA` restricted to some line is a squarefree binary form of formal degree `D`, then `B` is
squarefree and `gain(sA/C) = -N0 + sum_{T in S_s} (nu_T(B) - 1)`.

*Proof.* (a) For `T` with `a_T ≠ c_T`, `nu_T(B) = mu_T` and the initial forms have different degrees, so `T` is
in neither set. For `a_T = c_T`, `nu_T(B) > mu_T` iff the initial forms cancel, i.e. `in_T(C) = s in_T(A)`. If this
holds, the factors `eps(y_i - y_j)` coming from edges inside `T` must agree as multisets; the map
`{i,j} -> y_i - y_j` (with `y_{i0} = 0`) is injective on edges inside `T`, and `A`, `C` share no edge, so there are
none: `a_T = c_T = 0`. (b) Squarefreeness as in the proof of Thm 2.3. Then all `b_e = 1`, `m_Y = 1`, and (1.1) leaves
`sum_{T in S}(nu_T(B) - mu_T - 1)` with `mu_T = 0`. `[]`

*Implementation of (b).* The initial form of a unit polynomial is computed combinatorially (`common.red_init`: each
factor `x_i - x_j` is replaced by its initial factor, normalised with a sign); candidates are prefiltered by a
64-bit additive hash of the multiset of initial factors and then compared **exactly** (hash equality is necessary
for multiset equality, so the prefilter loses nothing). `nu_T(B)` for `T in S_s` is computed exactly by flint
substitution (`exact.nuT`). The squarefree test uses three fixed integer lines (seed 2026, coordinates in
`[-60, 60]`).

**Exactness.** All gains are integers computed from exact integer polynomial arithmetic (python-flint 0.9:
`fmpz_poly`, `fmpz_mpoly` factorisation, discriminant, composition). No floating point and no random
specialisation enters a gain value. Randomness enters only through the choice of the lines and the hash weights,
and neither can make a result wrong: a bad line can only make the squarefree test fail (then the full exact formula
is used) or make `Delta_ℓ1 ≡ 0` (this is detected and logged; it never happened), and the hash only prefilters.

**Chunking and resources.** One process per (orbit, block of 20000 `C`'s) for `D = 5` (100 chunks), 2 processes at
a time (`run.sh`, `xargs -P 2`), each under `prlimit --as=3500000000` and `timeout 3000`; each chunk appends one JSON
line (histograms, all records with gain `> -40`, every genuinely special `lambda`, counters) and `todo.py` lists
the chunks still missing, so a restart loses at most the running chunks (`D = 4`: one chunk per orbit; `D = 6`:
334 chunks of 50000; verification runs: blocks of 10000 / 5000). No chunk failed (no `FAIL` line in any `.err`
file); the `D = 5` verification was paused once and resumed from its output without loss.

## 4. Cross-checks against round-18 outputs

All of these were run **before** the `D = 5` search was interpreted.

| check | result |
| --- | --- |
| `D = 1, 2, 3`, `lambda = ±1`: histogram vs referee table (`referee_theory.md` §3.3) | identical: `D = 2`: `-32` (24), `-40` (12), `-48` (152), `-52` (144), `-56` (272); `D = 3`: `-32` (20), `-40` (4), `-44` (2), `-48` (208), `-52` (188), `-56` (876); 40 / 1180 / 19800 pencils |
| `D = 1..3`, fast path vs full exact formula on **every** pencil (`verify` mode) | 21020 / 21020 agree |
| `D = 2, 3`, special rational `lambda ≠ ±1` vs referee `ref_gain_lambda_D{2,3}.out` | identical counts: 6 at `D = 2`, 11 at `D = 3` (17 in all), every one of gain `-32`; e.g. `lambda = -1/4` for `A = (x_1-x_2)^2`, `C = (x_1-x_3)(x_2-x_3)` |
| `D <= 3`, generic value `max(-N0)` vs referee | `-32`, as the referee found |
| `D = 4`, `lambda = ±1`: full histogram vs referee `ref_gain_D4_chunks.jsonl` (66 chunks aggregated), `compare_D4.py` | **identical entry by entry**, all 43 distinct gain values, 256920 pencils (`out/compare_D4.out`) |
| `D = 4`, fast path vs full exact formula on every `lambda = ±1` pencil, and exact formula at every squarefree candidate `lambda` | 256920 / 256920 and 635 / 635 agree (`out_verify/D4v.jsonl`) |
| `D = 5`, fast path vs full exact formula on every `lambda = ±1` pencil and every squarefree candidate | 2556680 / 2556680 and 826 / 826 agree (§7.3) |
| `D = 5`, candidate `lambda` from a different set of lines | same 22 special values, same histograms (§7.3) |

The referee's numbering of pencils counts both signs `s = ±1` for each pair `(A, C)`; so does this search.

## 5. Run-time estimate and the searched class

| `D` | `A`-orbits | pairs `(A, C)` | `lambda = ±1` pencils | fast search (CPU / wall, 2 processes) | full exact re-verification |
| --- | --- | --- | --- | --- | --- |
| 4 | 22 | 128460 | 256920 | 13.7 s / 16 s | 759 s CPU |
| 5 | 60 | 1278340 | 2556680 | 223 s / 2 min | 13652 s CPU (all pencils, §7.3) |
| 6 | 173 | 11313931 | 22627862 | 3132 s / 27 min | not run (fast path only; proved by Prop 3.1) |
| 7 | 471 | 80687282 | 161374564 | estimated `>= 22000 s` CPU, `>= 3 h` wall | **not searched** |

The round-18 referee needed 66 processes for `D = 4` because every pencil went through the full formula (about
2.5 ms per pencil, plus a memory leak in python-flint). The fast path of Prop 3.1 spends about 0.09 ms per pencil,
and the full formula is only called for the few non-squarefree `B`. So the **entire** `D = 5` class was searched;
no sub-class restriction was needed. The searched class is exactly: all `S_7`-orbits of `A` (60), all `C` disjoint
from `A`, all `lambda in Q^*` (`±1` explicitly, the rest through Thm 2.3). The same complete class was also
searched at `D = 6` (§9). `D = 7` (estimate from the measured `D = 6` cost of 0.28 ms per pair) would need over
3 hours of wall clock on 2 processes and was **not searched**; nor was any `D >= 7`.

## 6. Results: D = 4, all rational lambda (exact; complete by Thm 2.3)

Files: `out/D4.jsonl` (22 chunks, one per orbit), `out/D4_analysis.json` (exact fibre decompositions).

* **Pairs:** 128460 coprime unit pairs `(A, C)` (22 `A`-orbits). `Delta_ℓ1 ≢ 0` and `t`-degree `D` on line 1 for
  every pair (counters `DISC_ZERO_line1 = 0`, `deg_drop_line1 = 0`).
* **Generic `lambda`** (all rational `lambda` outside `{±1} ∪ Λ_ℓ1`): gain `-N0`; histogram top
  `-32` (12 pairs), `-40` (7), `-48` (159), `-50` (5), `-52` (119), `-54` (3), `-56` (818); maximum `-32`
  (Prop 2.2 proves `<= -32` for all `D`).
* **`lambda = ±1`:** 256920 pencils, identical to round 18; maximum `-32` (24 pencils), next `-40` (14), `-48` (318).
* **Special `lambda ≠ ±1`:** `Delta_ℓ1` produced 654 rational candidates `≠ 0, ±1`; 635 of them are discarded
  rigorously (`B_lambda` squarefree on line 2 or 3; gain `-N0`, confirmed by the full formula for all 635). The
  remaining **19 are genuinely special** (`B_lambda` non-squarefree). **All 19 have gain exactly `-32`.** Up to
  `S_7` relabelling and the symmetry `f -> 1/f`, they are three Belyi compositions with the cross-ratio
  `z = (x_1 - x_2)/(x_1 - x_3)` (or a 4-point cross-ratio):

| `lambda` | `A` | `C` | `B_lambda` (interior factors, multiplicity, `Y.beta`) | `L.beta` | `F_red.beta` |
| --- | --- | --- | --- | --- | --- |
| `-27/256` | `(x_1-x_2)^4` | `(x_1-x_3)^3 (x_2-x_3)` | `(x_1 + 3x_2 - 4x_3)^2` (16) · quadric (32) | 64 | 96 |
| `1/16` | `(x_1-x_2)^4` | `(x_1-x_3)^2 (x_2-x_3)^2` | `(x_1 + x_2 - 2x_3)^2` (16) · quadric (32) | 64 | 96 |
| `-1/4` | `(x_1-x_2)^2 (x_3-x_4)^2` | `(x_1-x_3)(x_1-x_4)(x_2-x_3)(x_2-x_4)` | `(2x_1x_2 - x_1x_3 - x_1x_4 - x_2x_3 - x_2x_4 + 2x_3x_4)^2` (16) | 32 | 64 |

  (The full list of 19, with labels `x0..x6 = x_1..x_7` and exact `B`, is in `out/D4_analysis.json`. They are:
  `A = (x_1-x_2)^4` with `C = (x_1-x_k)^3(x_2-x_k)` or `(x_1-x_k)(x_2-x_k)^3`, `lambda = -27/256` (10 pencils,
  `k = 3..7`), and `C = (x_1-x_k)^2(x_2-x_k)^2`, `lambda = 1/16` (5 pencils); the two pencils with
  `C = (x_2-x_3)^4` and `A = (x_1-x_2)^3(x_1-x_3)`, `lambda = -256/27`, or `A = (x_1-x_2)^2(x_1-x_3)^2`,
  `lambda = 16` (these are `1/f` of the former up to relabelling); and the two 4-point pencils `lambda = -1/4`,
  `-4`.) In each case the fibre over 1 has exactly one double
  component, of `beta`-degree 16, and (1.1) reads `gain = -48 + 16 = -32`.
* **Structure of the maximum (Data).** Every pencil of gain `-32` at `D = 4` (24 with `lambda = ±1`, 19 special) is a
  pull-back `h ∘ chi_Q` along a forgetful map `chi_Q : Mbar_{0,8} -> Mbar_{0,4}`: `A ∪ C` is supported on 3 labels
  (the fourth point is `∞`), or on 4 labels with `A/C` invariant under `x -> 1/x` (checked symbolically,
  `analyze.py`). Such a pencil is a morphism of `Mbar_{0,8}`, so by theory.md §2.4(1)
  `gain = 16 deg h - 16 #h^{-1}{0,1,∞} <= -32`, with equality iff `#h^{-1}{0,1,∞} = deg h + 2`, i.e. `h` is a
  rational function branched only over `{0, 1, ∞}` (a genus-0 Belyi map), e.g. `h(z) = z^4` (`lambda = ±1`).

**Conclusion for `D = 4`.** Every unit pencil `lambda A/C` with `D = 4` and `lambda in Q^*` has gain `<= -32`.
This is a complete statement: for `lambda` outside the finite set `{±1} ∪ Λ_ℓ1(A, C)` the gain is `-N0 <= -32` by
Prop 2.1/2.2 and Thm 2.3, and every `lambda` in the finite set was evaluated exactly.

## 7. Results: D = 5, all rational lambda (exact; complete by Thm 2.3)

Files: `out/D5.jsonl` (100 chunks), `out/D5_analysis.json`.

### 7.1 Distribution

* **Pairs:** 1278340 (60 `A`-orbits); `DISC_ZERO_line1 = 0`, `deg_drop_line1 = 0`; no fallback to the full
  formula was needed for `lambda = ±1` (all 2556680 `B = C ∓ A` are squarefree, certified on a line).
* **Generic `lambda`:** gain `-N0`; top of histogram `-32` (10 pairs), `-44` (1), `-48` (167), `-52` (44),
  `-54` (14), `-56` (998), `-58` (113), `-60` (613). Maximum `-32`.
* **`lambda = ±1`:** 2556680 pencils; top of histogram `-32` (20), `-44` (2), `-48` (334), `-52` (88),
  `-54` (28), `-56` (1996), `-58` (226), `-60` (1226), `-61` (62), `-62` (860), `-63` (84), `-64` (10732),
  `-65` (648), `-66` (7838). Maximum `-32`. (Full histogram in `out/D5_analysis.json`.)
* **Special `lambda ≠ ±1`:** 848 rational candidates from `Delta_ℓ1`; 826 discarded rigorously as squarefree
  (gain `-N0`); **22 genuinely special, every one of gain exactly `-32`.**

### 7.2 The 22 special values

All 22 are Belyi compositions with a cross-ratio, with one double component of `beta`-degree 16 in the fibre over 1
(`N0 = 48`, excess 16, `L.beta = 80`, `F_red.beta = 112`):

| `lambda` | `A` | `C` (k = 3..7) | `B_lambda` | count |
| --- | --- | --- | --- | --- |
| `256/3125` | `(x_1-x_2)^5` | `(x_1-x_k)(x_2-x_k)^4` | `(4x_1 + x_2 - 5x_k)^2 ·` cubic (`Y.beta = 48`) | 5 |
| `-256/3125` | `(x_1-x_2)^5` | `(x_1-x_k)^4(x_2-x_k)` | `(x_1 + 4x_2 - 5x_k)^2 ·` cubic | 5 |
| `108/3125` | `(x_1-x_2)^5` | `(x_1-x_k)^3(x_2-x_k)^2` | `(2x_1 + 3x_2 - 5x_k)^2 ·` cubic | 5 |
| `-108/3125` | `(x_1-x_2)^5` | `(x_1-x_k)^2(x_2-x_k)^3` | `(3x_1 + 2x_2 - 5x_k)^2 ·` cubic | 5 |
| `3125/256` | `(x_1-x_2)^4(x_1-x_3)` | `(x_2-x_3)^5` | `(5x_1 - x_2 - 4x_3)^2 ·` cubic | 1 |
| `-3125/108` | `(x_1-x_2)^3(x_1-x_3)^2` | `(x_2-x_3)^5` | `(5x_1 - 2x_2 - 3x_3)^2 ·` cubic | 1 |

(exact cubics in `out/D5_analysis.json`; e.g. for `lambda = 256/3125`, `k = 5`:
`16x_1^3 - 88x_1^2x_2 + 40x_1^2x_5 + 203x_1x_2^2 - 230x_1x_2x_5 + 75x_1x_5^2 - 256x_2^3 + 565x_2^2x_5 - 450x_2x_5^2
+ 125x_5^3`.) In the cross-ratio `z`, these are `h(z) = c z^5/(z-1)^j`-type maps (`j = 1..4`) with the unique
critical value `1` of the polynomial part, i.e. genus-0 Belyi maps.

**No pencil with gain `> -32` exists at `D = 5`, and none with gain in `(-44, -32)`**: the second-largest value is
`-44`, attained by 2 pencils with `lambda = ±1` (generic value `-44` for 1 pair).

**Structure of the maximum (Data).** All 42 pencils of gain `-32` at `D = 5` (20 with `lambda = ±1`, 22 special)
are forgetful pull-backs `h ∘ chi_Q` with `h` a genus-0 Belyi map (`analyze.py`, flag
`all_gain_-32_records_are_forgetful_pullbacks = True`).

### 7.3 Full exact re-verification of D = 5

* **Every** `lambda = ±1` pencil at `D = 5` (2556680) was recomputed with the referee's full formula
  (`exact.exact_gain`: `Q`-factorisation of `B`, exact `nu_T` of every interior factor by substitution):
  **2556680 / 2556680 agree** with the fast path (`out_verify/D5v.jsonl`, 305 chunks, 13652 s CPU).
* **Every** squarefree candidate `lambda` (826) was recomputed with the full formula: **826 / 826** have gain `-N0`.
* **Line independence.** The whole `D = 5` search was re-run with three different lines (`LINESEED=777`,
  `out_verify/D5_seed777.jsonl`, `compare_seed.py`): 840 candidates instead of 848, but the **same 22 genuinely
  special `(A, C, lambda)`**, and identical `lambda = ±1` and generic histograms (`out_verify/compare_seed_D5.out`).
  Thm 2.3 already guarantees this; the re-run checks the implementation.

## 8. Stop rule

The round-18 judge's rule (`round18/README.md` §7.3): *if every unit pencil with `D = 5` (all `lambda`) and
`D = 4` (`lambda ≠ ±1`) has gain `<= -32`, close the abc lane and stop internal rounds; if one has positive gain,
the next task is the nested-chain extension.*

* `D = 4`, `lambda ≠ ±1`: every unit pencil has gain `<= -32` (§6; complete by Thm 2.3; maximum `-32`, attained by
  19 special values and by every pair with `N0 = 32` at generic `lambda`).
* `D = 5`, all `lambda in Q^*`: every unit pencil has gain `<= -32` (§7; complete by Thm 2.3; maximum `-32`).

**The stop rule fires.** No pencil with positive gain, and none with gain `> -32`, exists among unit pencils of
degree `<= 5` with any rational `lambda`; the additional complete run at `D = 6` (§9) gives the same maximum `-32`.

What this does and does not establish (labels as in the record):
* **Exact (finite exact computation, made complete by Thm 2.3 and Prop 2.1, 2.2):** `gain(lambda A/C) <= -32` for
  every unit pencil on `Mbar_{0,8}` with `deg A = deg C <= 6` and `lambda in Q^*` (for `D = 5` re-verified pencil by
  pencil with the referee formula, §7.3; for `D = 6` via the proved fast path).
* **Proved for all `D`:** by Prop 2.1-2.2 and Thm 2.3, for **every** degree `D` a unit pencil can beat `-32` only at one of the finitely many special `lambda` (`±1` or a rational root of
  `Delta_ℓ`), and only through a non-reduced fibre over 1.
* **Not established:** unit pencils with `D >= 7` (not searched), and non-unit
  pencils whose base locus lies in codimension-2 boundary strata (theory.md §2.4, which Prop 2.6 also allows).
  "Closing the abc lane" is the judge's decision rule, not a theorem that no positive-gain pencil exists.
* **Evidence (Data):** every pencil attaining `-32` at `D = 4, 5, 6` is a forgetful pull-back `h ∘ chi_Q` with `h`
  a genus-0 Belyi map, i.e. exactly the morphism extremals of theory.md §2.4(1). The next value after `-32` is
  `-40` at `D = 2, 3, 4, 6` and `-44` at `D = 5`, so the maximum is isolated and shows no drift upward with `D`: the top of the histogram is `-32` at every `D = 1, ..., 6`.
* **Consequence for the program:** none for growth or uniformity. `M <= (2/3 + eps) log R/loglog R` stands; the
  uniform bound remains conditional on Vojta (`V_8`, via `(L_8)`); the abc-conditional route via pencils produced no
  candidate, and even a positive pencil would have given only `(L_8^pc)` (round-18 README §3.3).

## 9. Beyond the stop rule: D = 6, all rational lambda (exact; complete by Thm 2.3)

Files: `out/D6.jsonl` (334 chunks of 50000 `C`'s), `out/D6_analysis.json`; orbit representatives
`reps_D6.json` from `orbits_aug.py` (augmentation; reproduces the brute-force lists for `D <= 5`).

* **Pairs:** 11313931 (173 `A`-orbits); `DISC_ZERO_line1 = 0`, `deg_drop_line1 = 0`; all 22627862 `B = C ∓ A`
  squarefree (fast path throughout, no fallback).
* **Generic `lambda`:** top `-32` (12 pairs), `-40` (8), `-44` (1), `-48` (254), `-50` (5), `-52` (220); max `-32`.
* **`lambda = ±1`:** top `-32` (24), `-40` (16), `-44` (2), `-48` (508), `-50` (10), `-52` (440), `-54` (26),
  `-56` (2854); max `-32`.
* **Special `lambda ≠ ±1`:** 1076 candidates, 1045 squarefree, **31 genuinely special, all of gain `-32`**:
  `lambda = -3125/46656` (10), `16/729` (10), `-1/64` (5) with `A = (x_1-x_2)^6`; their reciprocals
  `-46656/3125`, `729/16`, `-64` (1 each) with `C = (x_2-x_3)^6`; and three 4-point pencils `lambda = -4/27, 4/27,
  27/4` (e.g. `A = (x_1-x_2)^3(x_3-x_4)^3`, `C = (x_1-x_3)^2(x_1-x_4)(x_2-x_3)(x_2-x_4)^2`, `lambda = -4/27`:
  `B` = quadric (`Y.beta = 16`) · (quadric)^2 (`Y.beta = 16`), `L.beta = 48`, `F_red.beta = 80`).
* **Maximum gain at `D = 6`: `-32`**; every one of the 55 pencils with gain `> -40` has gain exactly `-32` and is a
  forgetful pull-back `h ∘ chi_Q` (Data, `analyze.py`).
* The `D = 6` gains rest on Prop 3.1 (fast path, proved) and on the referee formula at the 31 special values; the
  fast path was cross-checked against the full formula on every pencil with `D <= 4` (and at `D = 5`, §7.3), not
  on the `D = 6` pencils themselves.

## 10. Files

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/`:

* `search.md` (this write-up), `search_notes.md` (parked derivations).
* `search_checks/common.py`: edges, the 98 sets `Ts`, `S_7` canonical form, `N0`, combinatorial initial forms.
* `search_checks/exact.py`: the referee's full gain formula (exact; flint factorisation and substitution).
* `search_checks/search.py`: one chunk of the search (generic value, `lambda = ±1` fast path, special `lambda`
  via `Delta_ℓ1`); `VERIFY=1` re-computes every `lambda = ±1` pencil and every squarefree candidate with the full
  formula; `LINESEED` changes the lines.
* `search_checks/run.sh`, `todo.py`: resumable driver (2 processes, `prlimit --as=3500000000`, `timeout 3000`).
* `search_checks/orbits_aug.py`, `reps_D{1..7}.json`: orbit representatives.
* `search_checks/agg.py`, `analyze.py`, `compare_D4.py`, `compare_seed.py`: aggregation, exact fibre decompositions, and the
  entry-by-entry comparison with round 18.
* Outputs: `search_checks/out/D{4,5,6}.jsonl`, `out/D{4,5,6}_analysis.json`, `out/compare_D4.out`;
  verification runs `search_checks/out_verify/D4v.jsonl`, `D5v.jsonl` (+ `D5v_agg.out`), `D5_seed777.jsonl`,
  `compare_seed_D5.out`; test runs `search_checks/test2/`.
