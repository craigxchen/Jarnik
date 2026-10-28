# Referee report on round21/prove.md

Referee: adversarial pass on `round21/prove.md`. Scripts and outputs: `round21/referee_prove_checks/`.

## 0. Verdict summary

**The partial result survives, with its scope corrected.** No formal theorem, proposition or lemma is false. Every
proof I re-derived holds with its quantifiers. Every number I recomputed independently matches: the smooth scan
counts and table, the record-cluster algebra, `#v`, `E` and the Monte Carlo values, and all constants. The write-up
correctly claims **no** growth in the window.

What does not survive is part of the narrative that turns the formal results into "proved reasons" why the large
sieve fails, and one "proved" label.

| item | verdict |
| --- | --- |
| Thm 2.1, Cor 2.2, Cor 2.3 | correct; essentially `stability.md` Thm 3.3(3) plus `omega <= pi(Y;4,1)` (minimal novelty; the record's Prop 3.4 is sharper for `C <= 1`) |
| Remark 2.4 ("counting provably cannot pass exponent 1") | **mislabelled**: a heuristic, not a theorem |
| Lemma 3.1, Lemma 3.2, Prop 3.3 | correct (brute-force and eigenvalue checks pass) |
| Lemma 4.1, Prop 4.2, Cor 4.3 | correct as stated |
| (F1) "the large sieve is spreading and nothing else" | correct only up to an absolute factor (`< 23`); sharp-constant sieves, uncapped in the overview, are not covered |
| (F2) "the short-arc condition never forces concentration" | **overstated (major)**: proved only for unpinned profiles (`rank B <= s - 2`, pairwise independent) and only against realisation **without** the gaps. In the pinned (typical) case concentration is forced below a profile-dependent threshold, and Remark 4.4's "the gaps already give more" is a magnitude comparison, not an implication |
| Prop 5.2 | correct; the summary's error term `M^(-eta/2)` should be `M^(-eta/4)(log M)^(1/4)` |
| Prop 5.3, Cor 5.4 | correct as conditional statements; **(F3) "the window is not special for robust inputs" is overstated (major)**: inside the window `gamma > 1/(1+eta) > 1/2` is needed, so square-root-quality (GRH) prime sums and `x^(1/2)`-error box counts are **not** covered anywhere in the window; nor are primes in `[c M log M, M^(1+eta)]` |
| Prop 6.1(1)-(3), 6.2 | correct |
| Prop 6.1(4) | uses Cor 2.2 (actual clusters) for `A_0` fakes; transfer unproved (gap, likely repairable). "All refereed fakes have `sum p^(-1/2) <= 1/32`" is **false** (`construct.md` Thm 2 reaches `21^(-1/2) = 0.218`) |
| Prop 8.1 | correct, near-tautological; the §8.2 "equivalence" holds only for `A_0` alone |
| §7 evidence | numbers reproduced; reading "best `C_k` per decade increases beyond `10^9`" is false as worded (non-monotone, in the write-up's own table); "window populated" holds only at `delta = 0`, and no record cluster lies in the window for any `delta > 0` |
| Owner constraints | no hard violation; the pinned-case narrative leans on a template-dependent threshold (soft) |

## 1. Theorem 2.1 / Corollary 2.2 (counting bound below the window)

**Re-derivation.** Item 1: by `round8/stability.md` Thm 3.3(2) (checked in the record: proof lines 199-220, valid for
`M >= 4`, any `C`) each `e_p log p <= 2 beta'`; summing over the `s` primes gives `W' <= 2 beta' s`, and (3.0) gives
`2 beta' < 4W'/M + 8 log^+ C`. Item 2: divide `W' < s(4W'/M + 8 log^+ C)` by `W'` and rearrange:
`M < 4s(1 + 2M log^+ C/W') <= 4 pi(Y;4,1)(1 + ...)`. Correct. Item 3 needs `3 log M > K'` strictly; the
threshold `M_1(C)` of Cor 2.2 is defined with `>=`, but when `C > 1` the extra `8 log^+ C > 0` makes it strict, and
when `C <= 1` item 3 is not needed. Item 4: `W' <= 2 beta' s < Ds` with `W' = A M log M`; the limit
`(4A - 6)/A` and the value `2` at `A = 3` are right.

Cor 2.2: for `M >= M_1(C)`, `M < 5 s <= 5 pi(Y;4,1) <= 5 pi(Y) < 6.2753 Y/log Y`; the two cases
`Y <= (log R)^(1/2)` / `Y > (log R)^(1/2)` give `12.56 max((log R)^(1/2), Y)/loglog R`. Correct; `M_1(C)` depends on
`C` only. Cor 2.3: correct (the maximum must be `Y` since `eps log R > 12.56 (log R)^(1/2)`), and
`log R >= (3/2) M log M - (K'/2) M` is (I1) with `W' <= 2 log R`.

**Assessment.** The statements are correct, uniform in the arc's position, and the thresholds depend only on `C`
(and `eps`). Two remarks, neither an error:
* *Novelty is minimal.* Items 1-2 are `stability.md` Thm 3.3(3) (`omega(N') > W'/(4W'/M + 8 log^+ C)`) composed with
  the tautology `omega(N') <= pi(P_+(N'); 4, 1)`. For `C <= 1` the record's Prop 3.4
  (`omega(N') >= (M/2)(1 - M^(-1/3))^2`) already gives the sharper `M <= (2 + o(1)) pi(Y; 4, 1)` at every rate, not
  only at `A -> 3`. The restricted growth statement "`P_+(N') = o(log R)` implies `M = o(log R/loglog R)`" is
  correct and worth recording, but it should be labelled a corollary of `stability.md` §3.
* *Prop 3.3(3) cites Cor 2.2 for `s > M/5`*; Cor 2.2 states `M < 5 pi(Y;4,1)`, and `s > M/5` comes from its proof
  (`M < 5s`). Cosmetic.

**Remark 2.4 is mislabelled.** The table labels it "proved", and the summary says "counting provably cannot get past
exponent 1". It is not a theorem. To show that no statement depending only on the number and sizes of the primes
can exclude near-extremal clusters with primes in `[M^(1+eta), M^(2-delta)]`, one must exhibit objects (actual
clusters or models satisfying the other inputs) with those numbers and sizes. That is exactly the open existence
question of §§5-6 (no window fake is known). "There are `M^(1+eta-o(1))` times more admissible primes than a cluster
uses" is a heuristic about counting arguments, not a proof that they fail. **Relabel it as a heuristic remark.**

## 2. Lemma 3.2 / Prop 3.3 (large sieve = arc counts)

**Lemma 3.1.** Steiner's formula for the convex sector (`eta <= pi/2`) gives area `eta Y/2`, perimeter
`(2 + eta) Y^(1/2)`, `r = 2^(-1/2)`; `(2 + pi/2)/sqrt 2 = 2.52493 <= 2.53`, `pi/2 = 1.5708 <= 1.58`. Brute force
(`referee_prove_checks/consts_ls.py`, `Y` up to `2 * 10^4`, 2000 random and anchored arcs, `eta` from `1/Y` to
`pi/2`): the worst ratio count/bound is `0.978`, so the bound holds and is nearly tight at small `Y`. Correct.

**Lemma 3.2.** Re-derived. The class decomposition is right: in a sorted list, if `x_(i+m) - x_i < delta`, the
half-open arc `[x_i, x_i + delta)` holds `m + 1` points. The wrap-around gap holds `m + rho + 1 >= m + 1` points. There
are `m + rho < 2m` classes, and singletons are trivially `delta`-separated. Montgomery-Vaughan is quoted with its true
hypotheses: `delta`-separated points mod 1, `L` consecutive frequencies, constant `L - 1 + delta^(-1)`. With
`L = 2K + 1` this is `2K + delta^(-1)`. Cauchy-Schwarz over fewer than `2m` classes is fine. The converse uses
`|k(x_n - x_n')| <= 1/8`, so `cos >= 2^(-1/2)`. The factor `16 sqrt 2 = 22.63 < 23` follows from covering a half-open arc
of length `1/(2K)` by four closed arcs of length `1/(8K)`. Numerical check: 60 random, clustered and lattice point
sets, `n <= 150`, `K <= 30`. The top eigenvalue of the Dirichlet-kernel Gram matrix lies between the two bounds in
every case (0 violations). (When `m` is not an integer, as in Prop 3.3, use `floor(m)`. Cosmetic.)

**Prop 3.3.** Items 1-2: with `delta = 1/(2K)` in `x` (a `theta`-arc of length `pi/(4K)`),
`m = pi Y/(8K) + 2.53 Y^(1/2) + 1.58` and `8Km = pi Y + 20.24 K Y^(1/2) + 12.64 K`. Correct. Item 3 is correct
arithmetic (`s > M/5`, `Y^(1/2) <= M^(1 - delta/2)`).

**Scope of (F1).** Lemma 3.2 is a true and useful equivalence: the optimal dual large-sieve constant at
frequencies `|k| <= K` is comparable, within `[2^(-1/2), 32]`, to `K` times the maximal arc count at scale `1/(8K)`.
Two limits must be stated wherever (F1) is used as "the large sieve is spreading and nothing else":
1. *Constants.* The equivalence loses a factor up to 23. It says nothing about large-sieve statements with
   sharp constants, such as `K + delta^(-1)` exactly or `X + K` for `K >> X^(1/2)`. The overview lists the sharp
   `(X + K)` large sieve as **uncapped** (§5.5: "up to constants equivalent to the max-scale cross-gap"). Cor 4.3's
   comparison is likewise only "within a factor 3". So (F1) covers large sieves *with constants fixed up to an
   absolute factor*, not the large sieve as such.
2. *Which point set.* (F1) applies to whatever set of angles the sieve is run on. Prop 3.3 runs it only on the
   cluster's prime angles `{theta_j}`. On other sets built from the cluster (point angles, monomials `<v, theta>`),
   Lemma 3.2 again reduces it to arc counts. For monomial sets those counts are controlled, up to constants, by the
   pairwise spacing that the gaps (I5) already give, so I found no counterexample to the write-up's reading.
   But the write-up does not prove this; it treats only the prime-angle set.

## 3. Prop 4.2 / Cor 4.3 / Remark 4.4 (no forced concentration)

**Lemma 4.1.** (a) `Real(Delta)` descends to `T^s`, because `<2a_x - e, (pi/2) e_j> in (pi/2)Z`. The translation
identity is right and `theta_0` is kept free (not set to zero). (b) The annihilator of `H` is `B^T Z^(M-1)` by double
annihilation, since every subgroup of the discrete group `Z^s` is closed. So `H^ = Cok`, and the image of `h -> h_j` is
the subgroup of order `D_j`. The bound `2 eta/pi + 1/D_j` is right. (c) If `<e_i> cap <e_j> = 0`, the annihilator of
the image of `(h_i, h_j)` is `D_i Z x D_j Z`, so the image is a product and Haar pushes forward to the product
measure. Correct.

**Prop 4.2.** Re-derived. `Var N_I <= E N_I <= E_l` uses pairwise independence. Chebyshev with
`t_l^2 = 2L n_J(1 + 2^l/D_0)` gives `2^l E_l/t_l^2 = 1/(2L)`, and the union over levels `1..L` is `<= 1/2`. With
`2^l <= D_0`: `N_I < 2 n_J 2^(-l) + 2 (L n_J)^(1/2)`. An arc of length in `((pi/2)2^(-l-1), (pi/2)2^(-l)]` meets at
most 2 grid arcs, and `4 n_J 2^(-l) < (16/pi) n_J eta`. For `eta > pi/4` the trivial bound `n_J <= (16/pi) n_J eta`
covers it. The large-sieve consequence follows from Lemma 3.2 with `m = 4s/K + 4(Ls)^(1/2) + (s - n_J)`, giving
`8Km`. Correct.

**Cor 4.3.** The independence criterion is right: no nonzero `Q`-row-space vector supported on `<= 2` coordinates.
It needs `rank B <= s - 2`. Constants: `32 (1/ln 2)^(1/2) = 38.44 <= 38.5`, `4 (1/ln 2)^(1/2) = 4.80 <= 4.9`,
`40.4/pi = 12.86 <= 12.9`, `39 (1.26)^(1/2) = 43.78 < 44 < 60.72`. Correct as stated.

**What Cor 4.3 does not show; (F2) is overstated.** The formal statements are correct. The headline (F2), "the
short-arc condition never forces the angles to concentrate", and §8.4, "a restricted growth theorem in the window
cannot come from the large sieve", are not proved. Three gaps:

1. *Unpinned profiles only.* Cor 4.3 needs every pair of coordinates independent with `D_j = infinity`, hence
   `rank B <= s - 2`. The write-up itself says `rank B = s` is "typical when `s <= M - 1`". Near-extremal clusters
   have `M/4 (1 + o(1))^(-1) < s <= (2/eps + o(1)) M` (I3). Of the four records, two are pinned and one has every
   pair dependent (§7.2, reproduced below). So the proved case does not contain the typical case.
2. *In the pinned case concentration is forced.* By Remark 4.4, when `rank B = s` and `Delta < c(B)`, every
   `phi in Real(Delta)` lies within `||B^+|| Delta` of the finite group `H`. Each `phi_j` then sits within that
   distance of `(pi/2) Z/D_j`. With `D = lcm D_j`, all `s` angles lie in `D` arcs of width `2 ||B^+|| Delta`, and
   Prop 3.3(1) at `eta -> 0` then forbids `s > D (2.53 Y^(1/2) + 1.58)`. So in the pinned case the short-arc
   condition **does** force concentration. Whether this yields anything depends on how large `D` is forced to be
   (gaps) and on the profile-dependent threshold `c(B)`. That is exactly the territory the owner's constraint
   forbids treating uniformly in `R`. Remark 4.4 compares magnitudes: gaps force each `D_j >~ W/(2 log p_j)`, while
   the large sieve asks for `>= s/(2.53 Y^(1/2))` distinct points. This is **not** an implication: large `D_j` for
   each coordinate does not force the `s` angles onto many distinct points, because the cosets of different
   coordinates can overlap. So "in the pinned case the large sieve adds nothing to the monomial gaps" (§4.4, §8.1 F2)
   is unproved. It is correctly flagged "not used below", but the summary and §8.1 rely on it.
3. *Realisation alone vs realisation plus gaps.* `H` preserves `Real(Delta)` but **not** the gap conditions (I5). A
   point `phi + h` built from a gap-admissible `phi` generically violates gaps. So Cor 4.3 proves only that the large
   sieve (within factor 3) cannot contradict the realisation condition **alone**. A growth argument would use
   realisation, gaps and the large sieve together. Ruling that out needs `Real(Delta) cap Gaps cap Spread` to be
   non-empty at `W = O(M log M)` with window primes. That is the open "spreading proviso", i.e. the existence of
   spread window fakes (Cor 5.4's hypothesis). The write-up states the proviso in §8.2 but drops it in §8.1 (F2) and
   in §8.4.

Corrected form of (F2). *For profiles with all coordinates pairwise independent in `Cok` (in particular
`rank B <= s - 2`), every non-empty `Real(Delta)` contains points whose large-sieve constants are within a factor 3
of the genuine prime angles' bound. Hence the large sieve (with constants relaxed by a fixed factor), combined with
the realisation condition alone, cannot exclude such a profile. Nothing is proved for pinned or partially pinned
profiles, nor for the large sieve combined with the monomial gaps.*

## 4. Props 5.2-5.3 / Cor 5.4 (Theorem B extension, conditional cap)

**Prop 5.2.** Re-derived against `round10/analytic.md` Def 3.1, Lemma 3.2 and Def 3.3 (R1).
* Minkowski step: `D_k` is the grafted minus the genuine contribution of the `2 s_X` changed elements. Bounding
  each part by its own large-sieve norm needs no interaction with the unchanged primes. Correct.
* Case `K <= 2^(L-1)`: (i) applied to `psi` and to its reflection `pi/2 - psi` gives
  `m = 8s/K + 8(Ls)^(1/2)`, and Lemma 3.2 gives `U^2 <= 8Km = 64(s + K(Ls)^(1/2))`. Then
  `t = 256(s + K(Ls)^(1/2))/Lambda <= (256/c)(s/P_0 + (Ls/P_0)^(1/2))`, using `Lambda >= c(X + KX^(1/2))` and
  `X >= P_0`. Correct.
* Case `K > 2^(L-1)`: Montgomery-Vaughan with spacing `2^(-L)` gives `2K + 2^L < 4K`, and `t <= 16/(c P_0^(1/2))`.
  Correct.
* Hypothesis (ii) is automatic for any gap-admissible fake once `(pi/2) 2^(-L) <= 1/P`: the gaps (I5) at
  `v = e_i - e_j`, `e_i + e_j` and `2 e_j` give separations `>= arcsin((p_i p_j)^(-1/2))` and `arcsin(1/p_j)`. For the
  genuine angles it is the lattice cross-gap. So only (i) is a real hypothesis.
* *Quantitative slip.* The summary table (§0) says spread grafts satisfy the robust large sieve "up to
  `1 + O(M^(-eta/2)(log M)^(1/2))`". The proved factor is `(1 + t^(1/2))^2 = 1 + O(t^(1/2))`, i.e.
  `1 + O(M^(-eta/4)(log M)^(1/4))`. It still tends to 1. Minor.

**Prop 5.3.** Correct: `analytic.md` (R2), (R3) with `M^A` replaced by `P_0`. The change is
`<= 4 s ||f|| log X_1`, resp. `<= 2 s_J`, against `c X_1^gamma`, resp. `c x_B^gamma`, so it tends to 0 iff
`gamma (1 + eta) > 1` when `s = Theta(M)` and `P_0 = M^(1+eta)`.

**Scope of (F3); a major overstatement in the narrative.** Inside the window `eta < 1 - delta < 1`, so Prop 5.3
needs `gamma > 1/(1 + eta) > 1/(2 - delta) > 1/2`. **No square-root-quality input is covered anywhere in the
window**: GRH-quality Hecke prime sums (`E ~ X^(1/2) log^2`) and box counts with error `~ x^(1/2)`. For these the
graft change `~ s log X ~ M log M` exceeds the error `X^(1/2) <= M^(1 - delta/2)`. In Theorem B's own regime, by
contrast, every fixed `gamma` is covered by taking `A > 1/gamma`. Spreading might rescue (R2) at `gamma = 1/2` for
fakes whose angles have small Weyl sums, but the write-up does not prove this (it says "spreading is not used").
So the corrected (F3) is: *the large-sieve part of Theorem B extends into the window for spread grafts (Prop 5.2,
all `c`); the prime-sum and box parts extend only for `gamma > 1/(1 + eta)`, with primes `>= M^(1+eta)`.*
"The window is not special for robust analytic inputs once fakes are spread" (§8.1 F3) and "the window would be
capped for `A_0` + robust analytic inputs" are true only for that restricted class. Square-root-quality prime sums
and box counts are a **new uncapped input class inside the window**, conditional or not. This should go into the
§5.5 bookkeeping.

**Cor 5.4.** Correct as a conditional statement under its hypotheses: spread gap-admissible models, primes in
`[M^(1+eta), M^(2-delta)]`, `s = O(M)`, and `(pi/2)2^(-L) <= M^(-2)`. `s = O(M)` is in fact automatic from
`W <= kappa M log M` and `P_0 >= M^(1+eta)`, since `s (1 + eta) log M <= W`. Two further limits:
* it says nothing about window clusters whose primes lie in `[c M log M, M^(1+eta)]`, the lower part of the window
  of Cor 2.3;
* hypothesis "the genuine angles of `S` are `L`-spread" (Caveat 3) rests on an **unverified** upper-bound sieve for
  Gaussian primes in polar boxes of angle `>= P^(-1/2)` with short norm ranges. It is an extra assumption on top of
  fake existence, not just "a construction issue".

The cap is conditional on an open existence question, and the write-up says so (Caveat 1). Fine.

## 5. Props 6.1-6.2 / Heuristic 6.3 (monomial gaps)

**Prop 6.1(1)-(3).** Correct.
* Every `w(v) >= log 5`. The two excluded arcs are disjoint because `arcsin(5^(-1/2)) + arcsin(10^(-1/2)) =
  arctan(1/2) + arctan(1/3) = pi/4` (checked to machine precision).
* `(4/pi)(1 + 2^(-1/2)) = 2.17356`.
* The product formula `prod_j (1 + 2/(p_j^(1/2) - 1)) - 1` and `log(1 + x) >= x/(1 + x)` are right.
* Item 3 is right for Haar measure on `T^s` (independent coordinates).
* My independent enumeration (`gaps.py`) reproduces `E` against `1.087 (prod - 1)` for every row (§7 below).

**Prop 6.1(4): two problems.**
1. *Unproved transfer to fakes (gap, probably repairable).* Item 4 uses "`s > M/5` for `M >= M_1(C)` by
   Corollary 2.2". Cor 2.2 rests on `stability.md` Thm 3.3, which is proved for **actual** clusters. The record
   transfers it to models only in the class `A_1^>(Q_0)` with residue data (`stability.md` Prop 7.1, whose input 4
   is the residue pigeonhole). Item 4 is about constructing `A_0` fakes, and for `A_0` no transfer is proved. Without
   it, a fake with `s = o(M)` primes and large exponents is not excluded, and for such a fake the first-moment
   obstruction need not apply. Most likely the first form of (3.0), which needs only `F >= -O(s^2)` from the pair
   chords, gives `s >= (M/4)(1 - o(1))` for `A_0` with a worse constant. But this has to be proved or cited, not
   assumed.
2. *False side claim.* "All refereed fakes instead have `sum_j p_j^(-1/2) <= 1/32`." This holds for `lw.md` 6.1
   (`s < 2M` primes `>= (64M)^2`). It fails for the refereed fakes of `round5/construct.md` Thm 2: `r = 9M - 8`
   primes `>= P_* >= 21 r^2`, so `sum <= r/(21 r^2)^(1/2) = 21^(-1/2) = 0.218`, and this bound is approached at
   `P = P_*`. Correct statement: *all refereed fakes have `sum_j p_j^(-1/2) <= 0.22`*. The qualitative point
   (`O(1)` against `>> 1` in the window) is unaffected.

**Scope of (F4).** Prop 6.1 is about Haar measure on `T^s`. A fake lives on `Real(Delta)`, a neighbourhood of a
coset of the identity component of `H`. There, characters in the `Q`-span of the rows of `B` are deterministic, and
the single-prime events `v = e_j` are random only when `D_j = infinity`. So "the union bound needs
`sum_j p_j^(-1/2) = O(1)`" is right with the sum restricted to the free coordinates. For pinned or partly pinned
profiles the first-moment framing does not apply as stated. Also, the threshold `P >~ r^2 log r` of the known
equal-weight constructions has a second, purely counting cause the write-up does not mention: a window
`[P, P e^(1/(4r))]` holds `~ P/(8 r log P)` split primes, which is `>= r` only if `P >~ 8 r^2 log P`. So a window
fake cannot put its primes in one equal-weight short window at all, independently of the gaps. That supports
the write-up's thesis but should be stated.

**Prop 6.2.** A correct restatement of `analytic.md` Thm D (`Norm beta_j <= p_j`; the shadow has norm
`N'' <= N` and inherits `P_+(N'') <= max p_j`, so it stays in the window). Fine.

**Heuristic 6.3.** Correctly labelled a heuristic. Two cautions on its reading.
* The evidence (§7.3) stops at `E <= 10` (`sum p^(-1/2) <= 1.2`) on the **full torus**. The window needs
  `sum p^(-1/2) ~ M^(1 - a/2) -> infinity`, i.e. `E ~ exp(2 M^(1 - a/2))`, on the **realisation set**. That
  extrapolation is doubly exponential.
* Positivity of `exp(-kappa E)` in a Poisson model is not evidence of existence on a fixed low-dimensional
  subtorus at such intensities. "Window fakes plausibly exist" is an opinion, not a reading of the data.

## 6. Prop 8.1 and §8 (reduction)

**Prop 8.1(1)** is correct and close to tautological. A normalised actual window cluster is an angle model with
the same primes, `W' <= W` and chord constant `<= C`, so `W >= W_min^win(M)`. The case split at
`M = (log R)^(1/2)` is right.

**Prop 8.1(2)** is Cor 5.4. It inherits Cor 5.4's hypotheses: spread models, primes `>= M^(1+eta)`, and
`gamma > 1/(1+eta)` for prime sums and boxes.

**§8.2's "equivalence" is overstated.** "For every input class considered in this round, a restricted growth theorem
in the window is equivalent (up to the spreading proviso) to the non-existence of gap-admissible angle models with
small primes at `W = O(M log M)`" is right only for `A_0` alone, where it is the definition of a fake. For `A_0` plus
analytic inputs the converse direction (existence of models implies no proof) is proved only:
* with the spreading proviso;
* with the genuine-angle spreading assumption of Caveat 3;
* for prime-supported robust inputs with `gamma > 1/(1 + eta)` (not GRH-quality);
* for primes in `[M^(1+eta), M^(2-delta)]`, not in the lower part `[c M log M, M^(1+eta)]` of the window.

"The large sieve does not enter it" holds only for the large sieve with relaxed constants and with the spreading
proviso.

**§8.4.** "A restricted growth theorem in the window cannot come from the large sieve: (F1) and (F2) are proved" is
**not** established. See §3 above: (F2) is proved only for unpinned profiles and only against realisation without
gaps. See §2: sharp-constant sieves are excluded. The stated "`< 2%`" probability is an opinion. The proposed next
step (a sequential construction of spread window fakes) is reasonable. But by §5 above, any such construction cannot
use a single equal-weight short window, because there are not enough primes in it.

## 7. Evidence: numbers recomputed independently

All scripts are in `referee_prove_checks/`. None of them reuses the author's code.

**7.1 Smooth-norm scan** (`rscan.c`: my own DFS over split-only `B`-smooth `n`, long-double angles, `qsort`).
* *Counts.* `B = 200, X = 10^18`: 40,715,511 norms (plus `n = 1`) and 18,613,663,791 point classes. `B = 113`:
  4,656,319 norms. Both match the write-up exactly.
* *Table.* My run to `10^15` (`scan_B200_1e15.txt`) reproduces **every** entry of the §7.1 table for the decades
  `10^6, 10^8, ..., 10^14`, `k = 7..11`, in both the "all" and "`P_+ <= k^2`" columns, to the printed digits. The
  record constants `5.2699` (`k = 7`), `7.8804` (`k = 9`) and `8.2393` (`k = 10`) reappear.
* *The 1e15-1e16 rows.* My `X = 10^17` run (`scan_B200_1e17.txt`; 21,615,255 norms, 7.40e9 classes) reproduces the `10^15` and `10^16` rows exactly. I did not re-run the `10^17` decade, which needs `X = 10^18`.
* *Reading 1 is false as worded.* "For every `k >= 7` the best `C_k` per decade **increases** with `n` beyond
  `10^9`" is contradicted by the write-up's own table, which my run confirms:
  * `k = 7`: `7.12` (1e10), then `6.61` (1e11);
  * `k = 10`: `12.97` (1e11), then `11.68` (1e12);
  * `k = 11`: `16.17` (1e9), then `16.14` (1e12), which is below the `10^9` value.

  Correct wording: *for `k = 7..10` no decade beyond `10^9` improves on the `10^9` value; for `k = 11` the `10^12`
  decade is marginally lower (`16.14 < 16.17`); the trend is upward but not monotone.*
* *Reading 2 ("the window is populated") holds only at `delta = 0`.* I added a column `P_+(n) <= k^(7/4)`, i.e.
  `delta = 1/4`. It costs much more than `P_+ <= k^2`: at `10^9`, `k = 9` it is `13.43` against `9.35` (and `7.88`
  unrestricted); at `10^13`, `k = 7` it is `39.6` against `13.9` (and `7.12`); at `10^16`, `k = 7` it is `103.4` against `32.6` (and `8.40`). Every record cluster has
  `P_+ > k^2` (ratios 1.13-3.64; §7.2), so **no record cluster lies in the window for any `delta > 0`**. At `M <= 11`
  the difference between `k^(2-delta)` and `k^2` is a bounded factor, and "the window is populated" says nothing
  about the asymptotic window. The write-up's own caveat (reading 3) is right and should be put first. (Also the
  phrase "because all scanned primes are `<= 200 < 11^2 + 80`" makes no sense: primes `137..197 > 121` are scanned.
  For `k = 11..14` the best clusters simply do not use them.)

**7.2 Record clusters** (`records.py`, `pairs2.py`; exact Gaussian arithmetic, my own `hnf` for the `D_j`, sympy for
Smith forms). All four rows reproduced exactly:
* `N'`, `s`, `P_+`, `P_+/k^2`, and `sum p^(-1/2)` (`1.2469, 1.2469, 0.9893, 1.4342`);
* `rank B`; Smith invariants `[2,2,2,2,24]`, `[2,2,2,2,24]`, `[2,2,2,2,0]`, `[2,2,2,2]`;
* `D_j = (24, 8, 6, 4, 24)` for the 7- and 10-point records, `inf` otherwise;
* arc constants `8.2393, 5.2699, 5.8978, 2.2584` and chord constants `8.2386, 5.2695, 5.8974, 2.2584`;
* golden-8: 10 of 10 pairs dependent; 5-point: exactly one dependent pair, `(5, 37)`.

Observation (a) of §7.2 ("record clusters sit where first-moment fakes fail") is true but weak: any cluster
using the prime 5 has `sum >= 0.447`.

**Thm 2.1 on windows** (`thm21.py`). On all 724 windows (`M >= 5` consecutive points, span `< pi`, chord constant
`<= 9`) of the four record norms, `W' < s(4W'/M + 8 log^+ C)` and `M < 4 pi(Y;4,1)(1 + 2M log^+ C/W')` hold, with
exact Gaussian gcds. 0 violations. (A first version of my script without the span filter counted wrap-around
"windows" spanning more than `pi` and reported spurious violations. Those are not arcs.)

**7.3 Gap Monte Carlo** (`gaps.py`: my own enumeration of `v` up to sign, exact `Prob(Bad_v)`, direct Monte Carlo
with early elimination, and the genuine angles checked in 40-digit arithmetic):

| primes | `#v` (mine) | `E` (mine) | `log P_good` (mine, direct) | write-up |
| --- | --- | --- | --- | --- |
| 101..157 (6), WMAX 24 | 1039 | 2.090 | `-2.350 +- 0.005` (4e5 samples) | 1039, 2.09, `-2.33`/`-2.35` |
| 53, 61, 73, WMAX 22 | 115 | 1.259 | `-1.465 +- 0.003` | 115, 1.26, `-1.48` |
| first 5 `>= 53` | 759 | 2.463 | `-2.914 +- 0.004` | 759, 2.46, `-2.91` |
| first 7 `>= 53` | 2381 | 4.107 | `-4.906 +- 0.008` | 2381, 4.11, `-4.96` |
| 29..97 (8), WMAX 22 | 6786 | 8.522 | `-11.34 +- 0.10` (8e6 samples, 95 hits) | 6786, 8.52, direct `-10.8 +- 0.7`, SMC `-11.26` |

All agree within the stated errors. The genuine angles pass every enumerated gap, with ratio exactly `1.000000` at the
equality cases. I did not rerun the SMC rows for 9-12 primes. They are evidence only, and nothing depends on them.

**Constants** (`consts_ls.py`): every numerical constant in §§2-6 checked (output in `consts_ls_out.txt`).

## 8. Owner constraints

* *Uniform exceptional sets from fixed-`S` Subspace Theorems:* none used. OK.
* *Gaussian gcd or denominator costs:* Thm 2.1 works with `N'` after exact division by `G`. The passage from
  `P_+(N)` to `P_+(N')` runs in the safe direction. Prop 6.2 keeps `N''`. My `records.py` divides by the exact
  Gaussian gcd (it is a unit for all four records). OK. (The unproved transfer of Thm 3.3 to `A_0` fakes, §5
  above, is not a gcd issue: models have no content.)
* *Drifting common phase replaced by zero:* `theta_0` is existential throughout, and Lemma 4.1(a) moves it explicitly.
  OK.
* *Template-dependent threshold treated as uniform in `R`:* the formal statements comply (`M_1(C)`, `R_0(eps, C)`;
  Remark 4.4's `c(B)` is flagged). **But** the narrative (F2) and §8.1 use Remark 4.4 to dispose of the pinned case
  ("the monomial gap already forces more than the large sieve does"), and that claim lives only below the
  profile-dependent threshold `Delta < c(B)`. Using it as part of a "proved reason" in §8.1 is a soft violation in
  spirit: a template-dependent comparison used to settle a uniform question. It must be removed from the list of
  proved reasons (see §3 above).

No hard violation in any formal statement.

## 9. Corrected statement

*Status: partial. No growth theorem in the window `P_+(N') <= M^(2-delta)`; the best bound there is still
`(2/3 + eps) log R/loglog R`.*

**Proved.**
1. **Thm 2.1, Cor 2.2-2.3** (a corollary of `stability.md` Thm 3.3(3)). Every cluster with `M >= 5` points
   satisfies `M < 4 omega(N')(1 + 2M log^+ C/W') <= 4 pi(P_+(N'); 4, 1)(1 + 2M log^+ C/W')`. Hence, for fixed `C`,
   clusters with `P_+(N') = o(log R)` satisfy `M = o(log R/loglog R)`, uniformly in the arc's position. This is a
   restricted growth theorem **below** the window. Near-extremal clusters have
   `P_+(N') > (eps/12.56) log R >= c_eps M log M`.
2. **Lemma 3.2.** For any finite set of angles, the optimal dual large-sieve constant at frequencies `|k| <= K` lies
   between `(2K + 1) m''/sqrt 2` and `32 K m''`, where `m''` is the maximal count in an arc of length `1/(8K)`. For a
   cluster's prime angles (indeed for any non-associate Gaussian integers of norm `<= Y`) the counts are at most the
   lattice count `eta Y/2 + 2.53 Y^(1/2) + 1.58`, which in the window is non-trivial by a factor `>> M^(delta/2)`.
3. **Prop 4.2 / Cor 4.3.** `Real(Delta)` is invariant under `H = ker(B mod pi/2)`. If all coordinates are pairwise
   independent in `Cok` (which needs `rank B <= s - 2`), every non-empty `Real(Delta)` contains points whose
   large-sieve constants and arc counts are within a factor 3 of the genuine prime-angle bounds. So the large sieve,
   with constants relaxed by a fixed factor, **together with the realisation condition alone**, cannot exclude such
   profiles. Nothing is proved for pinned or partly pinned profiles, nor for the large sieve combined with the
   monomial gaps.
4. **Props 5.2-5.3, Cor 5.4 (conditional).** Grafts whose angles and whose primes' genuine angles are spread satisfy
   every `(c, gamma)`-robust prime-supported large sieve up to `1 + O(M^(-eta/4)(log M)^(1/4))`, even with
   `s >> X^(1/2)`. They satisfy robust prime sums and box counts only for `gamma > 1/(1 + eta)`. Hence, **if** spread
   gap-admissible angle models exist with `W = O(M log M)` and all primes in `[M^(1+eta), M^(2-delta)]` (and the
   genuine angles of the chosen primes are spread, an unverified assumption), then `A_0` plus these restricted
   robust inputs cannot prove growth for such clusters. Square-root-quality (`gamma = 1/2`) prime sums and box
   counts are **not** covered anywhere in the window, nor are primes below `M^(1+eta)`.
5. **Prop 6.1.** For Haar-random `phi` the expected number of violated monomial gaps is
   `>= 1.08 (exp(2 sum_j (p_j^(1/2) + 1)^(-1)) - 1)`, and `Prob(all single-prime gaps) <= exp(-2.17 sum_j p_j^(-1/2))`.
   So first-moment constructions on the free coordinates need `sum_j p_j^(-1/2) = O(1)`. All refereed fakes have
   `sum_j p_j^(-1/2) <= 0.22`. For a window object with `s >= cM` primes the sum is `>> M^(delta/2)`.
   (`s >= cM` is proved for actual clusters; for `A_0` fakes it still needs a transfer of `stability.md` Thm 3.3.)
   A single equal-weight short window `[P, P e^(1/(4r))]` cannot even hold `r ~ M` primes when `P <= M^(2-delta)`.
6. **Prop 6.2** (`analytic.md` Thm D). Lattice-angle realisations are shadows of actual window clusters.
7. **Prop 8.1.** For `A_0` alone, growth in the window is equivalent to `W_min^win(M)/(M log M) -> infinity` (a
   restatement of the definition of a fake). With analytic inputs only the conditional half (item 4) is available.

**Heuristic / not proved.** Remark 2.4 ("counting cannot pass exponent 1"); "in the pinned case the monomial gaps
already give more than the large sieve" (Remark 4.4); "the large sieve plays no part / a window theorem cannot come
from the large sieve"; Heuristic 6.3 ("window fakes plausibly exist").

**Evidence.**
* The smooth scan to `10^18` and the record table are reproduced exactly (re-run through `X = 10^17`). Best
  constants per decade beyond `10^9` trend upward but not monotonically, and never beat the records for
  `k = 7, 9, 10`.
* Small clusters reach `P_+ <= k^2` (`delta = 0`). No record cluster has `P_+ <= k^2`. With `P_+ <= k^(7/4)` the
  best constants are worse by factors from 1 (small decades) up to about 12 (`k = 7`, `10^16`).
* The gap-admissible Haar measure is about `exp(-1.1 to 1.3 E)` for `E <= 10` on the full torus only.

**New open item for the barrier map (§5.5).** Square-root-quality Hecke prime sums and polar-box counts
(`gamma <= 1/(1 + eta)`, in particular GRH-quality) are not capped inside the window by any proved or conditional
result. Neither are sharp-constant large sieves, nor the pinned case of (F2).
