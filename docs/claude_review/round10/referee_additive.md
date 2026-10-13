# Referee report on round10/additive.md

Adversarial referee, 2026-10-13. Scripts: `referee_additive_checks/`. They run from that directory with
`python3` and need `mpmath`, `sympy` and `python-flint`. Outputs are the `*_out.txt` files there. Labels as in
the write-up: **proved**, **exact** (integer arithmetic, or mpmath at 60-200 digits for angles and logs),
**evidence**, **heuristic**.

## 0. Verdict

**The verdict "partial, no growth, no uniform bound" stands.** Every theorem I re-derived is correct as a
statement about its own hypotheses: Thms 2.1-2.3, 2.4 (in its literal box-height form), 2.6, Prop 2.7,
Prop 2.8 (literal form), Prop 3.1, Thm 3.2, Cor 3.3, Lemmas 4.1-4.2, Prop 4.3 and Remark 4.4. Independent
exact checks (4,602 runs, 1,790 of them with `M >= 4`) find no counterexample. Prop 3.1 and Thm 3.2(b)
are sharper on data than the write-up reports: the slacks are `-0.010` and `-0.025`.

**The headline "fake-violated input" does not survive as stated.** It has three defects.

1. **(major, logic) Generic, not universal.** `AG(M^A)` is shown false for *almost every* fake of a
   positive-measure family. The barrier is an existence statement: one fake that satisfies an argument's
   inputs at `W = O(M log M)` caps the argument. A.e. violation does not exclude such a fake, so it does
   not escape the barrier. The universal statement is the open Lemma (AG) of §7. That lemma is
   **strictly stronger than the growth theorem**: actual clusters are angle models satisfying `AG(infinity)`,
   so its restriction to them *is* growth at `c_0`. It also contains the norm-form problem for every
   imaginary quadratic order. So `AG` is an input that the *existing* barrier constructions violate for
   generic positions. It is not shown to be fake-violated in the sense the task requires, and using it is
   not shown to be easier than proving growth.
2. **(major, quantitative) The heights are wrong, by an additive `1` and by a factor 2 in `c_0`.**
   * Thm 2.4(1) bounds the *box* height `H_1 = M^(c_0/2+o(1))`. The `AG` height (`l_1`) is `2(M-1)H_1 = M^(1+c_0/2+o(1))`.
     The summary, Cor 2.5, Prop 2.8's reading, Heuristic 2.9 and the §6 table all say `M^(c_0/2+o(1))`.
   * The pigeonhole also ignores curvature. The values `sum c'_x (w_x - w_1)` lie in a rectangle of radial
     width `~ H M C^2` and tangential width `~ H M L`, not in a disc of radius `H M L`. With this correction,
     **every relation-free configuration on the arc violates `AG` at `l_1` height `M^(1+c_0/4+o(1))`**
     (referee Thm R2, proved). This is also the typical height: the heuristic, corrected for the
     `(M-1)!` that Heuristic 2.9 omits, gives the same exponent. Exact LLL on random arcs and on explicit
     skew double core fakes confirms it (§2.2, §6.2).
   * Corrected ladder: relation-free threshold `A = 1 + c_0/4`; Prop 2.8 becomes
     `rho <= 2 + c_0 M/(4(A-1)) + O(M/log M)`; the dimension-count ceiling for imposed relations is
     `A = 1 + c_0/2` (heuristic), not `c_0`; full lattice-ness is forced at height `e^(W/4+O(M))`.
   * So the claim "violated at `M^(c_0/2+o(1))`" holds (through Thm R2) only for `c_0 >= 4`. For
     `3 <= c_0 < 4` it is unproved, and the heuristic says it is false. The qualitative picture survives,
     because `A` may be chosen after `c_0`.
3. **(minor) Cor 2.5's measure space.** Positions are `delta = A~ phi - theta_0 1` (mod `pi/2`). For a
   rank-deficient profile they fill a null subset of `R^M`, and Cor 2.5 as stated is then vacuous. The
   statement holds in `phi`-space, since `F_c o delta` is a nonzero exponential polynomial when the rows
   `a_x` are distinct. For the skew double core, rank `M`, the hypothesis holds as stated. Its admissible
   position set has positive measure but **empty interior** (referee Prop R3). This answers the write-up's
   open point 7.2(a).

**What the round adds, after correction.**
* The collapse theorems (§2.1, §2.3, §2.4) are correct and clarifying.
* Prop 3.1 is correct and nearly sharp, and capped (it uses monomial gaps only).
* Thm 3.2 is correct and nearly sharp. It is a restricted (R1-arch) in low position.
* One new referee observation (Prop R1, proved) locates `AG`'s content. The monomial gaps, which fakes
  satisfy, already give `|sum c w| >= sqrt2 e^((W - w(v(c)))/2) - o(1)` for every `c` with `v(c) != 0`.
  So `AG` adds nothing on combinations of multiplicative weight `w(v(c)) <= W + log 2 - o(1)`. Its whole
  content is on combinations with `w(v(c)) > W`. Every exact relation has `w(v(c)) >= 3W/2 - O(log(n C^3))`.
  On actual data the minimum is `1.414 W`.

## 1. What was checked and how

I re-derived every proof in §§2-4 line by line (§§2-5 below). I wrote new code; nothing is imported from
`additive_checks/`.

| script | what it does | result |
|---|---|---|
| `glib_ref.py` | Gaussian toolkit: Euclidean gcd, primes above `p = 1 mod 4`, all points of a circle, maximal windows, normalisation, profiles | library |
| `actual_checks.py` | On every maximal window of width `<= C N^(-1/4)`, normalised. Records `1176852625`, `537606725`, `235370525` at `C in {3,4,5.3,5.9,6,7,8}`; random split `N <= 10^13` (150 values, `C <= 4`) and `N <= 10^15` (400 values, `C <= 6`). Checks Prop 3.1, Thm 3.2(a)(b)(c), round-9 Thm 1.1 dichotomy, Prop 2.7(1)(2)(3) (index by HNF, independently of the minors), `rank K = M - 3` (integer kernel by exact LLL), the exact covolume identity `covol(K)^2 I^2 = sum (2 Area)^2`, Remark 4.4, and referee Prop R1(ii) on LLL-reduced relations | `actual_checks_out.txt` (2,455 runs, 1,454 with `M >= 4`) and `actual_checks_out2.txt` (4,602 runs, 1,790 with `M >= 4`); the record circles are in both. **0 failures** |
| `ag_fake.py` | Explicit skew double core fakes (round5 `exact.md` Thm 5.1; `M = 8` in classes `A_0` and `A_1(2M)`, `M = 12` in `A_0`) with random 664-bit positions. Pre-checks the Thm 5.1 gap requirement on all characters with `||c||_1 <= 6` (`M = 8`) or `<= 4` (`M = 12`). Exact LLL sweep for zero-sum `c` with `0 < |sum c w| < 1`. Compares with the author's disc pigeonhole, the referee's rectangle pigeonhole and the corrected heuristic. Part B: the same on generic arcs, `M in {6, 8}`, `log L in {20, 40, 80, 160}` | `ag_fake_out.txt` |
| `p27_qdivN.py` | Prop 2.7(2) for odd `q | N`: `q | I` iff exactly 2 residues mod `q` | `p27_qdivN_out.txt`: 13,024 pairs, 0 failures |
| `r1_check.py` | Referee Prop R1(i) on the skew core fakes, every `c` with `||c||_1 <= 6` (`M = 8`) or `<= 4` (`M = 12`) | `r1_check_out.txt`: 0 violations |

A trap that matters for anyone who runs the write-up's proposed finite test (§7.2). My first run of
`ag_fake.py` drew positions with `random.random()`, i.e. 53-bit dyadic rationals. Such positions carry
rational structure. Many `c` then make the first-order term `c . delta` vanish *exactly*, and spurious
`AG` violations appear at heights `e^5`-`e^8`, far below the true ones. That output is kept as
`ag_fake_out_53bit_artifact.txt`. All reported numbers use 664-bit random positions.

## 2. Section 2 (collapse theorems, AG)

### 2.1 Thms 2.1-2.3: correct

* **Thm 2.1.** The map `Lambda_0 -> Gamma` is onto with kernel `K(w)`. A rank-2 subgroup spanning `R^2`
  is discrete. Correct.
* **Thm 2.2.** `K(z) ⊂ K(w)` makes `z_x - z_1 -> w_x - w_1` a well-defined homomorphism on the rank-2
  lattice `Lambda`. It extends uniquely to an `R`-linear `T`, invertible since three points of a circle
  are not collinear. Two irreducible conics sharing `>= 5` points coincide, so `T` is a (possibly
  orientation-reversing) similarity. Correct.
* **Thm 2.3.** The coefficient vector `X` spans the rational kernel of an integer matrix. A second kernel
  vector would be a conic sharing `>= 5` points with the irreducible conic `X`, so it is proportional to
  `X`. The centre then solves a rational system. Correct.
* **Reading of Thm 2.3.** Correct, and in my view the most useful conceptual point of the write-up:
  existence of the relation lattice cannot tell `Z[i]` from any other imaginary quadratic order.

### 2.2 Thm 2.4 and Cor 2.5: correct as literally stated; the asymptotic readings are wrong

**(a) Literal statements.**
* Thm 2.4(1): squares of side `2/3` have diameter `0.943 < 1`, and the difference vector is zero-sum
  with `||c||_1 <= 2(M-1)H`. Correct.
* Thm 2.4(2): the fibres of `c' -> V(c')` are translates of `K' ∩ box`, with at most `(2H+1)^kappa`
  points, and points pairwise `>= 1` apart number at most `(2r+1)^2` in a disc of radius `r`. Correct.
* Cor 2.5: at a zero of `F_c` one of `Re F_c`, `Im F_c` has nonzero gradient, so the zero set is locally a
  `C^1` hypersurface. A countable union of null sets is null. Correct.

**(b) Slip 1: box height versus `AG` height.** `H_1` in Thm 2.4(1) is the *box* height, the bound on each
`|c'_x|`. The `AG` height is `||c||_1 <= 2(M-1) H_1`. With `H_1 = M^(c_0/2+o(1))` this is
`M^(1 + c_0/2 + o(1))`; the factor `2(M-1)` is a full power of `M`, not `M^(o(1))`. The same slip affects:
* the summary ("violates AG at height `M^(c_0/2+o(1))`");
* the "Fake-violated input" paragraph ("`A > c_0/2`");
* Prop 2.8's reading: `AG(M^A)` means box height `M^(A-1)/2`, so the bound is
  `rho <= 2 + c_0 M/(2(A-1))`, not `c_0 M/(2A)`;
* Heuristic 2.9 and the §6 table.

**(c) Slip 2: curvature is discarded.** Put `s_x = delta_x - Delta/2 in [-Delta/2, Delta/2]` and rotate by
the arc's mid-angle. Then `w_x - w_1 = R(e^(i s_x) - e^(i s_1))` has
* radial part `R(cos s_x - cos s_1)`, of absolute value `<= R Delta^2/8 = C^2/8`;
* tangential part `R(sin s_x - sin s_1)`, of absolute value `<= R Delta = L`.

The disc of radius `(M-1)HL` used in Thm 2.4 and Prop 2.8 can therefore be replaced by a rectangle of
sides `(M-1)H C^2/4` and `2(M-1)HL`.

> **Theorem R2 (referee; proved).** Let `w_x = R e^(i(theta_0 + delta_x))`, `delta_x in [0, Delta]`,
> `M >= 4`, `L = R Delta >= 1`, `C^2 = R Delta^2`. Put `a = (M-1)H C^2/4` and `b = 2(M-1)HL`.
> 1. If `(H+1)^(M-1) > (3a/2 + 1)(3b/2 + 1)`, some nonzero `c in Lambda_0` with `||c||_1 <= 2(M-1)H`
>    has `|sum c_x w_x| < 1`.
> 2. If `w` satisfies `AG(2(M-1)H)`, then
>    `(H+1)^(rho-2) <= 2^(M-1-rho) (4/pi)(a/H + 1)(b/H + 1)`. In particular
>    `(H+1)^(rho-2) <= 2^(M+2) M^2 L max(C^2, 1)`.
>
> For `L = C N^(1/4)`, `W <= c_0 M log M` and `C` fixed:
> * the least admissible `H` in (1) is `M^(c_0/4 + o(1))`, so `AG` is violated at `l_1` height
>   `M^(1 + c_0/4 + o(1))` by every relation-free configuration;
> * `AG(M^A)` forces `rho <= 2 + c_0 M/(4(A-1)) + O(M/log M)`;
> * `AG` at `l_1` height `e^(W/4 + O(M))` forces a lattice coset.

*Proof.*
1. Cover the rectangle by squares of side `2/3`. This needs at most `(3a/2 + 1)(3b/2 + 1)` squares, and
   then the proof of Thm 2.4(1) applies.
2. The proof of Thm 2.4(2) gives at least `(H+1)^(M-1)/(2H+1)^kappa` distinct values, pairwise `>= 1`
   apart, in an `a x b` rectangle. The disjoint discs of radius `1/2` around them lie in an
   `(a+1) x (b+1)` rectangle, so there are at most `(4/pi)(a+1)(b+1)` values. Then use
   `(2H+1)^kappa <= 2^kappa (H+1)^kappa` and `kappa = M - 1 - rho`.
3. Asymptotics: `log(a b/H^2) = log L + 2 log C + O(log M) = W/4 + O(log M + log C)`, and
   `W/(4(M-3)) = (c_0/4 + o(1)) log M`.

QED

This is sharp in order, not just an upper bound. The expected number of `c in Lambda_0` with
`||c||_1 <= h` and `|sum c w| < 1`, for uniform random positions, is about
`(2h)^(M-1)/(M-1)! * M^2/(h^2 L C^2)` (heuristic). Heuristic 2.9 replaces `#{||c||_1 <= h} ~ (2h)^(M-1)/(M-1)!`
by `h^(M-1)`. The expectation reaches 1 at `h = M^(1 + c_0/4 + o(1))/(2e)`. **Evidence** (`ag_fake.py`
Part B; exact LLL, 664-bit random positions, 4 samples per row, natural logs):

| `M` | `log L` | author's disc bound | rectangle bound (Thm R2) | corrected heuristic | least violation found (min / median / max) |
|---|---|---|---|---|---|
| 6 | 40 | 30.77 | 16.75 | 11.63 | 11.86 / 12.58 / 12.91 |
| 6 | 80 | 57.44 | 30.08 | 24.96 | 25.13 / 25.41 / 25.78 |
| 6 | 160 | 110.77 | 56.75 | 51.63 | 51.72 / 52.40 / 52.59 |
| 8 | 40 | 19.86 | 11.44 | 7.34 | 7.55 / 7.85 / 7.96 |
| 8 | 80 | 35.86 | 19.44 | 15.34 | 15.69 / 15.95 / 16.36 |
| 8 | 160 | 67.86 | 35.44 | 31.34 | 31.41 / 31.77 / 32.28 |

"Found" is an upper bound for the true least height, since LLL can miss shorter vectors. Least-squares slopes of the median `log h` against `log L`, over `log L in {20, 40, 80, 160}` (rows for 20 are in the output file), are `0.332` (`M = 6`) and `0.198` (`M = 8`). The curvature-aware
prediction is `1/(M-3)`, i.e. `1/3` and `1/5`. The disc pigeonhole has slope `2/(M-3)`.

**Consequences for the write-up.**
* "Every configuration either has an exact relation or violates `AG` at height `M^(c_0/2+o(1))`" is true
  for `c_0 >= 4`, via Thm R2, since then `1 + c_0/4 <= c_0/2`. The write-up's own proof does not give it
  for any `c_0`.
* For `3 <= c_0 < 4` (and `two_thirds` gives `c_0 >= 3`) it is unproved, and the corrected heuristic says
  that typical relation-free configurations satisfy `AG(M^A)` for `c_0/2 < A < 1 + c_0/4`.
* The "relation-free" thresholds of the ladder should read `1 + c_0/4`, both proved (upper) and typical
  (heuristic, lower).

**(d) Heuristic 2.9, imposed relations.** Plug `rho >= (M-1)/2` into Thm R2(2) instead of Prop 2.8. The
attainable `l_1` height is `M^(1 + c_0/2 + o(1))`, not `M^(c_0 + o(1))`. So the heuristic ceiling for
fakes with independently imposed relations is `A = 1 + c_0/2`.

**(e) Cor 2.5, measure space.**
* For a profile of rank `< M`, the realisable positions `delta = A~ phi - theta_0 1 - (pi/2)k` form a null
  subset of `R^M`, and Cor 2.5 says nothing.
* *Corrected statement (proved).* For Lebesgue-a.e. `phi`, the points have no exact relation. Indeed
  `F_c(delta(phi)) = sum_x c_x e^(-i theta_0 - i pi k_x/2) e^(i <2a_x - e, phi>)` is an exponential
  polynomial in `phi` with pairwise distinct frequencies, because distinct points of an arc shorter than
  `pi/2` have distinct rows `a_x`. It is therefore not identically zero, and its zero set is a proper real-analytic
  subset, hence null.
* So the corollary holds for every family that has positive measure in `phi`-space.

### 2.3 Thm 2.6 and Prop 2.7: correct

* **Thm 2.6(1).** Nonzero elements of `N Z[i]` have modulus `>= N`. The four-point relations have
  `l_1` norm `<= 4 max|2 Area| <= 2 C^3 N^(1/4)`. Correct.
* **Thm 2.6(2).** The norm lies in `(0, 2N)` and is `0 mod N`. Correct.
* **Prop 2.7.**
  * (1) Gcd of the maximal minors. The parity holds because `N` is odd.
  * (2) A line meets the nondegenerate conic in at most 2 points.
  * (3) `K` maps onto the kernel of `Z^(M-3) -> Lambda/Lambda_123`, which is surjective.

  All three are correct. Checked exactly on 1,790 clusters with `M >= 4` (the write-up checked 468). For
  (1) I computed `I` by HNF and compared it with the gcd of all minors. For (3) I used an independent
  integer kernel.
* *Remark.* Item (2) also holds for odd `q | N`, where the reduction is the pair of lines `pi | z`,
  `conj(pi) | z`. After normalisation both lines carry residues, so `q | I` iff there are exactly 2
  residues. The restriction "`q` not dividing `N`" is harmless. Exact check (`p27_qdivN.py`): 13,024
  pairs (cluster, `q | N`), 180 with `q | I`, no exception.

**New exact identity (referee).** Let `B` be the `3 x M` matrix with rows `1`, `Re z`, `Im z`. The kernel `K`
is saturated, and `det(B B^T) = sum_triples (2 Area)^2` by Cauchy-Binet. The gcd of the `3 x 3` minors is `I`.
Hence

```text
covol(K)^2 * I^2 = sum over triples {x,y,w} of (2 Area_xyw)^2      (exact; 1,790/1,790 clusters).
```

So `covol(K) <= C(M,3)^(1/2) C^3 N^(1/4)/(2I)`. By Minkowski's second theorem, `K` has `M - 3` independent
relations whose `l_2` heights have geometric mean at most `sqrt(M-3) covol(K)^(1/(M-3)) = M^(c_0/4 + 1/2 + o(1))`. This corrects §4's "geometric mean
`M^(c_0/2 + 1/2 + o(1))`". Round-9 `L1L2.md` Prop 5.2 bounded both rows of `Phi` by the chord length,
which loses the same curvature factor; its bound is valid but not sharp.

## 3. Section 3 (Prop 3.1, Thm 3.2, Cor 3.3)

**Prop 3.1: correct, capped, nearly sharp.**

Re-derivation:
* Step 1 is the per-column inequality `|lo + hi - e| <= e - (hi - lo)`. Its two cases are `2 lo >= 0` and
  `2 hi <= 2e`.
* Step 2 adds the two realisation equations.
* Step 3 uses the `(pi/4) Z` form of Lemma 0.1 with `dist <= D = C e^(-W/4)`, which gives
  `w(v) >= W/2 - 2 log C - log 2`. Together with `w(b_P) <= W/2 - mu_P` this gives
  `mu_P + mu_P' <= W/2 + 2 log C + log 2`.
* Step 4: `a_y = a_y'` forces `delta_y = delta_y'`.

All correct. Exact check (`actual_checks.py`, 4,602 runs): no failure, and 3,062 pairs above the threshold. The largest value of
`mu - (W/4 + log C' + (1/2) log 2)` over pairs outside the class `b*` is `-0.0100`, so the constant is nearly
optimal on data. The write-up correctly labels Prop 3.1 as class `A_0` and capped. It uses only the
monomial gaps, which the fakes satisfy by construction.

**Thm 3.2: correct.**
* (a) is `outside.md` Prop C.2 with `L^2/R = C^2`. When `eta > 0` the bound is
  `|v| C N^(1/4) eta + |v| C^2/2 + 1 <= C^2 h(v) + 1`.
* (b): the angle between the chord and `i v*` is at most `eta + D`. The cross term `Im(u conj(i v))` is a
  nonzero integer. `|u_xy| <= sqrt2 |g_xy|` because `B_xy` is a rational integer. The `|B|/t` bound uses
  `nu_s >= 2`.
* (c) is Minkowski's theorem in a box of area 4, with `eta <= alpha <= (pi/2) sin alpha`. The step
  `|v| <= Q + 1` needs `Q >= 1`, which holds for `N^(1/4) >= 2C/pi`. This is harmless, but the
  condition should be stated.

Exact check, 4,602 runs:
* (a): max of `(M-2)/(C'^2 h*)` is `0.45`;
* (b): max slack of `mu - 2 log(sqrt2 C' h*)` is `-0.0248` (the write-up reports `-0.061`), and the
  exceptional pairs always form a matching, empty when `eta(v*) > 0`;
* (c): max of `h*/((2pi/C')^(1/2) N^(1/8) + 1)` is `0.43`;
* every pair satisfies the round-9 Thm 1.1 dichotomy (`q/2` or `q/4` a square, never both).

**Cor 3.3: correct.** One remark: (R1-arch) is a statement about *every* near-extremal cluster. Cor 3.3
gives it only on clusters in low position, and only off one matching. The write-up says "restricted",
which is right. Its remark that fakes in generic position meet every conclusion of Thm 3.2 is also right,
so Thm 3.2 is not a fake-violated mechanism. (It uses an input that is false for fakes, but its conclusions
are not.)

## 4. Section 4 (Newton polygon lemmas)

* **Lemma 4.1.** `v_pi(c_x) = v_p(c_x)` for a rational integer, and the minimum of a vanishing sum is attained twice. Correct.
* **Lemma 4.2.** For `x in S_m`, `c_x z_x/pi^m = c'_x conj(pi)^(v_x) (z_x pi^(-a_x))`, and the other terms
  vanish mod `pi`. Correct.
* **Prop 4.3.** `M - 3` independent relations of sup-norm `< p_j` all vanish at `x*`. They then lie in a
  rank-`(M-4)` lattice, a contradiction. Correct.
* **Remark 4.4.** Every 4-subset of the column must show `{a,a,b,b}` or `{a,a,a,a}`. For `M >= 5` the
  column is then constant, so `p_j <= max |2 Area| <= C^3 N^(1/4)/2`. Correct; exact check on 648
  clusters with `M >= 5` (from the record circles).
* **The "geometric mean" sentence of §4** ("M - 3 independent relations whose heights have geometric mean
  `M^(c_0/2 + 1/2 + o(1))`") should be an upper bound, and the exponent should be `c_0/4 + 1/2` in `l_2` (§2.3
  above). The heuristic paragraph after it ("random level data almost never admit even one relation") is
  correctly labelled heuristic. Its count `exp(M log H - 2W)` does not depend on that exponent.
* **Referee addition (proved; Prop R1(ii), §6.1).** Every exact relation `c` with `v(c) != 0` has
  `w(v(c)) >= 3W/2 - 6 log C - 2 log ||c||_1 + 2 log(96/sqrt2)`. Exact check on 2,750 LLL-reduced relations
  of actual clusters: min slack `3.70`, min `w(v(c))/W = 1.414`. 80 basis relations have `v(c) = 0`, where
  the statement is void.

## 5. Section 5 (reductions to capped classes)

I checked each reduction for correctness; none is load-bearing for the verdict.
* §5.1: the parallelogram law, and the fact that two chords with the same midpoint `m != 0` coincide.
  Correct.
* §5.3: `sigma` cocycle with defect `2 Area`. Correct.
* §5.4: `K` is saturated, and three residues on a nondegenerate conic span the plane. Correct.
* §5.6: correct.
* §5.5 cites earlier refereed caps. I did not re-derive them.
* The open point §7.3 (`MOD_N`) is stated correctly. The bound `M <= d(N)` holds because a box of side
  `L < R/sqrt2` meets each principal ideal of norm `N` at most once.

## 6. The fake-violated input: does AG really separate?

### 6.1 Where AG's content lies (referee Prop R1; proved)

> **Proposition R1.** Let `w_x = R e^(i(theta_0 + delta_x))` be the points of an angle model (actual or fake)
> with `delta_x in [0, D]`, `D = C N^(-1/4)`. Let `c in Lambda_0`, `n = ||c||_1 <= pi/(2D)`, and suppose
> `v(c) := sum_x c_x a_x != 0`. Then
> 1. `|sum_x c_x w_x| >= sqrt2 e^((W - w(v(c)))/2) - n C^3 N^(-1/4)/48`;
> 2. if `sum_x c_x w_x = 0`, then `w(v(c)) >= 3W/2 - 6 log C - 2 log n + 2 log(96/sqrt2)`.

*Proof.*
1. Put `s_x = delta_x - D/2`; for zero-sum `c`, `c . s = c . delta`. By (1.1),
   `<v(c), phi> = c . delta/2 + (pi/4) c . k`.
2. Since `|c . delta/2| <= nD/4 <= pi/8`, the distance of this value to `(pi/4) Z` is `|c . delta|/2`. The
   `(pi/4) Z` form of Lemma 0.1 gives `|c . s| >= sqrt2 e^(-w/2)` with `w = w(v(c))`.
3. Up to a rotation, `|sum c_x w_x| >= R |sum c_x sin s_x| >= R(|c . s| - sum_x |c_x| |s_x|^3/6)`
   `>= sqrt2 R e^(-w/2) - n R D^3/48`. Moreover `R e^(-w/2) = e^((W - w)/2)` and `R D^3 = C^3 N^(-1/4)`.
   This proves (1).
4. If the sum vanishes, then `sqrt2 e^(-w/2) <= |c . s| <= n D^3/48`, which is (2) after taking logs.

QED

*Reading.*
* The monomial gaps (class `A_0`, satisfied by every fake) already give `AG` on every `c` with
  `w(v(c)) <= W + log 2 - o(1)`.
* So every `AG`-violating `c` (`0 < |sum c w| < 1`) has `w(v(c)) > W + log 2 - o(1)`, and every exact
  relation has `w(v(c)) >= (3/2 - o(1)) W`.
* `AG`'s new content therefore lies entirely in additive combinations whose associated monomial has weight
  above `W`. That is the range where the gap inequalities are vacuous at the scale of the arc.

*Checks.*
* **exact**, actual clusters: (2) holds on 2,750 LLL-reduced relations, min slack `3.70`, min `w/W = 1.414`.
* `r1_check.py`, skew core fakes: (1) holds for all 3,318 (`M = 8`) and 2,277 (`M = 12`) characters
  tested, including 349 and 1,551 with `w(v(c)) <= W`. The least ratio of `|sum c w|` to
  `sqrt2 e^((W-w)/2)` is `18.2`.

### 6.2 Test on the round-5 fakes and the skew double core

**All barrier constructions draw positions at random.** This covers:
* `outside.md` Thm A.3;
* `construct.md` (§5, step 3);
* `lower.md` Lemma 5.3;
* `exact.md` Thm 5.1 (the skew double core);
* `round6_sharp/sharp.md` Thms 5.1 and 6.4.

Each draws `delta` uniformly on `[0, Delta]^M` and succeeds with probability bounded below (`>= 0.586`
for Thm 5.1). Each has a rank-`M` profile, so the realisable positions are all of `[0, Delta]^M`. Rank `M` is
stated in `outside.md` Thm A.3 step 4, `construct.md` ("`rank S = M`"), `lower.md` Lemma 3.1 and `exact.md`
Thm 5.1 step 5. `sharp.md` reuses the twin-column profiles. Hence
Cor 2.5 applies **as stated**: almost every member has no exact relation. By Thm 2.4(1), or Thm R2,
almost every member violates `AG` at `l_1` height `M^(1+c_0/4+o(1))`. For these families the author's
claim holds, with the corrected height.

**Explicit instances** (`ag_fake.py` Part A; skew double core, `C = 1`, primes from `exact.md` Thm 5.1;
664-bit random positions passing the Thm 5.1 gap requirement for all characters with `||c||_1 <= 6`
(`M = 8`) or `<= 4` (`M = 12`); natural logs):

| instance | `W` | `W/(M log M)` | disc bound (author) | rectangle (Thm R2) | heuristic | least `l_1` violation found (3 samples) |
|---|---|---|---|---|---|---|
| `M = 8`, `A_0` | 284.2 | 17.1 | 32.28 | 17.65 | 13.55 | 13.96, 14.21, 14.37 |
| `M = 8`, `A_1(2M)` | 409.9 | 24.6 | 44.85 | 23.94 | 19.84 | 20.20, 20.53, 20.60 |
| `M = 12`, `A_0` | 482.3 | 16.2 | 30.66 | 17.03 | 13.63 | 13.92, 14.02, 14.10 |

In every sample:
* no exact relation appears among the candidates;
* every violating `c` has `w(v(c)) - W` of order `10^7` to `10^10`, as Prop R1 predicts.

So **the skew double core fakes, as constructed, violate `AG`**, at `l_1` height about `e^(14)` for
`M = 8` in class `A_0`. That is about `M^(6.7)`, against the write-up's `M^(c_0/2) = M^(8.5)`. At such small `M`
the asymptotic exponents are not informative: the main term of Thm R2 is `W/(4(M-3))`, not
`W/(4M)`. The finite-`M` formulas in the table are what the data match.

### 6.3 The admissible position set has empty interior (referee Prop R3; proved)

> **Proposition R3.** For any profile of rank `M` (in particular the skew double core), the set of
> positions `delta in [0, Delta]^M` satisfying Lemma 0.1 for every integral character has empty interior.
> For the barrier families it has positive measure.

*Proof.* For each `c in Lambda_0 \ {0}`, the violation set
`{delta : dist(c . delta/2 + (pi/4) c . k, (pi/4) Z) < arcsin(e^(-w(v(c))/2)/sqrt2)}` is open and contains
the hyperplane `{c . delta = 0}`; here `v(c) != 0` by rank `M`. Fix any ball `B` in the cube, centred at a point `delta_0` off the
diagonal. Choose a unit vector in `{sum = 0}` orthogonal to `delta_0`. Rational directions `c/|c|` are
dense in the unit sphere of that `(M-1)`-space, so some `c` has `|c . delta_0| < radius(B) |c|`, and its
hyperplane meets `B`. So the union of the violation sets is open and dense, and its complement has empty
interior. Positive measure is the union bound of the constructions. QED

This answers the write-up's open point 7.2(a). It also refines the picture in Heuristic 2.9.
* Fakes are a fat Cantor-type set in position space. "Generic" must mean "Lebesgue-almost every", never
  "open".
* Imposing relations means intersecting with a codimension-`2k` submanifold. Whether that intersection
  is nonempty is an additional question, which the dimension count does not address. Prop R1(2) at
  least says that the imposed relations must have weight `w(v(c)) >= 3W/2 - O(log(nC^3))`.

### 6.4 Barrier logic: generic violation is not violation by the fakes

The task asks for an input that is "provably FALSE for the fakes (or for every fake of size
`O(M log M)`)". There are two readings.
1. **"The fakes" are the existing constructions.** Then `AG(M^A)`, `A > 1 + c_0/4`, qualifies in a weak
   sense. The constructions as proved (random positions) give members that violate `AG` with
   probability 1. But the constructions are existence statements. Their good sets have positive measure
   and empty interior (Prop R3). Nothing excludes a null subset of `AG`-satisfying members. **An argument
   that uses `AG` escapes these constructions only in this weak sense, and is not known to escape them
   in any other.**
2. **Every fake of size `O(M log M)`.** This is exactly Lemma (AG), and it is open. Three facts about it
   should be in the headline.
   * *It is strictly stronger than the target.* A normalised actual cluster is an angle model with
     `AG(infinity)`, so Lemma (AG) restricted to actual clusters is the growth statement at `c_0`. The
     implication "(AG) implies growth" is therefore a specialisation, not a reduction.
   * *Its `rho = 2` case is the Gaussian-profile version of the circle problem for every imaginary
     quadratic order at once* (Thm 2.3). The write-up notes this. Hybrids (an order-`O` cluster with a
     gap-admissible Gaussian profile) are not excluded. A heuristic check:
     * Suppose the `O`-cluster has an `O`-profile with the same split primes and exponents. Then each
       character value is the argument of an `O`-integer of norm `e^(w(v))`, so the axis (`pi/2`) gaps
       hold by `O`-integrality.
     * The diagonal (`pi/4`) gaps are not automatic, since `O` has no unit `i`. If the relevant arguments
       behave like random angles, these gaps fail with probability about `sum_v e^(-w(v)/2)`, which is
       small when the primes are large.
     * So hybrids are not ruled out by the gaps. They need large `O`-clusters, which is the same open
       problem.
   * *Its quantitative threshold* must be `A(c_0) > 1 + c_0/2` if the dimension-count heuristic is right
     (§2.2(d)), not `A > c_0`.

So the corrected status of `AG` is **generic-fake-violated, not fake-violated**. It is one of a large
family of such inputs. By Cor A.4(IV) of `outside.md`, Haar-null properties of `phi` are allowed in the
barrier classes, but positions are not covered. Any property of positions that fails on a null set of
`delta` is equally "generic-fake-violated", for example "the points have an exact relation of height
`<= e^(W)`". The write-up does not show that `AG` is better placed than these to become a proof. Its
evidence that `AG` is "the cleanest" input is that it interpolates, through Prop 2.8 / Thm R2(2),
between local data and lattice-ness. That is a fair description, and I would keep it as such, not as a
separation theorem.

## 7. Corrected statement

> **Status: partial. No growth improvement and no uniform bound. `M <= (2/3 + eps) log R/loglog R` stands.**
>
> 1. **Collapse (proved).** Thm 2.1 (rank `K = M - 3` iff lattice coset), Thm 2.2 (rigidity, `M >= 5`),
>    Thm 2.3 (lattice cosets on a circle are norm-form problems for an arbitrary imaginary quadratic
>    order), Thm 2.6 (residues mod `N` plus small lifts are the actual problem) and Prop 2.7 (the
>    invariants of `K` are local) hold as stated. Prop 2.7(2) also holds for odd `q | N`. Moreover
>    `covol(K)^2 I^2 = sum_triples (2 Area)^2` exactly, so `K` has independent relations of `l_2` height
>    with geometric mean `<= M^(c_0/4 + 1/2 + o(1))`.
> 2. **The additive gap (proved, corrected heights).** For `M` points on an arc of length
>    `L = C N^(1/4)` with `W <= c_0 M log M`:
>    * every relation-free configuration violates `AG` at `l_1` height `M^(1 + c_0/4 + o(1))`
>      (Thm R2(1); Thm 2.4(1) as written gives `M^(1 + c_0/2 + o(1))`);
>    * `AG(M^A)` forces `rho <= 2 + c_0 M/(4(A-1)) + O(M/log M)` (Thm R2(2));
>    * `AG` at height `e^(W/4 + O(M))` forces a lattice coset.
>
>    The positions with an exact relation are null in `phi`-space for every profile, and in
>    `delta`-space when the profile has rank `M`. All barrier constructions of rounds 4-6 draw `delta`
>    uniformly, so almost every member violates `AG(M^A)` for `A > 1 + c_0/4 + o(1)`. Explicit skew double
>    core instances do (exact LLL, §6.2). Heuristically, typical relation-free configurations satisfy `AG`
>    below `M^(1 + c_0/4)/(2e)`, and fakes with independently imposed relations cannot reach `A > 1 + c_0/2`.
> 3. **Barrier status of `AG` (corrected).** `AG(M^A)`, `A > 1 + c_0/4`, is true for actual circles and
>    false for almost every member of each existing barrier family. It is **not** shown false for every
>    fake of size `O(M log M)`. That statement is Lemma (AG), which is open. Lemma (AG) is strictly
>    stronger than the growth theorem (its restriction to actual clusters *is* growth) and contains the
>    Gaussian-profile circle problem for every imaginary quadratic order. So `AG` is generic-fake-violated
>    only, and no argument is known that uses it short of proving growth.
> 4. **Where `AG`'s content lies (referee Prop R1, proved).** The monomial gaps give
>    `|sum c w| >= sqrt2 e^((W - w(v(c)))/2) - o(1)` for every angle model and every `c` with
>    `v(c) != 0`. So `AG` is automatic (class `A_0`) on combinations of weight `w(v(c)) <= W + log 2 - o(1)`.
>    Every exact relation has `w(v(c)) >= 3W/2 - O(log(n C^3))`.
> 5. **Prop 3.1 (proved, class `A_0`, capped).** As stated. Nearly sharp on data (slack `-0.010`).
> 6. **Thm 3.2 and Cor 3.3 (proved, actual circles only).** As stated, with `N^(1/4) >= 2C/pi` in (c).
>    Thm 3.2(b) is nearly sharp on data (slack `-0.025`). It is a restricted (R1-arch) in low position;
>    fakes in generic position meet its conclusions, so it is not a fake-violated mechanism.
> 7. **§4 (proved).** Lemmas 4.1-4.2, Prop 4.3 and Remark 4.4 as stated, but the "geometric mean"
>    sentence should read `<= M^(c_0/4 + 1/2 + o(1))` (`l_2`).
> 8. **Admissible positions (referee Prop R3, proved).** For rank-`M` profiles the gap-admissible
>    positions have empty interior and positive measure. This answers open point 7.2(a). The finite test of
>    7.2(b) must use high-precision random positions: 53-bit positions create spurious relations.

## 8. Issue list

| # | severity | location | issue | fix |
|---|---|---|---|---|
| 1 | major | §0 item 2, §2.2 last paragraph, §6, "Fake-violated input" | "Provably false for the fakes" is a.e. violation. The barrier needs universal violation, which is Lemma (AG): open, and strictly stronger than growth | Call `AG` generic-fake-violated. State that Lemma (AG) restricted to actual clusters is the growth statement |
| 2 | major | Thm 2.4 asymptotics, Cor 2.5, summary, Prop 2.8 reading, Heur. 2.9, §6 table, §7.1 | `H_1` is the box height. The `AG` (`l_1`) height is `2(M-1)H_1 = M^(1+c_0/2+o(1))`, not `M^(c_0/2+o(1))`. Prop 2.8 at `AG(M^A)` gives `c_0 M/(2(A-1))` | Restate with `l_1` heights |
| 3 | major | Thm 2.4(1)(2), Prop 2.8, Heur. 2.9 | The disc of radius `HML` ignores curvature. The correct region is a `HMC^2 x HML` rectangle. The true threshold is `M^(1+c_0/4+o(1))` (Thm R2; data slopes `0.332`, `0.198` vs `1/(M-3)`). "Violated at `M^(c_0/2+o(1))`" is then proved only for `c_0 >= 4` and heuristically false for `3 <= c_0 < 4`. Heur. 2.9 omits `(M-1)!` and its ceilings `c_0/4`, `c_0` should be `1 + c_0/4`, `1 + c_0/2` | Replace with Thm R2 and the corrected ladder |
| 4 | minor | Cor 2.5 | The positive-measure hypothesis in `R^M` is vacuous for rank-deficient profiles | State it in `phi`-space (exponential polynomial with distinct frequencies). It holds as stated for the barrier families (rank `M`) |
| 5 | minor | §4 "geometric mean `M^(c_0/2+1/2+o(1))`"; round-9 Prop 5.2 | Loose: same curvature loss | Use `covol(K) I = (sum (2 Area)^2)^(1/2)` to get `M^(c_0/4+1/2+o(1))` |
| 6 | minor | Thm 3.2(c) | `|v| <= Q + 1` needs `Q >= 1` | Add `N^(1/4) >= 2C/pi` |
| 7 | minor | §7.2 | The proposed finite test needs high-precision random positions, and the admissible set has no interior | See §1 and Prop R3 |
| 8 | remark | Prop 2.7(2) | Also true for odd `q | N` (proved in §2.3; 13,024 exact cases) | Optional |
| 9 | remark | §6 table, row "lattice lines" | Fine. Thm 3.2 is not fake-violated in its conclusions, as the write-up says | none |

**Bottom line for the coordinator.**
* Keep as proved: §2.1 (Thms 2.1-2.3), §2.3 (Thm 2.6), §2.4 (Prop 2.7), §3 and §4.
* Keep Thm 2.4 and Cor 2.5 in their literal form, with the asymptotic readings corrected as in §2.2.
* Replace the `AG` ladder by the corrected one (Thm R2, thresholds `1 + c_0/4` and `1 + c_0/2`).
* Record `AG` as *generic-fake-violated only*, with Lemma (AG) open and strictly stronger than growth.

Nothing in this write-up moves the growth rate.
