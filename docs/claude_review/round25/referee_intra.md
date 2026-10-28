# Referee report on round25/intra.md (intra-cluster n-point method)

*Referee: adversarial pass on `round25/intra.md`. Scripts: `referee_intra_checks/` (independent code, `rlib.py`;
nothing imported from the author's `glib.py`). All runs under `prlimit --as=3500000000`, one process at a time.*

## 0. Verdict

**The status `reduces_to_route5` and the conclusion "no growth bound, no uniform bound; `M <= (2/3+eps) log R/loglog R`
stands" survive.** Lemma 2.1, Lemma 2.2, Lemma 3.2, Lemma 3.3, Theorem 3.1, Lemma 4.0 and Theorem 4.1 are correct.
I re-derived each with all quantifiers and checked the exact identities with independent code on three actual
clusters. Theorem 5.1 is correct **as stated**, that is, for clusters whose point `Z` is simple and `eta`-balanced
on all `2^{M-1}-1` cuts of `[M]`. The arithmetic of Cor 5.2 and Prop 6.1 is correct.

Three claims do not survive as written. They are all claims about scope.

1. **(major) "Settles round 24's open point 1 for intra-cluster configurations; the cap is `15/4` even with maximal
   prime sharing."** This is false as a statement about intra-cluster configurations. Theorem 5.1 needs `Z` to be
   balanced on **every** cut of `[M]`, so `r >= 2^{M-1}-1` primes. Subtuple balance is much weaker. Linear 5-wise
   uniform cut designs (§5 below) make **every** 5-subtuple exactly balanced with only `r = 2^m - 1` primes. For
   the extended double-error-correcting BCH design (`M = 2^t`, `r = 2M^2 - 1`) with `t >= 13`, route 5's
   **profile-level** LP theorem (`round10/determinant.md` Thm 2.1(d), valid for arbitrary profiles) applies. So at
   these data the intra-cluster `n`-point inequality with star blocks is violated, although every subtuple used is
   exactly balanced at `mu < 4`. Prime sharing therefore **does** beat the one-point cap at the data level. What it
   produces is exact relations (route 5 rigidity), not contradictions. Open point 1 for intra-cluster configurations
   has a negative answer in general: there is no data-level `15/4` cap. The positive answer holds only on the cube
   class. "Maximal sharing" is undefined. The cube is the profile with the most **distinct cuts**. It is not the
   profile with the most sharing.
2. **(major, repairable) Cor 5.2(2) and Cor 5.3 have a template-dependent threshold.** The threshold on `w` contains
   `|C(F)|/idx F`. Cor 5.3 then asserts consistency of "every inference" with one family of data, so its
   quantifiers run in the order "for every `F` there is `w`". The owner's constraint forbids exactly this. The
   repair: only the negative part `C^-` matters, and for the isotropic Taylor remainder bound `C^-(G) <= 1.39 ord_e G`
   uniformly (proved in §5 below). The Dyson-index version of this repair is not proved.
3. **(major) Archimedean scope.** Def 3.0 admits **any** archimedean upper bound `U`. Theorem 5.1 and Cor 5.3 hold
   only for Taylor/Dyson bounds, which see the data only through `idx`/`ord`, the `lambda`'s and constants. "Every
   inference of `L_intra`" and "no intra-cluster inequality excludes it" are therefore overclaims. Route 5's
   direction-sensitive leading form `Phi_c(u)` is one bound they do not cover. For the same reason Cor 4.2
   ("the classes coincide") holds only after "route 5's framework" is read as (TL_M) with an **arbitrary**
   archimedean bound. Lemma 4.0 says route 5's own bound is isotropic, and `L_intra` uses anisotropic Dyson weights.

There are also minor issues (§9). Evidence is labelled as evidence throughout, except for two phrasings: "strict in
practice (47/47)" and "nothing supercritical at `M = 10`".

## 1. What was checked and how

* Every proof in §§2–7 of `intra.md` re-derived line by line (notes: `referee_intra_notes.md`).
* `r1_cluster.py` → `r1_cluster.out` (exact): Lemma 2.1 (1)–(4) on all 5-subsets, Theorem 3.1(1) `F(y) = G(Z)` in
  `Q(i)`, width monotonicity, Lemma 3.3 integrality, sharing gain, on the record 10-point cluster
  (`N' = 1176852625`, recomputed `C = 8.2386`, `lambda(Z) = 3.1127`, both matching the author), its 9-point
  subcluster, and an independent 8-point window of `N = 5·13·17·29·37·41·53` (not the author's third cluster).
* `r2_regB.py`, `r2b_structured.py` → `.out` (exact rank over `Q`): the zero-sum-slice mean-width bound (5.1) on 95
  random supports (`M = 5,6,7`) and on structured near-extremal supports (`R_5`, `2R_5`, `R_6`, `R_7`, slices of
  `{-1,0,1}^5`, `{-1,0,1}^6`). Also `reg < 2 max|k|_1`, and Theorem 4.1(1)–(2) at `M = 9, 13`, `T = 1,2,3`.
* `r3_design.py` → `r3_design.out` (exact counting): 5-wise uniform linear cut designs, exact balance of all
  5-subtuples, number of cuts of `[M]` used, and the route-5 regime arithmetic for the extended BCH designs.

## 2. Lemma 2.1 / 2.2 (projection, gcd reduction): correct

* (1) Since `sum s_j = 0`, `g_S^{-sum s}` and `z_1^{-sum s}` both drop out. Correct for every representative, so no
  phase is fixed.
* (2) `v_pi(g_S) = min_S a(p)`, because the `pi`-part of each `z_j` is `pi^{a_j}` and the `pibar`-part is
  `pibar^{e-a_j}`. So `v_pi(zeta_j) = a_j - min`, and (2.1) is the pairing `<iota_S s, v> = <s, v|_S>`. Correct.
* (3) `N(g_S) = prod p^{min a + e - max a}`, so `h_B = log(N/N(g_S)) = sum (max-min) log p`. The lower bound uses
  `g_S | z_a - z_b`, which gives `N(g_S) <= L^2` and `h_B >= W - 2 log L`. The bound `>= W/2` needs `L <= sqrt R`,
  that is `C <= 1`. For `C > 1` it becomes `h_B >= W/2 - 2 log C`. This matters only for the constant `kappa <= 2`
  quoted in Prop 6.1(d), where `kappa <= 2 + O(log C / W)` is the honest form (minor).
* (4) `dist(y_S,e) = max_{a,b in S}|z_a-z_b|/R`. Correct.
* Independent exact check (`r1_cluster.out`): (1)–(4) hold on all 252 / 126 / 56 five-subsets. `h_B/W` lies in
  `[0.7737, 1]` on the record cluster, matching the author, and in `[0.8286, 1]` on my 8-point window.
* Lemma 2.2. Monomial by monomial `m(z) = prod g^{a_i} gbar^{b_i} m(zeta)` (Lemma K1). The toric bound
  `|F(zeta)|_w <= max_m |m(zeta)|_w <= 1` implies integrality. By the product formula the toric bound gives a lower
  bound for `log|F(zeta)|` that is at least as large as integrality's `0`. After the common factor
  `z_1^{sum(a-b)}N^{sum b} prod g^{-a} gbar^{-b}` is divided out, it is (I1) for `G` at `Z`. The common factor is
  tracked, so no gcd or denominator cost is discarded. Correct.

## 3. Theorem 3.1 (reduction to one-point Liouville on T_M): correct, and close to a tautology

* Lemma 3.2(a): `G(Z) = F(y)`, monomial by monomial. (b): in log coordinates `G = F∘dPhi`, so Taylor parts of degree
  `< idx F` vanish. (c): if `dPhi` is invertible, the degree-`m` part vanishes iff all `F_J` with `|J| = m` vanish.
  Correct.
* Lemma 3.3 is the product formula over `Q(i)`, with the inert primes and `1+i` units because `N` is odd and has only
  split primes. Exact check: `Norm(G(Z)) prod p^{width}` is an integer `>= 1` in every random case (`r1_cluster.out`).
* Theorem 3.1(2): both widths are max–min of the same numbers `<Phi^*s, v_pi(Z)>`, and `supp G ⊆ Phi^*(supp F)`.
  Correct.
* Reading. The reduction is exact, but it says little more than this: `F(y)` and `G(Z)` are the **same** element of
  `Q(i)`, and every valuation of the `y_S` is a linear image of `v_pi(Z)`. The theorem is right. Its content is the
  observation, not the inequality.
* "Strict in 47/47 exact cases (`idx < ord`)": the author's mode 1 **builds** `F` by imposing `T_M`-order on
  overlapping blocks, where `dPhi` is not injective. Strictness is a property of that construction, not something
  seen "in practice" (minor wording).
* Sharing gain (Remark 3.4). The numbers depend on the sample of `F`: my independent random sample gives
  `[0, 134]` against the author's `[3.5, 64.3]`. The gain is nonnegative by definition, and comparing it with `W`
  has no meaning because both sides scale with block degrees. Evidence only. The structural statement (the gain is
  the width of one polytope at one vector) is correct.

## 4. Theorem 4.1 / Cor 4.2 (converse, star partition): correct, with one scope correction

* Lemma 4.0. `G_k(z) = N^{(D-|k|_1)/2} prod z^{k^+} zbar^{k^-}`. With `zbar = N/z` and `|k^+|_1 = |k^-|_1 = |k|_1/2`
  this is `N^{D/2} Z^k`. So (I1) is Lemma 3.3 for `G_c` at `Z`, and `V(A,m)` is `ord_e G_c >= m`. Correct.
* Theorem 4.1(1)–(2). `Phi^*` maps the basis `(e_x - e_o)_{x in B_j}` to itself. Complex conjugation acts by `-1` on
  both character lattices, so `Phi` is a `Q`-isomorphism. Also `|s_j|_1 <= 8T`, which means `s_j in 4T R_5`. Exact
  check: 500 random box exponents for each `T = 1,2,3` at `M = 9, 13` pull back uniquely and land in `4T R_5`
  (`r2_regB.out`). Note that route 5's refereed slab box has `T = n` (`round10/determinant.md` Thm 2.1), so the block
  degrees in Cor 7.1 are `4n`. Harmless.
* Theorem 4.1(4), wording (minor). After padding, `Phi^*` is surjective and not injective, so `Phi` is a closed
  immersion onto a subtorus, not a "surjection". `(Phi^*)^{-1}` in Cor 7.1 then denotes a chosen section. For a
  clean statement one can instead take `M ≡ 1 (4)` (any `4q+1` points of a design remain a design).
* **Cor 4.2 (scope).** The chain is `L_intra ⊆ {(TL_M) at Z, any G} ⊆ L_intra`. It is correct only if the archimedean
  side is unrestricted on both ends, and the middle term is (TL_M) with an arbitrary upper bound for `|G(Z)|`. Route
  5's own archimedean bound is the **isotropic** Taylor remainder (Lemma 4.0; `determinant.md` (4.1)). `L_intra` uses
  the anisotropic Dyson bound `min_J sum j_i lambda(y_{S_i})`, which can be strictly stronger when some subtuples are
  much tighter than the cluster. So the honest statement is: route 5's inferences are `L_intra` inferences (star
  blocks); `L_intra`'s non-archimedean side is (TL_M) for `G = Phi^*F`; and the archimedean refinements are
  refinements of the Taylor expansion of `G` at `Z`, which route 5 has not used. "The two methods coincide" holds
  for the polynomial class and the non-archimedean side, not for the archimedean bookkeeping. This does not affect
  the verdict.

## 5. Theorem 5.1 / Cor 5.2–5.3 (cube barrier): correct on its hypotheses; three scope errors

**5.1 What is correct.** I re-derived each step.
* (5.1): for `S ⊂ M_M`, `w_P = w_{P^c}` and `w_∅ = 0`, so `sum_P w_P = 2 sum_cuts`, and Theorem B gives
  `reg <= 2^{2-M} sum_cuts width_S(1_A)`.
* `reg` on the slice equals `reg` on `T_M`, since the quotient `G_m^M -> T_M` is smooth.
* Hence `RHS_G >= 2^{M-3}(1-eta) w ord_e G` on simple, lower-balanced profiles.
* `LHS(*) <= idx F · max_i lambda_i`.
* Cor 5.2(1): each of the 15 cuts of `S` is induced by exactly `2^{M-5}` cuts of `[M]`.
* The identity `2^{M-3}(1-eta) - (2^{M-1}-1)(1+eta)/4 = (1-(2^M-1)eta)/4` holds.
* The window `mu(y_S) >= (4-2^{3-M})(1-eta) - log C/(2^{M-5}w)` is correct.

Exact checks (`r2_regB.out`, `r2b_structured.out`) agree:
* (5.1) holds on 95 random supports and on all structured ones (`R_5`: `reg 3 <= 3.75`; `2R_5`: `6 <= 7.5`, which is
  `t_max(2) = 6`; slice of `{-1,0,1}^6`: `6 <= 8.25`).

**5.2 Major: Theorem 5.1 does not settle open point 1 for intra-cluster configurations.** Round 24's open point 1
asks for a cap on (*) when the **points** share primes. Here the points are the subtuples `y_S`, and the round-24
cap is about points that are each `eta`-balanced. Theorem 5.1 assumes much more. It asks that **`Z`** be balanced on
**all** `2^{M-1}-1` cuts of `[M]`, so `r >= 2^{M-1}-1`. Subtuple balance only asks that the cut family be 5-wise
uniform.

*Proposition R1 (referee; exact part (a)–(b), proved part (c) from refereed route-5 theorems).* Let `c_1..c_M in
F_2^m` have every 5 of them linearly independent. Take one prime `p_u` for each `u in F_2^m \ {0}` (so
`r = 2^m-1`), with simple profile `a_x(p_u) = <u, c_x> in {0,1}` and `log p_u in [(1-eta)w_0, (1+eta)w_0]`.
* (a) For every 5-subset `S`, the map `u -> (<u,c_x>)_{x in S}` is onto `F_2^5` with fibres of size `2^{m-5}`.
  Hence each of the 15 cuts of `S` carries exactly `2^{m-4}` primes, and the gcd class carries `2^{m-4}-1`. Every
  `y_S` is simple and `eta`-balanced at scale `w_S ≈ 2^{m-4} w_0`.
  If some `u` has `<u,c_x> = 1` for all `x` (for ext-BCH, `u = e_1`), its prime divides every `z_x`. That violates
  `gcd(z) = 1`, so drop it (`r = 2^m - 2`). It lies in the gcd class of every `S`, so (a) is unchanged except that
  the gcd class carries `2^{m-4}-2`.
* (b) `Z` carries at most `r <= 2^m - 1` of the `2^{M-1}-1` cuts of `[M]`. For ext-BCH it carries `2^{m-1}`
  including the trivial one, since `u` and `u + e_1` give complementary patterns. So for `M > m + 1` Theorem 5.1
  does not apply.
  Exact (`r3_design.out`):

  | design | `M` | `m` | cuts of `[M]` carried |
  | --- | --- | --- | --- |
  | greedy | 7 | 6 | 32 / 63 |
  | greedy | 9 | 7 | 64 / 255 |
  | ext-BCH, `t = 3` | 8 | 7 | 64 / 127 |
  | ext-BCH, `t = 4` | 16 | 9 | 256 / 32767 |

  All of these have every 5-subtuple exactly balanced.
* (c) Take the extended double-error-correcting BCH design: `M <= 2^t`, `m = 2t+1`, `r = 2^{2t+1}-1`. That every 5
  columns are independent is the minimum distance 6 of that code. This is standard coding theory, cited from memory
  (MacWilliams–Sloane, Ch. 9), and checked exactly for `t = 3, 4`. Take `t >= 13` and `M ≡ 1 (4)` with
  `M - 1 >= 300 ln(4r)`, for example `t = 13`, `M = 8189`, `8188 >= 6030`. Then `round10/determinant.md` Theorem
  2.1(d), which is stated for **arbitrary** columns `beta_p` and weights `omega_p`, gives the slab box with
  `m' >= 4w`. Its proof of Theorem 4.1 shows that, once `W >= 16 ln n + 8 log^+ C + 14`, the Liouville comparison
  `(m'/2 - w)W <= 2 ln|c|_1 + ...` fails for every `c in V_H \ {0}` unless `F_c(z) = 0`. By `intra.md` Theorem 4.1
  this is the `n`-point comparison for `F = (Phi^*)^{-1} G_c` at the `q = (M-1)/4` star subtuples.
* (d) Let `eta = 2^{-(2t+4)}`, which the prime number theorem in progressions allows once `w_0 >= (4t+5) log 2 +
  O(1)`. Let the cluster have diameter `C sqrt R` and every star subtuple diameter `>= C^{-1}·C sqrt R`, that is
  `kappa <= C`. Then
  `mu(y_S) <= 4(1-2^{-m})(1+eta)/(1-eta) + log(kappa/C)/w_S < 4(1-eta)`,
  because `-2^{1-2t} + 2^{-2t-1} + 2^{-2t-2} < 0`.

So at these data **every subtuple used is exactly balanced and strictly inside the one-point window `mu < 4(1-eta)`**.
For pairwise disjoint primes, round 24's Prop 4.3(iii) guarantees that (*) holds there. Yet the intra-cluster
`n`-point comparison is violated unless `F(y) = 0`. Shared primes therefore **do** beat the one-point cap at the
data level, through route 5's LP supercriticality. The output is an exact relation, not a contradiction, so nothing
changes for growth. But the correct answer to open point 1 for intra-cluster configurations is **no data-level cap
in general**. Round 24's assessment ("suggests the same cap") is wrong for these configurations. `intra.md` already
contains the ingredients (Cor 7.1). It never checks that route-5-regime subtuples can be balanced, and so its §0
item 4 and summary state the opposite. "Even with maximal prime sharing" has no definition. In both the cube and
the design, a prime separating a pair is shared by about `15/16` of the 5-subsets. What separates the two cases is
the number of distinct cuts, `log r` against `M` (route 5's `n - 1 >= 300 ln 4r`).

*Status of Prop R1:*
* (a) and (b) are exact.
* (c) is proved from two inputs: the refereed route-5 theorems, and the cited BCH fact. The BCH fact is used only
  for `t >= 13`. Any 5-wise independent family with `log r = o(M)` would do, and random families also work by a
  union bound (not written out).
* (d) is arithmetic.
* Like Theorem 5.1, R1 is a statement about the method on data. It does not assert that such clusters exist on short
  arcs.

**5.3 Major (repairable): template-dependent threshold in Cor 5.2(2) / Cor 5.3.** The proof of Theorem 5.1 needs
`LHS - C(F) <= RHS_F`. When `C(F) >= 0`, `LHS <= RHS_F` already suffices, so only `C^-(F) := max(0, -C(F))`
matters. The author writes `|C(F)|`. With `C(F) << log H(F) + sum d_i log(2+d_i)` (round 24), the threshold
`w >= 4(log(kappa/C) + |C(F)|/idx F)/(1-(2^M-1)eta)` grows with the family. Cor 5.3 then claims that *every*
inference is consistent with one family of data. That is "for every `F` there is `w`" presented as "there is `w` for
every `F`", which is the owner's forbidden step. The repair for the isotropic Taylor bound is the following lemma.

*Lemma R2 (referee; proved).* Let `0 != G = sum_{k in S} c_k x^k`, `c_k in Z[i]`, `S ⊂ M_M`, `m = ord_e G >= 1`,
`K_1 = max_S |k|_1`. Let `Z = exp(i theta)` with `delta = max|theta_a - theta_b| <= 1`. Then
`log|G(Z)| <= -m log(1/delta) + C_T(G)`, where `C_T(G) = log|c|_1 + m log(K_1/2) - log m! + K_1 delta/2`, and
`C_T(G) >= -(log 4 - 1)m - 1 - log m >= -1.39 m`.
*Proof.*
1. Since `ord_e G >= m`, `sum c_k <k,i theta>^j = 0` for `j < m`.
2. `|e^x - sum_{j<m} x^j/j!| <= |x|^m e^{|x|}/m!`.
3. `|<k,theta>| <= |k|_1 delta/2`: centre `theta` at the midpoint, using `sum k = 0`.
4. This gives the upper bound.
5. For the lower bound: `|c|_1 >= 1`.
6. `width_S(1_A) <= K_1`, because `|k(A)| <= |k|_1/2`. So (5.1) gives `m <= reg S <= 2^{2-M}(2^{M-1}-1)K_1 < 2K_1`,
   that is `K_1/2 > m/4`.
7. `log m! <= m log m - m + 1 + log m`. `□`

*Exact check:* `reg < 2K_1` holds in all cases (`r2_regB.out`, `r2b_structured.out`; equality `reg = K_1` occurs on
the `{-1,0,1}^6` slice).

Since `log(1/delta) = lambda(Z) + O(delta^2)`, Theorem 5.1(2) and the TL_M half of Cor 5.3 hold **uniformly in
`G`** under `w >= 4(log(kappa/C) + 1.39 + o(1))/(1-(2^M-1)eta)`. That threshold depends on no template. For the
Dyson form (Theorem 5.1(1)) the analogous blockwise bound `C^-(F) <= c·idx F` is plausible by the same Stirling
argument but is **not proved**, here or in `intra.md`. Until it is, Cor 5.3 for `L_intra` holds only for families
with `C^-(F) = o(idx(F) w)`.

**5.4 Major: archimedean scope of Cor 5.3.** Def 3.0 allows "any upper bound `U`". Theorem 5.1 bounds only the Dyson
or Taylor forms, which see the data only through `idx`/`ord`, the `lambda`'s and a constant. An auxiliary
polynomial chosen **after the positions**, as Def 3.0 explicitly allows, can have a small leading form at the actual
direction `theta/|theta|`. That is route 5's `Phi_c(u)` mechanism. Its archimedean bound is then not of Taylor type,
and Theorem 5.1 says nothing about it. "Every inference of `L_intra` … any zero estimate" and "no intra-cluster
inequality excludes it" must be restricted to Taylor/Dyson-type archimedean bounds. With that restriction, Cor 5.3's
negative conclusion holds (with 5.3's repair): it gives no uniform bound and nothing better than
`M ≈ log_2 log R` on cube data. Direction-sensitive bounds are where route 5's open MIN^pair / (R_H) problem lives,
so the restriction does not open a new route. It only stops the barrier from being stated more strongly than it is
proved.

*Minor.* Cor 5.3 writes `W >= 2^{M-1}(1-eta)w`. It should be `(2^{M-1}-1)(1-eta)w`.

## 6. Prop 6.1 (lacunary degrees): algebra correct; the dichotomy is not a theorem

* (a) Pick `s, s'` realising `width_{P_1}(v_1)`. Each other block changes the pairing by at most
  `d_i width_{R_5}(v_i)`. Correct. The author checked it on 300 random `F`; I re-derived it.
* (b) `lambda(y) <= h_B(y)/2` holds because `|zeta_a - zeta_b| >= 1` and `|y_a/y_b - 1| = |zeta_a-zeta_b|/R''`.
  Also `R'_i <= floor(10 d_i/3)`: my exact ranks give `reg(R_5) = 3` and `reg(2R_5) = 6`. Correct.
* (c) and (d): `(1 + 5/3) = 8/3`, `(8/3)·15 = 40`, and `RHS_{P_1} >= 4(1-eta)w_1 R'_1` by (5.1) at `M = 5`. Correct.
  The constant `kappa <= 2` needs `C <= 1` (§2).
* `R'_1 >= rho d_1` is an **assumption** about the construction. It is labelled as such in (d), but the summary
  states the conclusion without it (minor).
* The "dichotomy" has two halves.
  - The lacunary half is proved, given the assumption above.
  - The other half, "with comparable degrees no product zero lemma applies, and the zero step is route 5's ZE", is
    not a theorem. Dyson-type lemmas in the Esnault–Viehweg / Nakamaye form hold for all degree vectors (they are
    weak, not inapplicable, when degrees are comparable). Those lemmas are on `(P^1)^n` and the product-theorem
    analogues are cited from memory. Also, "`G(Z) != 0` is a base-locus statement on `T_M`" is a tautology.
  The summary's item 5 should say "heuristic dichotomy; the lacunary half is proved". This does not affect the
  verdict.

## 7. Cor 7.1 (growth regime, route 5): correct

* The hypotheses quoted (`n >= 64`, `n-1 >= 300 ln 4r`, `W >= 16 ln n + 8 log^+ C + 14`) are exactly those of
  `round10/determinant.md` Thm 4.1. Its proof (`(m/2 - w)W <= 2 ln|c|_1 + 2m ln(DC/2) - 2 ln m!`, with
  `m/2 - w >= m/4`) uses Lemma 3.3's bound and the isotropic Taylor remainder, so transporting it through Theorem
  4.1(2) is legitimate.
* Route 5's box has `T = n`, so the block degrees are `4n` (§4). For `M` not `≡ 1 (4)` the map `(Phi^*)^{-1}` is a
  section (minor).
* The template-dependence point about `pr_S(B)` is correct and is an honest application of the owner's constraint.
* Cor 7.1 together with Prop R1 is precisely the reason the open-point-1 claim of §0 item 4 fails (§5.2 above).

## 8. Independent computations (all exact unless marked)

| script → output | result |
| --- | --- |
| `r1_cluster.py` → `.out` | Record 10-pt cluster: `C = 8.2386`, `lambda(Z) = 3.1127` (author's values reproduced); 9-pt subcluster; independent 8-pt window `N = 2576450045`, `C = 15.01`. Lemma 2.1(1)–(4): all True. `F(y) = G(Z)` exactly in `Q(i)` (60 random `F` per cluster, 2–4 blocks): True. `width_{supp G} <= width_{Q_F}`: True. Lemma 3.3 integral `>= 1`: True. Sharing gain `[0, 134]` (evidence; sample-dependent) |
| `r2_regB.py` → `.out` | (5.1) on 95 random supports, `M = 5,6,7`: 0 violations; `reg < 2K_1`: 0 violations. Star partition `M = 9, 13`, `T = 1,2,3`: unique preimage, `s_j in 4T R_5`: True |
| `r2b_structured.py` → `.out` | `R_5`: `reg = 3`; `2R_5`: `reg = 6`; `{-1,0,1}^5` slice: 5; `R_6`: 3; `{-1,0,1}^6` slice: 6; `R_7`: 3. All satisfy (5.1) and `reg < 2K_1` |
| `r3_design.py` → `.out` | 5-wise linear designs: any 5 columns independent, every 5-subtuple exactly balanced, cuts carried 32/63, 64/255, 64/127, 256/32767. Ext-BCH: route-5 condition `n-1 >= 300 ln 4r` first holds at `t = 13` (`M = 8192`, `8191` vs `6030`) |

On the author's evidence (§8 of `intra.md`), check3b's `reg` is a mod-`p` rank. That gives an **upper** bound on
`reg`, and dropping `-m log K - log|c|_1` also only raises `E`, so "`max E < 0`" is conservative. It covers only the
scanned family (`k' in {-1,0,1}^9`, boxes `b_p <= 2`, 190 boxes). "Nothing is supercritical at `M = 10`" must read
"nothing in the scanned family" (minor; the table already says evidence).

## 9. Issues list

| # | severity | issue | repair |
| --- | --- | --- | --- |
| M1 | major | §0 item 4 and the summary say Theorem 5.1 "settles round 24's open point 1 for intra-cluster configurations" with cap `15/4` "even with maximal sharing". Theorem 5.1 needs `Z` balanced on all `2^{M-1}-1` cuts. 5-wise uniform designs make every subtuple exactly balanced with `r ≈ 2M^2`. In route 5's regime (ext-BCH, `t >= 13`) the intra-cluster comparison is violated at data where every subtuple used is balanced at `mu < 4(1-eta)` (Prop R1). | Restrict the claim to cube-balanced `Z`. State that open point 1 has a **negative** answer for intra-cluster configurations in general (sharing beats the one-point cap at the LP level), with output exact relations whose use is route 5's ZE. Drop "maximal sharing". |
| M2 | major (repairable) | Cor 5.2(2) and Cor 5.3: the threshold on `w` contains `|C(F)|/idx F` (template-dependent), while Cor 5.3 asserts one family of data defeats every inference. | Use `C^-` only. Lemma R2 gives `C_T^- <= 1.39 ord` for the isotropic Taylor bound, which is uniform. The Dyson-form analogue is unproved, so either prove it or restrict to families with `C^-(F) = o(idx F · w)`. |
| M3 | major | Def 3.0 admits any archimedean `U`, but Theorem 5.1 and Cor 5.3 hold only for Taylor/Dyson-type bounds. "Every inference of `L_intra`" and "no intra-cluster inequality excludes it" are therefore overclaims. Cor 4.2 equates `L_intra` (anisotropic) with route 5 (isotropic, Lemma 4.0). | Restrict Cor 5.3 to Taylor/Dyson-type archimedean bounds. Read Cor 4.2 as a statement about the polynomial class and the non-archimedean side. |
| m1 | minor | Prop 6.1 "dichotomy": the zero-lemma half is cited from memory and is not a theorem. Esnault–Viehweg/Nakamaye Dyson lemmas are not inapplicable at comparable degrees. "Base-locus statement on `T_M`" is tautological. `R'_1 >= rho d_1` is an assumption. | Label it heuristic. |
| m2 | minor | `kappa <= 2` (Prop 6.1(d)) and `h_B >= W/2` need `C <= 1`. | `kappa <= 2 + O(log C/W)`. |
| m3 | minor | "strict in practice (47/47)": strictness is built into mode 1 (overlapping blocks). The sharing-gain range `3.5–64` is sample-dependent (mine: `0–134`), and comparing it with `W` has no meaning. | Rephrase as construction-dependent evidence. |
| m4 | minor | check3b: "nothing supercritical at `M = 10`" holds only within the scanned slab family. | Add "in the scanned family". |
| m5 | minor | Theorem 4.1(4): with padding, `Phi` is a closed immersion, not a surjection, and `(Phi^*)^{-1}` is a section. Cor 5.3 has `2^{M-1}` where it should have `2^{M-1}-1`. Route 5's box has `T = n` (block degree `4n`). | Wording. |

No issue touches Lemma 2.1, Lemma 2.2, Lemma 3.2, Lemma 3.3, Theorem 3.1, Lemma 4.0, Theorem 4.1(1)–(3), Theorem 5.1
as stated, or Cor 7.1. The bookkeeping of shared-prime divisibility is exact: no Gaussian gcd or denominator cost is
discarded, no common phase is set to zero, and every evaluation is of a degree-0 character. No fixed-`S` Subspace
exceptional set is used.

## 10. Corrected statement

**Status: reduces_to_route5 (confirmed). There is no growth bound and no uniform bound;
`M <= (2/3+eps) log R/loglog R` stands.**

1. *(proved)* The `binom(M,5)` gcd-reduced subtuple points are the projections `pr_S(Z)` of the one cluster point
   `Z in T_M(Q)`. Their heights, gcds, valuations and shared primes are functions of `v_pi(Z)` (Lemma 2.1). The gcd
   reduction is invisible to forced divisibility (Lemma 2.2, multi-point K1).
2. *(proved)* Every intra-cluster `n`-point Liouville comparison is implied by toric Liouville (TL_M) at `Z` for
   `G = Phi^*F`, with the same value, `RHS_G <= RHS_F` and `ord_e G >= idx F` (Theorem 3.1). Conversely, a star
   partition gives `T_M ≅ T_5^{(M-1)/4}`, and route 5's slab-box forms and (I1) are intra-cluster forms and (TL_M)
   (Theorem 4.1). The two coincide **as classes of auxiliary polynomials with their non-archimedean bounds**.
   `L_intra`'s anisotropic Dyson weights are an archimedean refinement that route 5 does not use.
3. *(proved)* **Cube barrier.** Suppose `Z` is simple and `eta`-balanced on **all** `2^{M-1}-1` cuts, with
   `eta < 1/(2^M-1)`. Then no comparison with an isotropic Taylor upper bound excludes it, for any `G`. This holds
   uniformly in `G` once `w >= 4(log(kappa/C) + 1.39 + o(1))/(1-(2^M-1)eta)` (Lemma R2). The same holds for Dyson-
   form intra-cluster comparisons, for families with `C^-(F) = o(idx F · w)`. Any zero estimate leaves this
   unchanged. The subtuples sit at `mu >= 4 - 2^{3-M} - O(eta)`. So within Taylor/Dyson-type archimedean bounds the
   class cannot give the uniform bound, nor anything below `M ≈ log_2 log R`. Direction-sensitive archimedean bounds
   (route 5's `Phi_c(u)`) are not covered.
4. *(proved from refereed route-5 theorems; R1 above)* **Round 24's open point 1, intra-cluster case: no data-level
   cap in general.** Consider 5-wise uniform cut designs, for example extended BCH with `t >= 13` and `M ≡ 1 (4)`,
   `M >= 8189`. Every 5-subtuple is exactly balanced and lies strictly inside the one-point window `mu < 4(1-eta)`.
   Yet the intra-cluster `n`-point comparison on the star subtuples is violated (route 5 LP supercriticality)
   unless `F(y) = 0`. So shared primes **do** beat the one-point cap `15/4`, but only by forcing exact relations.
   Using those relations needs route 5's zero estimate (ZE-col*; MIN^pair / (R_H)), which is open.
   - The cap holds on the cube class (item 3).
   - Round 24's "suggests the same cap" is false for intra-cluster configurations.
   - The disjoint-primes hypothesis in the expert note's anisotropy statement is therefore essential.
5. *(proved for lacunary degrees under `R'_1 >= rho d_1`; heuristic otherwise)* With lacunary degrees the `n`-point
   comparison is the leading-block one-point comparison up to `O(kappa omega/rho)`. The claim that comparable
   degrees leave only route 5's ZE as a zero step is heuristic.
6. *(proved)* In route 5's regime the method is LP-supercritical and outputs exact relations among star subtuples
   (Cor 7.1). Growth needs route 5's zero estimate, which is open. `T_5` counting cannot supply it, because the
   base locus is template-dependent.
7. *(evidence)* At `M = 10`, within the scanned slab family, nothing beyond the monomial and pair relations is
   forced. The monomial relations are fake-blind.
