# Referee report on round26/ze.md

Referee: round-26 adversarial referee. Date: 2026-10-28. Target:
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round26/ze.md` (status claimed:
obstruction_found).

## 0. Verdict

**The result survives in part; four major corrections are needed.** No growth or uniformity claim is made, and
none is needed; there is nothing to refute there.

**What I could not break:**
* Lemma 3.1 (Liouville at `Z^j`, all gcd costs kept). I re-derived it.
* The spectral LP cap (Prop 3.3).
* The degree-level vacuity of Philippon's theorem in both pictures (Props 3.5, 4.4). All margins were reproduced
  independently.
* The structure of the dual obstruction subgroups (Lemma 4.2, Prop 4.3).
* The nondegeneracy of `Z` (Prop 5.1).

**What fails:**
1. **M1.** "Route 5 certifies vanishing only at the powers `Z^j`" is false. The same Liouville argument certifies
   `P_c(Z^M) = 0` for integer endomorphisms `M` of low cost. For example, every profile automorphism costs exactly
   `Lambda(c)`, and the Sylvester profile has `>= |AGL(k,2)| ~ e^{107}` of them at `n = 4096`.
2. **M2.** So Thm 3.6's application ("for every point, a form carrying **all** route-5 data, every
   `mu <= m_max+1`") is not established. With the `AGL` orbit included, at ze.md's own orders, the exact count fails
   at `mu = m_max` for `n = 4096, 8192`. It holds only to `mu = 2559` and `5186`.
3. **M3.** Even for powers-only data the barrier holds only for `lambda >= lambda*(n,mu)` (`0.0168` at
   `n = 2966`, `m_max`). The summary states it without this template-dependent qualifier.
4. **M4.** The dual picture concerns `Bs(V(A,m))`, but route 5 only has `Z in Bs(V_H)` with
   `dim V_H >= |A|/2`. So "ZE^all ⇒ ZE" and "ZE-tor^glob ⇒ ZE" are unproved. The right target is ZE^all_H.

The negative verdict on the two natural framings stands. The universal "no single-form estimate with any constants"
barrier is downgraded to the restricted statement of §5.

## 1. Section-by-section check

Labels: **[confirmed]** (I re-derived it or recomputed it independently), **[gap]**, **[error]**,
**[overclaim]**. Independent code: `ref_ze/recompute.py`, output in `ref_ze/recompute_out.txt` (paths are under
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round26/`).

### 1.1 Lemma 3.1 (forced vanishing at `Z^j`): [confirmed]

I re-derived Step 1 from scratch.
* Write `G_k(y) = N^{j(D-|k|_1)/2} prod y^{k+} conj(y)^{k-}`. At a split `p` the `pi`- and `conj pi`-valuations
  are `j(e_p D/2 +- <k,a_p>)`. Inert primes and `1+i` divide every `G_k` to the full power `N`-share; I checked
  this by hand for `q^f || z_x`, which gives `q^{jfD}` in every `G_k`.
* The minimum over `supp c'` gives `Norm F >= N^{jD} prod_p p^{-j h_{supp c'}(a_p)} = N^{jD} e^{-jW Lambda(c')/2}`.
  Here `F` is a nonzero Gaussian integer divisible by that ideal, so no gcd cost is discarded. The minimum is
  taken over the support, not over `A`. That is a sharpening, and it is legitimate.
* The upper bound uses the Taylor remainder after `mu'` vanishing moments, with
  `|<k, j delta>| <= |k|_1 |j| Theta/2`.
* Comparing the two bounds gives exactly (3.1) at `s = 0`.
* Step 2 is also right. `theta^alpha P_c = P_{c k'^alpha}`, order drops by `|alpha|`, and `|c^(alpha)|_1` grows
  by at most `T_A^{|alpha|}`.
* Degree-0 characters cancel the common phase, so the owner's constraints are respected here.

### 1.2 The claim that route 5 certifies vanishing **only** at the powers `Z^j`: [gap, major]

§0 item 1, §3 (opening), §3.4 ("all vanishing data route 5 can certify") and §6.3 ("Other points ... So
`Sigma ⊂ <Z>` is forced") all assert that, for one form, the certified points are exactly `Z^j`. Lemma 3.1 proves
that these points are certified. It does not prove that nothing else is, and the claim is false. The proof of
Lemma 3.1 works verbatim for any tuple `y_x = prod_y z_y^{M_xy}` (with `conj z_y` for negative entries, and padding
by powers of `N` to equalise norms). Here `M` is an integer `n x n` matrix. The certified point is the image `Z^M`
of `Z` under the endomorphism `M` of `G_m^d`, and its cost is
`Lambda_M(c) = sum_p omega_p h_{supp c}(M beta_p)`. Two explicit families:

* **(a) Profile automorphisms.** Let `sigma` be a permutation of the `n` points that maps every column `beta_p` to
  `+- beta_{p'}` and preserves the weights. Then `Lambda_sigma(c) = Lambda(c)` exactly, because `h_S` is even. So
  route 5 certifies `P_c` vanishing at `Z^{j sigma}` with the **same** orders `tau_j` as at `Z^j`. The derivative
  step can be run in the coordinates that drop `sigma(n)`, so `T_A = 1` is kept.

  For the Sylvester profile that ze.md itself uses in Thm 3.6 (`n = 2^k`, columns `beta_u = ((-1)^{u.x})_x`,
  `u != 0`, uniform weights), the affine group `AGL(k,2)`, acting by `x -> Ax + b`, maps `beta_u` to
  `+- beta_{A^T u}`. So `|Aut| >= |AGL(k,2)|`, with `ln|AGL(12,2)| = 106.89` and `ln|AGL(13,2)| = 124.91`
  (`recompute_out.txt`).
* **(b) Light endomorphisms.** Suppose that some prime `p` exists for which the other columns together with `1` do
  not span. Take `v` in the integer left kernel of `[a_q (q != p), 1]` and set `M = 1 m^T + u v^T`. Then `Z^M`
  equals `(Z^v)^{u'}`, a point on the one-dimensional subtorus `t -> t^{u'}`, and it costs
  `omega_p |<v,beta_p>| h_S(u)`, i.e. one light prime. Whether such `p` exist in regime profiles is
  profile-dependent. I record (b) only to show that the set of certified points is not `<Z>` in general.

Consequences:
* Prop 3.5 is correct **as stated**, because (3.4) is written with `Sigma(d) ⊂ {Z^j}`.
* For family (a), the real (sumset) theorem is still vacuous. Every `s = Z^{M}` with `s != e` needs `s^d` to be
  certified, and `Lambda_{dM} = d Lambda_M = d Lambda(c) > ord(c)` for `M = sigma`. For family (b), a proof of
  singleton-ness would need `Lambda_M(c) >= ord(c)/d` for every certifiable `M` with `Z^M != e`. ze.md does not
  prove this.
* Under ze.md's own convention for Thm 3.6 (orders bounded independently of `c`), the **application** ("carries all vanishing data route 5 can certify") fails at `mu = m_max` for `n = 4096` and `8192` once the `AGL(k,2)` orbit is added. The barrier is therefore **not established** there; it is not disproved. See §1.5 and §2.

### 1.3 Prop 3.3, Cor 3.2, Cor 3.4: [confirmed]

**Prop 3.3.**
* `Omega(S) = E_eps h_S(eps)` is (I2) (`threshold/construct.md` Cor B1, which I re-read).
* Step 3 is the identity `G G^+ P eps = sum_p omega_p bt_p <G^+ bt_p, eps>`. It uses `G^+` symmetric with range
  `1^perp`.
* Steps 4–6 are a seminorm bound, Jensen, and `|G^+|_op = 1/lambda`, `|beta_p|_2 <= sqrt n`. All correct.
* Sylvester `lambda = n/(n-1)` is exact: `sum_{u != 0} beta_u beta_u^T = nI - 11^T`.
* Note: `lambda = n/(n-1)` is the **largest** value `lambda` can take, because `tr G <= sum omega_p |beta_p|^2 <= n`.
  So every numerical statement made "at `lambda = n/(n-1)`" is a best case.

**Cor 3.4.** At `s = 0`, `2m ln(D/2) - 2 ln m! >= -2m ln 2` for `m <= D`. With `W >= 40` this gives (3.3). The
full-rank bound `Lambda(c) >= 2 log 5/W` is correct.

**Cor 3.2 table.** All entries are reproduced, including `E(m_max,1)`:

| `n` | `E(m_max,1)` | `J(m_max)` | `J(ceil 4w)` |
| --- | --- | --- | --- |
| 2966 | 1.690 | 2 | 1 |
| 3000 | 1.685 | 2 | 1 |
| 4096 | 1.541 | 2 | 1 |
| 8192 | 1.296 | 3 | 1 |

I did not rerun `n = 16384` and `10^6`; I checked them by hand from `m/(2w)`, and they are consistent.

### 1.4 Prop 3.5 (torus picture, Philippon vacuous): [confirmed as stated], with a scope remark

* `Sigma(d) ⊃ {d s}` forces `d <= 2J` for a non-singleton `Sigma` inside a progression.
* At `T = 1` the condition `d T_0 < tau <= m_max < d` forces `T_0 = 0`. So the left side of (2.2) is
  `|Sigma| = 1`, and the margin is `ln(d! 2^d) = 22799.2`. ze.md's `22784` subtracts `ln binom(2T+d,d)`; it is
  conservative and correct.
* GV: `S_d = 17142.1` and the `log_10 W` bound of `7441.1` are reproduced (`n = 2966`).
* Scope: the premise `Sigma ⊂ <Z>` is the gap of §1.2. For the real (sumset) theorem the conclusion survives the
  automorphism points (family (a)), because `s^d` costs `d Lambda`. For light endomorphisms (family (b)) no
  argument is given.

### 1.5 Thm 3.6 and its application: theorem [confirmed]; application [overclaim, major]

The theorem is correct linear algebra: `binom(tau+d-1,d)` conditions per point, and order at `e` built into
`V(A,mu)`. It is also monotone in `mu`, so checking the endpoints plus the upward search does cover "every
`mu <= m_max + 1`" at the four `n` that were run. I reproduced every number in the ze.md table exactly.

The application has three defects.

1. **Profile restriction hidden in the summary.** The numbers are for `lambda = n/(n-1)` (the extremal value).
   §0 item 4 and the task summary state the barrier with no profile qualifier. My threshold scan (float search,
   then exact checks at `0.99` and `1.01` times the threshold) shows that at `mu = m_max` the count (3.5) holds iff
   `lambda >~ lambda*`:

   | `n` | `lambda*` at `m_max` | `lambda*` at `ceil(4w)` |
   | --- | --- | --- |
   | 2966 | 0.01684 | `<= 0.001` |
   | 3000 | 0.01978 | `<= 0.001` |
   | 4096 | 0.00943 | `<= 0.001` |
   | 8192 | 0.00613 | `<= 0.001` |

   For `lambda = 0` (exponent vectors not spanning `1^perp` with positive weight) only the `log 5` bound
   `J < 1.07 dTW/log 5` is available, and `tau_j >= mu + 1`. At `m_max` each `j` therefore contributes
   `>= 2 binom(m_max+d-1,d) ~ 3^d/4`. Already `J = 2` gives `>= 3^d/2 > avail ~ 3 * 3^d/8`, so the count fails. The barrier is therefore proved
   only for `lambda >= lambda*(n, mu)`. This is a **template-dependent threshold**, and the summary presents it as
   uniform.
2. **Omitted certified data (§1.2).** Add the `AGL(k,2)` orbit for the Sylvester profile, each orbit point with
   ze.md's own `tau_j`:

   | `n` | `mu` | `ln(count x |AGL|)` | `ln avail` | holds |
   | --- | --- | --- | --- | --- |
   | 4096 | 2256 (`ceil 4w`) | 4197.7 | 4498.1 | yes |
   | 4096 | 2631 (`m_max`) | 4565.8 | 4498.0 | **no** |
   | 8192 | 3302 (`ceil 4w`) | 6967.1 | 8998.0 | yes |
   | 8192 | 5260 (`m_max`) | 9067.6 | 8997.9 | **no** |

   The largest `mu` for which the orbit-inclusive count holds is `2559` at `n = 4096` and `5186` at `n = 8192`
   (`m_max = 2631`, `5260`).

   A subtler point: ze.md bounds the certified orders independently of `c`, through the spectral cap. Under that
   convention **every** permutation `sigma` has the same cap, since `sigma G sigma^T` has the same spectrum. So a
   `c`-independent count would have to include up to `n!` points, and it fails at every `mu` in the range. A
   correct barrier must be stated per form: for each `c`, its certified set `D(c)` (endomorphisms `M` with
   `Lambda_M(c) < ord(c) - E`) has count `< dim V(A, mu)`. ze.md neither states nor proves this.
3. **Small-value estimates.** The witness `c` lies in `V(A,mu) ⊗ C`. At a transcendental `Y` it is not integral, so
   no height is controlled. Small-value estimates (Roy, Laurent–Roy) take integer polynomials of bounded height as
   input, so Thm 3.6 says nothing about them. The phrase "zero, multiplicity **or small-value** estimate, with any
   constants" overclaims. It is correct for geometric zero and multiplicity estimates only.

Also, "`n = 2966..8192`" means the four values `{2966, 3000, 4096, 8192}`, not the whole range. No Hadamard matrix
exists at `n = 2966 ≡ 2 (mod 4)`, so "Sylvester" there means only "the parameter `lambda = n/(n-1)`".

### 1.6 Dual picture (§4.2): [error, major] in the link to ZE; the propositions themselves [confirmed]

* Lemma 4.1, Lemma 4.2 and Prop 4.3 are correct as statements about `Bs(V(A,m))`.
* Prop 4.4's numbers are reproduced: `22370.7` and `19113.3` at `n = 2966`, `T_0 = 0`.

The error is this. Route 5's rigidity gives only `Z in Bs(V_H)`. Here `V_H` is the span of the integer
`c in V(A,m')` with `|c|_inf <= H = U^2` (R12 README l. 96). The record guarantees only `dim V_H >= |A|/2` ([D]
l. 142, 287), while `dim V(A,m) >= |A| - binom(m+d-1,d) ~ (3/4)|A|` at `m_max`. So in general
`Bs(V(A,m)) ⊊ Bs(V_H)`, and:

* The dual-picture hypothesis "`Q = X_{d+1} - p(X)` vanishes on `phi_Z(A)`" is equivalent to `Z in Bs(V(A,m))`.
  Route 5 does **not** know this; ze.md's own Lemma 4.1 says only `⇐` transfers to `V_H`. A dual zero estimate
  therefore cannot even be fed the actual point.
* ZE needs `Z notin Bs(V_H)`, which is **stronger** than `Z notin Bs(V(A,m))`. So the chain in §4.2 (l. 548),
  "ZE-tor^glob ⇒ ZE^all, hence ZE", breaks at the second arrow. The same holds for §0 item 2, §6.2 and the map
  statement in §6.4. The arrow needs `V_H ⊗ Q = V(A,m) ⊗ Q`. That is not on record, and Siegel/Bombieri–Vaaler
  heights for a basis of `V(A,m)` are far beyond the Liouville budget.
* The correct target is **ZE^all_H(A,m): `Bs(V_H)` contains no nondegenerate point**. It is strictly stronger than
  ZE^all(A,m), and it is not group-structured.

The negative conclusion (the dual picture cannot close ZE) survives a fortiori. The map statement must be
corrected.

### 1.7 Prop 4.5 (HS false): [confirmed for the stated `(n,m)`], one small gap

`ref_ze/hs_check.py` reruns every `e = 0..d`, `L = 3`, in exact integers:

| `n` | `m` | tightest `ln(LHS/RHS)` (at `e = d-1`, `V x {1}`) |
| --- | --- | --- |
| 2966 | 1887 | 0.2008 |
| 2966 | 1906 | 0.1969 |
| 4096 | 2256 | 0.2545 |
| 4096 | 2631 | 0.1971 |

All hold. The tight case uses "`<k,a_p>` takes `L + 1 >= 4` distinct values on `A`". `h_A = L` alone does not give
this, and with only 3 values the `e = d-1` case fails (ratio `1.5d/(m-1+d) < 1`). For sign columns with two 1's
and two 0's off position `n`, `k = e_{x1} + e_{x2} - e_{y1} - e_{y2}` and similar vectors give the values
`-2..2`, which repairs it. The general `L >= 3` statement needs this added. HS is a hypothetical strengthening
under GV hypotheses, so its falsity only shows that a lossless count is unavailable.

### 1.8 §5 (actual point): [confirmed] up to two constant/proof slips

* **Prop 5.1.** The proof says `|<v,delta>| < 2pi`. That is false for large `v = i N_0 u`, since
  `|<v,delta>| <= |v|_1 Theta/2` is unbounded in `v`. The conclusion is still true: if `Z^{v} = 1` with
  `v' in Z^n_0`, unique factorisation gives `<v', a_p> = 0` for all `p`, so `v' = 0` by affine independence. Minor.
* **Prop 5.2.** From `|<u,delta'>| >= 2 e^{-X_u/2}` ([D] Prop 5.3(iii) proof) and `|delta'|_inf <= Theta` (these
  are differences), the constant is `2/(C'|u|_2)`, not `4/(C'|u|_2)`. Minor; the "`/poly`" reading is unaffected.
* **Uniformity.** Prop 5.1 holds for every actual cluster on a full-rank subcluster, uniformly in `R`. Prop 5.2 is
  uniform with explicit constants. Neither is used to claim growth. Correct.

### 1.9 §6: [overclaim] in the verdict wording

* The "factor `d!`" deficit is the degree normalisation (`d!(2T)^d` against `(2T+1)^d` lattice points) at `T = 1`.
  It is not intrinsic to the data. At the Hilbert level, the honest single-form margin is Thm 3.6's: `33` nats at
  `n = 2966`, `mu = m_max`. Once symmetric data are included it can turn negative (§1.5). The verdict "short by
  `~ d!`" is accurate for Philippon's theorem and misleading as a description of the obstruction.
* §6.3 "Other points ... So `Sigma ⊂ <Z>` is forced" is false (§1.2). The parenthetical "heuristic" label covers
  only "all conceivable arithmetic inputs". The specific sentence "Within one cluster these are the powers `z^j`"
  is wrong: `z_{sigma(x)}` and `prod_y z_y^{M_xy}` are tuples within one cluster with a common norm and tiny
  relative angles.
* The ZE-col\* transfer is labelled "argued, not recomputed". That is accepted as labelled, but it inherits §1.2
  and §1.6.

## 2. Independent recomputations

Scripts were written from scratch; ze.md's scripts were read only to identify conventions.
* `ref_ze/recompute.py` → `ref_ze/recompute_out.txt`.
* `ref_ze/hs_check.py` → `ref_ze/hs_check_out.txt`.

Exact big-integer binomials are used for every holds/fails verdict. `lgamma` is used only for the displayed logs
and for the `lambda*` bisection, which was then confirmed exactly on both sides.

| quantity | ze.md | referee | status |
| --- | --- | --- | --- |
| `m_max` (`8 binom(m+d-1,d) <= 3^d`), `n = 2966/3000/4096/8192` | 1906/1928/2631/5260 | same | exact, agree |
| `w`, `ceil(4w)`, `n = 2966` | 471.7, 1887 | 471.72, 1887 | agree |
| `E(m_max,1)`, `J(m_max)`, `J(ceil 4w)` | Cor 3.2 table | §1.3 table | agree |
| Prop 3.5(1) margin, `n = 2966` | 22784 | 22799.2 (`T_0 = 0`); 22783.9 with ze.md's `binom(2T+d,d)` | agree |
| GV `S_d`, `log_10 W` bound, `n = 2966` | 17142, 7441 | 17142.1, 7441.1 | agree |
| Prop 4.4, `n = 2966` | 22370, 19113 | 22370.7, 19113.3 | agree |
| Thm 3.6 `ln` count / `ln` avail, all 8 cells | table | identical to 0.1 | agree |
| Thm 3.6 largest `mu` | 1907/1929/2633/5262 | not rerun (monotone; endpoints agree) | — |
| Prop 4.5 tightest margin | 0.197 | 0.1969 (`n = 2966`, `m_max`), 0.1971 (`n = 4096`, `m_max`) | exact, agree |
| **new:** `lambda*` for Thm 3.6 at `m_max` | not stated | 0.01684 / 0.01978 / 0.00943 / 0.00613 | exact either side |
| **new:** Thm 3.6 with the `AGL(k,2)` orbit, `m_max` | not considered | fails at `n = 4096, 8192`; holds up to `mu = 2559, 5186` | exact |

## 3. Literature check of the quoted zero estimate

* **Bibliography confirmed** by web search: P. Philippon, *Lemmes de zéros dans les groupes algébriques
  commutatifs*, Bull. SMF 114 (1986) 355–383, DOI 10.24033/bsmf.2060 (numdam; eudml doc 87515), with *Errata et
  addenda* in Bull. SMF 115 (1987). The search also surfaced the Nakamaye-refined form and Philippon's *Nouveaux
  lemmes de zéros* (1996).
* **Statement not re-read.** In this session too, arxiv.org and eudml.org fail DNS (network policy). So the exact
  normalisation is **unverified by me as well**. To read the source, numdam.org, eudml.org or arxiv.org has to be
  allowed in the session's network settings.
* **From memory**, the shape quoted in ze.md (2.1) agrees with the form in Waldschmidt, *Diophantine approximation
  on linear algebraic groups* (2000), Thm 5.1:
  * vanishing to order `> nT` along `W` on `Sigma(n)` (`n = dim G`);
  * conclusion `binom(T+l*, l*) card((Sigma+G*)/G*) H(G*;D) <= H(G;D)`;
  * `H(G_a^{d_0} x G_m^{d_1}; D_0, D_1) = ((d_0+d_1)!/d_0!) D_0^{d_0} D_1^{d_1}`.

  That gives `d! D^d` and `(d+1) D_0^d D_1` in the two cases used, matching ze.md.
* **Robustness.** Whatever the exact constants are, the vacuity results are robust. They need only (2.2), and the
  margins are `>= e^{1.7*10^4}`, as Remark 2.2 says. GV (vanishing on `Sigma`, order `> T`) also covers the
  Nakamaye/Philippon-1996 type refinements I am aware of, which weaken the sumset hypothesis.
* **Not robust.** The Hilbert-level statements (Thm 3.6, Prop 4.5) do not depend on the quotation at all. Their
  thin margins (33–56 nats) are the real quantitative content, and they are affected by §1.2 and §1.5.

## 4. Issues

**Major**

* **M1 (§1.2).** "Route 5 certifies vanishing only at `Z^j`" is false. Lemma 3.1's proof certifies `P_c(Z^M) = 0`
  for every integer endomorphism `M` with `Lambda_M(c) < ord(c) - E`. This includes profile automorphisms at cost
  exactly `Lambda(c)` (the Sylvester profile has `|Aut| >= |AGL(k,2)| ~ e^{107}` at `n = 4096`), and possibly
  light endomorphisms. §6.3's "`Sigma ⊂ <Z>` is forced" is wrong.
* **M2 (§1.5).** Thm 3.6's application ("every point carries a witness with all route-5 data, every
  `mu <= m_max + 1`") is not established.
  * With the `AGL` orbit at ze.md's own orders, the count fails at `mu = m_max` for `n = 4096, 8192` (exact). It
    holds only up to `mu = 2559` and `5186` respectively.
  * Under ze.md's `c`-independent convention, all `n!` permutations would have to be counted.
  * A per-form statement is required, and it is absent.
* **M3 (§1.5).** The single-form barrier is proved only for `lambda >= lambda*(n,mu)` (e.g. `0.0168` at `n = 2966`,
  `m_max`), yet it is stated without that qualifier in §0 item 4 and §5 item 4(b). This is a template-dependent
  threshold presented as uniform, which the owner's constraint forbids.
* **M4 (§1.6).** The dual picture is built on `Bs(V(A,m))`, but route 5 knows only `Z in Bs(V_H)`, and the record
  has only `dim V_H >= |A|/2`. "ZE^all(A,m) ⇒ ZE" (l. 548, §0 item 2, §6.2, §6.4) is unproved. The relevant target
  is ZE^all_H, which is stronger and not group-structured.

**Minor**

* **m1.** The phrase "small-value estimate" in Thm 3.6's reading: the witness is complex, so no height control.
* **m2.** Prop 4.5 needs "`>= L+1` distinct values of `<k,a_p>` on `A`". It is repaired for sign columns with two
  0's and two 1's off `n`.
* **m3.** The proof of Prop 5.1 uses a false bound `|<v,delta>| < 2pi`. The conclusion holds by unique
  factorisation.
* **m4.** Prop 5.2's constant should be `2/(C'|u|_2)`, not `4/(C'|u|_2)`.
* **m5.** "`n = 2966..8192`" is four values. "Sylvester" at `n = 2966` is only the parameter `lambda = n/(n-1)`;
  no Hadamard matrix exists there.
* **m6.** The "`~ d!` short" headline describes Philippon's degree normalisation, not the intrinsic (Hilbert-level)
  deficit, which is 33–56 nats and fragile (M2).
* **m7.** Philippon's statement is still a memory-quote. The vacuity conclusions are robust to this; nothing else
  depends on it.

**Confirmed (no issue):**
* Lemma 3.1, including no discarded gcd costs and the common phase cancelling.
* Prop 3.3 and Cor 3.4.
* The Cor 3.2 table.
* Prop 3.5 as stated, with every margin reproduced.
* Lemma 4.1, Lemma 4.2, Prop 4.3 (as statements about `V(A,m)`) and Prop 4.4.
* The Prop 4.5 numerics.
* The Prop 5.1 conclusion.
* The absence of any growth or uniformity claim. I found no step anywhere that asserts `M = o(log R/loglog R)` or a
  uniform bound. ze.md correctly claims neither.

## 5. Corrected statement

> **Round 26 (ZE via zero estimates), corrected.** No growth improvement and no uniform bound;
> `M <= (2/3+eps) log R/loglog R` stands. ZE-col\*, MIN^pair and `(R_H)` are not derived.
>
> **(1) Proved.**
> * Route 5's Liouville step certifies, for a form `c`, the vanishing of `P_c` at `Z^M` for every integer
>   endomorphism `M` of the cluster torus with `Lambda_M(c) < ord(c) - O(E)`. This includes the powers `Z^j`
>   (`|j| < ord(c)/Lambda(c) <= 1.07 sqrt n/lambda`, by the new spectral cap
>   `ord(c) <= (sqrt n/lambda) Lambda(supp c)`) and the profile-automorphism images `Z^{j sigma}`.
> * On `T = 1` regime boxes the powers guaranteed for every `c in V_H` are `j = +-1` at `ceil(4w)` and `|j| <= 2..4`
>   at `m_max`.
>
> **(2) Proved, given the shape of Philippon's theorem** (bibliography verified; statement a memory-quote, with
> margins `>= e^{1.7*10^4}`).
> * On `G_m^d`, with data in `<Z>`, Philippon's conclusion is satisfied by the trivial subgroup whenever
>   `lambda > 2.14 sqrt n/d`.
> * A generous variant without the sumset loss is vacuous for every full-rank profile with
>   `ln W <~ d ln(d/11)`.
> * Profile automorphism points do not change this, because their `d`-th powers are not certified.
> * Light endomorphisms are not treated.
>
> **(3) Hilbert-level single-form barrier, proved only in this restricted form.**
> * For every `Y`, some `c in V(A,mu) ⊗ C` carries the vanishing data at the **powers** `Y^j` with ze.md's
>   `c`-independent orders. This holds exactly at `n in {2966, 3000, 4096, 8192}` for `mu <= m_max + 1` when
>   `lambda = n/(n-1)`, and at `mu = m_max` iff `lambda >= lambda*` (`0.0168`, `0.0198`, `0.0094`, `0.0061`).
> * Including the `AGL(k,2)` orbit of the Sylvester profile, it holds only for `mu <= 2559` (`n = 4096`) and
>   `mu <= 5186` (`n = 8192`).
> * A barrier covering **all** route-5 data needs a per-form count over the certified endomorphisms `D(c)`. That
>   count is open.
> * The barrier says nothing about small-value estimates, which need integral polynomials.
>
> **(4) Proved.** In the dual picture (`G_a^d x G_m`) the obstruction subgroups are `V x {1}` and `V x G_m`, and
> for nondegenerate points the counts are point-independent. The degree-level estimate is vacuous, and the
> Hilbert-sharp GV form is false (exact at `n = 2966, 4096`). However, the dual picture concerns
> `Bs(V(A,m)) ⊂ Bs(V_H)`. Route 5 gives only `Z in Bs(V_H)` and needs `Z notin Bs(V_H)`. So ZE^all(A,m) and
> ZE-tor^glob do **not** imply ZE without `V_H ⊗ Q = V(A,m) ⊗ Q`, which is not on record. The relevant global
> target is ZE^all_H.
>
> **(5) Proved.** `Z` on a full-rank subcluster lies on no proper algebraic subgroup (by unique factorisation).
> Its relative distance to `K_U` is `>= (2/(C'|u|_2)) e^{(W/2 - X_u)/2}` for each defining character `u`.
>
> **Verdict.** The two natural zero-estimate framings of route 5 (one form on `G_m^d` with data in `<Z>`; one
> exponential polynomial on `G_a^d x G_m`) are proved useless: vacuous at the degree level, with the trivial
> subgroup as obstruction. The general claim that *no* single-form zero/multiplicity estimate can help is proved
> only for powers-data and `lambda >= lambda*`, with thin Hilbert-level margins (33–56 nats) that symmetric
> profiles overturn at `mu` near `m_max`. Status: **obstruction_found (partial)**. The open lemma is unchanged.
