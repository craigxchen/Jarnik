# Referee report on `round16/antitheory.md` (route 5, `T = 1`, polynomial anti-concentration versus AC(u))

Setting of every statement below, unless a line says otherwise: `T = 1`, `d = n - 1`, `w = sqrt(8 d ln 4r)`,
`j = floor(w/2)`, `m = m' = ceil(4w)`, regime `d >= 316 ln 4r`, scale `delta_0 = d^m eps_0`,
`eps_0 = C d e^{-W/4}/(m+1)`, `abar = W/(2d)` (`C = C' = 1` wherever numbers are given). Profiles: arbitrary unless
stated; "Sylvester" means `n = 2^k`, `r = n - 1`. All my numbers come from my own scripts in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_antitheory_checks/`
(50-digit mpmath arithmetic, exact integer central trinomial coefficients, exact integer `c` in the tail tests). I did
not import or run the author's code. I only *read* the author's `act_bits_theta.py` header and output, to find where a
discrepancy came from (§1.3).

## 0. Verdict

**The write-up survives as a negative/structural report. It contains no growth claim and no uniformity claim, and
none is warranted.** `M <= (2/3 + eps) log R/loglog R` stands. `AC_A` on `G` (MIN^pair at `H = infinity`) and
MIN^pair itself remain open.

Every proved statement I re-derived is correct up to the corrections below: Prop 1.2, Remark 1.3's inequality,
Theorem 3.1, Lemma 4.1, Prop 4.3, Theorem 4.4, Lemma 5.1, Cor 5.2 and Lemma 5.3. Every numerical table I recomputed
agrees with mine to the printed digits: §3.1, Thm 3.1(4), Prop 4.2, Thm 4.4 including the `W_min` negative, and
Table 6.1.

No fatal issue. There is one **major** issue, a headline overreach. "A true anti-concentration statement must see the
slabs" (summary item 3, §4.4) is shown only at finite `n`, and for `G` not at all. My own extension (§2.3) pushes the
`(MIN)`-region failure of `AC_T(.;1/2)` beyond the author's range: it holds at `k = 26, 30, 31`. It then **stops at
`k = 32`**, and the necessary condition for `G` dies at `k = 36`. So the failure is a finite-`n` phenomenon of this
family, and that is consistent with the author's own asymptotic remark, which contradicts the headline.

The minor issues are:
* a wrong closed form for `theta_0` (the code is right);
* "never at `theta_S`", which the author's own output contradicts;
* Prop 4.2 quantified with the wrong inequality;
* "maximally regular" at `b = 1`;
* "strictly stronger" without a regime witness;
* the fake-test remark on ILO_m names the wrong tangent cone (`TC_H` instead of `TC_V`);
* ILO_m's MIN^pair implication omits the 32-class projection-cheap spans, which push the exponent-`m` threshold to
  `n ~ 2.2 * 10^6`;
* an unproved meta-claim about how anti-concentration theorems are proved.

## 1. Statement-by-statement

### 1.1 Prop 1.2 (exact restatement, dual form): correct

* (1) This is nearstrata Prop 2.1 with `V = V(A,m)`. The infimum is attained (finite dimension), so the strict
  inequalities match.
* (2) and (3) are immediate. I checked against `local_smoothing/steelman.md` S3: the steelman's `AC(u)` is
  `Pr_{k uniform on T}(|Q| <= d^m eps_0) < 1/2`, a measure on the **full** cube. So `AC = AC_T(.;1/2)` exactly as
  the author says.
* (4) `1 in Poly_{<m}` forces `sum c = 0`, and `kappa_c d^m = |Phi_c|/|c|_1 = |E_{nu+} l_u^m - E_{nu-} l_u^m|/2`.
  The normalisation matches the record (`kappa_c = |Phi_c|/(|c|_1 d^m)`, README §1). Correct.
* **"Strictly stronger" (summary item 1, §4.4) is not proved for any regime box.** To show `AC_T(.;1/2)` strictly
  stronger than `AC_A` one needs a `u` where `AC_A` holds and `AC_T` fails. At every failing point of Thm 4.4,
  `AC_A` is open (that is MIN itself). A witness does exist outside the regime:
  * Take the uncut cube `A = Cube` and `u = (1, ..., 1)`.
  * nearstrata Thm 2.4 gives `kappa > 0` there, since `nu(u) = 2 |supp u| = 2d >= m`.
  * `Q = S^{m-2c-1} prod_{|t|<=c} (S - t)`, with `S = k([d])` and `c ~ sqrt(d)`, vanishes on half the cube.

  So the correct wording is: "strictly stronger as a statement form (witness: the uncut cube). For regime boxes the
  converse is merely unproved."

### 1.2 Remark 1.3 (the `H`-gap): inequality correct, title overstated

* I re-derived the chain:
  1. `Lambda = V cap Z^A` is primitive.
  2. `covol(Lambda) = covol(V^perp cap Z^A)`.
  3. That is at most the covolume of the integer row lattice, which Hadamard bounds by `U^{E/2}`.
  4. Vaaler gives `vol_N(V cap [-1,1]^A) >= 2^N`.
  5. Minkowski II then gives `prod lambda_i <= covol`.
  6. Since `lambda_i >= 1`, `#{lambda_i > H} <= (E/2) ln U/ln H`.

  All correct.
* The rigidity tolerance `ln|c|_1 <~ wW/2 - m ln(dC') + ln m!` is exactly nearstrata Lemma 1.1's inequality
  `(m'/2 - w) W <= 2 ln|c|_1 + 2m' ln(dC') - 2 ln m'!` at `m' = 4w`. So "H up to `exp(Theta(n^{3/2} ln^{3/2} n))`"
  is right (one should add "for `W > W_1`").
* **Overstatement.** The remark is titled "the `H`-gap is real", but it proves only that *this* count leaves a
  codimension bound of size `E^{1-o(1)}`. It proves neither `V_H != V` nor `kappa_H < kappa_V` anywhere.
* **A sharper exact statement**, which I add and which the author may use. The maximum of the linear functional
  `c -> Phi_c(u)` on `V cap {|c|_1 <= 1}` is attained at a vertex. The vertices of a subspace cut with the
  cross-polytope are normalised elementary vectors (circuits) of `V = ker M`. A circuit has an integer representative
  whose entries are `rank(M) x rank(M)` minors of a `{-1,0,1}` matrix, so by Hadamard
  `|c|_inf <= E^{E/2}` and `|c|_1 <= (E+1) E^{E/2}`. Hence

  ```text
  kappa_V(u) = kappa_{H*}(u)   for every u,   H* = E^{E/2},   ln H* = (E/2) ln E .
  ```

  The `H`-gap is therefore exactly the window `ln H in [2 ln U, (E/2) ln E]` against a rigidity tolerance of
  `~ wW/2 = Theta(n^{3/2} ln^{3/2} n)`. Since `ln E = Theta(m ln(d/m))` is superpolynomial in `n`'s logarithm
  (`E >= e^{c sqrt n}`), the circuit bound does not close it either. This confirms the author's conclusion from a
  different side.

### 1.3 Remark 1.4, §2 table, §3.1 (Carbery-Wright proxy)

* §3.1 table: reproduced exactly with `sigma = sqrt(m!)` (`tables_out.txt`):
  `+3.639, +3.347, +1.982, +3.129, -6.364, -59.02`.
* With the ternary-proxy variance `sigma >= (2/3)^{m/2} sqrt(m!)` every entry rises by `+0.088` in log10, so the
  crossover is unchanged. Evidence about a proxy only, as the author says.

### 1.4 Theorem 3.1 / Cor 3.2 (degree barrier): correct; one framing issue

* (1)-(3) re-derived. `Var_{k_x} Q` has mean `(2/3)/(2b/3) = 1/b`. `mu{Q != 0} = (1 - p_b)^m` by independence.
* `p_b >= 0.4/sqrt b` holds for all `2 <= b <= 3000`. I checked it with exact central trinomial coefficients, via the
  recurrence `(b+1)T_{b+1} = (2b+1)T_b + 3b T_{b-1}`, which I validated on `1, 1, 3, 7, 19, 51`.
* Table (4) reproduced to all printed digits: `-328.23, -397.26, -648.54, -955.07, -1105.44, -2088.47`. One row was
  not printed: `n = 10^4`, `r = n^2` gives `b = 1` and `-886.62`.
* **Framing (minor).** At `n = 2900` and `n = 4096` we have `b = floor(d/m) = 1`. There `Q` is a *single monomial*:
  MNV rank 1, every influence equal to `1 = Var Q`. Calling it "maximally regular" is true only relative to the lower
  bound `max_x Inf_x >= m/d` that every homogeneous top-chaos degree-`m` polynomial obeys, and `m/d = 0.64` at
  `n = 2900`. So at the record's thresholds no degree-`m` polynomial is regular in the MOO sense at all. The barrier
  is then vacuous rather than informative: the hypotheses it "defeats" cannot even be met. The rows with `b >= 2`
  (`n >= 10^4`) carry the content.
* Cor 3.2 as written, with `F(m, 1/b, 0, b, 0)`, is correct.

### 1.5 Lemma 4.1, Props 4.2-4.3, Theorem 4.4

* **Lemma 4.1:** re-derived and correct. The interpolants `R_i o pi` have degree `<= rho`, the subtracted terms have
  degree `<= rho + i < m`, and `|Q| <= d^m sum_{i>=o} binom(m,i) eta^i <= 2^m eta^o d^m`.
* **Prop 4.2: wrong inequality in the hypothesis.** "Disjoint subsets of size `>= N`" must read "`|X_i| <= N`" (the
  script uses `|X_i| = floor(d/q)` exactly). The tail `2 exp(-3t^2/(4M))` *increases* with `M`, so a lower bound on
  `|X_i|` does not give the stated `mu(Y)`.
  * Counterexample to the literal statement: `q = 1`, `X_1 = [d]`, `N = 1`, `c_1 = 2`. The hypothesis
    `1 - 2e^{-27/4} >= 1/2` holds, but `mu(|k([d])| <= 2) -> 0`.
  * The *conclusion* happens to hold there by another choice of `c`.
  * The table (`q_max = 49, 49, 52, 89, 58, 69, 86`) is reproduced exactly with `|X_i| = floor(d/q)`; I also get
    `60` at `n = 10^5` and `94` at `n = 10^9`.
* **Prop 4.3:** correct. `X_i` avoids the base point `0`, `|X_i| = n/2`, and the tail is
  `2 exp(-3(c+1)^2/(2n))`. The sub-Gaussian MGF bound `(1 + 2cosh s)/3 <= e^{s^2/3}` holds termwise
  (`2 * 3^{i-1} i! <= (2i)!`).
* **Wording, l. 265:** "polynomials that vanish on half of the cube for collision-free `u`" is wrong. Exact vanishing
  holds only on `L_q` for admissible `q`, and admissible `q < k` (e.g. `q <= 11` at `k = 12`). The `L_q` with
  `q < k` have classes of size `>= 2`, i.e. collisions. At collision-free `u` one gets only `|Q| <= delta_0`.
* **`theta_0` has the wrong closed form (minor; code right, text wrong).**
  * The text has `theta_0 = 1 - (1/4)(4r)^{-1/2}` (l. 155, l. 346).
  * The sub-Gaussian union bound gives `r * 2 exp(-3w^2/(16d)) = 2r (4r)^{-3/2} = (1/2)(4r)^{-1/2}`.
  * The author's script uses `beta0 = 2r(4r)^{-3/2}`, so the numbers are right.
  * With the text's value, `AC_T(theta_0)` would not even be implied by a proved bound on `mu(A)`.
  * My recomputation (`bits_t0_out.txt`) with the **correct** value reproduces the author's table: `k = 16` fails at
    `q = 6` (`W = 3n ln n`), and `k = 17` at `q = 6` (`W = 4n ln n`). With the **text's** value, `k = 16` and
    `k = 17` would *not* fail.
* **Theorem 4.4:** reproduced independently. With `theta = 1/2` and `W = 3n ln n`, the failing `q` are:

  | `k` | failing `q` |
  | --- | --- |
  | 12 | 2..7 |
  | 13 | 3..8 |
  | 14 | 3..9 |
  | 15 | 4..9 |
  | 16 | 4..10 |
  | 17 | 5..11 |
  | 18 | 6..12 |
  | 19 | 6..12 |
  | 20 | 7..13 |
  | 21 | 7..13 |
  | 22 | 8..13 |

  At `theta_S`: `q_max = 4` for `k = 12, 16, 22, 30, 60`. At `W = W_min`: no failure, for `k = 12, 14, 16, 20, 22`
  and both `theta = 1/2` and `theta_0`.

  I also checked the geometric inputs:
  * `|u|_inf in [1, 1+eta]`;
  * within-class gaps `eta |Delta rank|/s`;
  * cross-class gaps `>= 1/(K-1) - eta`;
  * the `ln^+` bound (all bounds `>= 1`, so `ln^+ = ln`);
  * pairs with the base point included (class 0 contains the base point, and `N_2 = binom(n,2)`).
* **"never at `theta_S`" is false.** The author says (l. 388) that the necessary condition `K o_q > d` holds "never at
  `theta_S`". The author's own `act_bits_theta_out.txt` prints `K*o>d at q in [2, 3, 4]` for `k = 12` and `[3]` for
  `k = 13`. I get the same (`necK_out.txt`): at `k = 12`, `q = 4`, `c = 246`, `o = 288`, `K o = 4608 > 4095`.
  * It is harmless: these are 2-4-column bit strata, i.e. the record's own cheap range, where MIN\* is already known
    to fail at `n = 4096`.
  * The sentence should read "only for `k <= 13`".
* The other asymptotic claims check out: `K o_q > d` for some `q` holds exactly up to `k = 36` (`theta = 1/2`) and
  `k = 19` (`theta_0`).

### 1.6 §4.4 and summary item 3: the headline lesson overreaches (MAJOR)

The summary says "So a true anti-concentration statement must see the slabs"; §4.4 adds "A profile-free statement of
the steelman type is false at regime sizes". What is actually proved:
1. The failures are at finite `n` only. By the author's own asymptotic test, bit strata satisfy the averaged
   necessary condition only for `k <= 36` (`theta = 1/2`) and `k <= 19` (`theta_0`). Partition strata are
   collision-rich. Growth needs only `n >= n_1`, and `n_1 > 4096` is already forced (MIN\*). So a profile-free
   `AC_T(.;1/2)` on `G`, or even on the `(MIN)` region for `n >= 2^37`, is **not refuted** by anything in the write-up.
2. The failing points are in `{Gamma <= abar}`, not in `G`. For genuine Sylvester pair data `G` is empty.
3. They need `W = 3n ln n`, not `W_min`.

The correct lesson is narrower: "at the finite sizes where route 5 would have to start (`n = 2^12 .. 2^36` at
`theta = 1/2`), a density form that ignores the slabs is false on the averaged region". That is still useful,
because it pushes any density-form `n_1` past `2^36`. But it does not show that the slabs must be seen
asymptotically. I extended the `(MIN)`-region test beyond the author's `k <= 22` (§2.3).

### 1.7 §5: Lemma 5.1, Cor 5.2, Lemma 5.3, §5.3, ILO_m

* **Lemma 5.1:** correct. On `F_t`, `l_y = <tau, t>` is constant, and `l_u^m = eta^m l_v^m mod Poly_{<m}|_{F_t}`.
  Restriction only lowers the distance, and zero extension keeps entries, so the height carries over to `V_H(A)`.
  One addition the author should make:
  * the fibre quantity depends only on `v mod L` (`(l_v + const)^m = l_v^m mod Poly_{<m}` on `F_t`);
  * so the decomposition should take `eta = dist_inf(u, L)` and `v` with `dist_inf(v, L) = 1`.

  Otherwise "regular relative to `F`" is ill-posed.
* **Cor 5.2:** a correct conditional. The loss `m(Q_F abar/2 + ln 2)` uses the per-pair separation lemma, which holds
  for every `u in G`, not only for the actual point (its proof uses only the per-pair gaps and
  `sum eps <= N_2 abar`). The "`n ~ 10^6`" threshold is my crossover `n ~ 5.27 * 10^5` (`r = n`, `W_min`, `Q = 16`,
  exponent `m`; `ledgerQ_out.txt`).
* **Lemma 5.3:** correct. `floor((m-1)/(2j)) = 4` unless `w` is an even integer, when it is 3; for irrational `w` it
  is 4. The mixed difference `prod Delta^{alpha_p}` has `|c|_1 = 2^m` and value `m! prod z_p^{alpha_p}`.
* **§5.2's meta-claim (minor).** "Every located anti-concentration theorem is proved by (i) differencing ... or (ii)
  continuous comparison" is unproved and not accurate as stated:
  * Halasz/Esseen-type bounds (Nguyen-Vu) are Fourier-analytic;
  * Kwan-Sauermann's quadratic theorem is algebraic/combinatorial;
  * Fox-Kwan-Sauermann is Sperner/LYM-type.

  It should be labelled a heuristic reading. It does not affect anything proved.
* **ILO_m (§5.4), three corrections.**
  1. *Fake test, wrong tangent cone.* ILO_m is a `kappa_V` (`H = infinity`) statement. It is refuted by
     nondegenerate real points of `TC_V = {kappa_V = 0}`, not by points of `TC_H` as the author says (l. 519-520).
     A nondegenerate point of `TC_H \ TC_V` is a fake obstruction to the *route*, but not a counterexample to ILO_m.
     There ILO_m holds and the `H`-gap of Remark 1.3 is what fails. So even a proof of ILO_m for all profiles would
     not by itself decide the barrier question for route 5; the "equivalent in difficulty" phrasing needs this
     qualification.
  2. *ILO_m quantifies over "every regime profile".* Fake pattern profiles (sign columns, slab half-width `j`) are
     regime profiles. So the conjecture as stated predicts that no regime (Mar) fake box has a nondegenerate real
     `TC_V` point. That should be said explicitly. The alternative is to state ILO_m only for actual profiles, but
     then it is not a box statement and fake Thm 2.1 says such a statement cannot be distinguished by box data alone.
  3. *The implication to MIN^pair needs all projection-cheap `L`, and the known ones already move the threshold.*
     Round 13 put 32-class projection-cheap spans `L_{2,1}` in `TC` for Sylvester `k >= 5`. At exponent `m` (forced by
     ILO_m's `1/m`), the ledger `L0 - Lc - m(Q abar/2 + ln 2) > 0` holds (`ledgerQ_out.txt`; `r = n`, `W_min`,
     `C = 1`) only from:

     | `Q` | from `n ~` | margin at `n = 10^6` (nats) |
     | --- | --- | --- |
     | 16 | `5.3 * 10^5` | `+1.78e6` |
     | 32 | `2.2 * 10^6` | `-3.07e6` |
     | 64 | `9.4 * 10^6` | |
     | 128 | `4.0 * 10^7` | |
     | 256 | `1.7 * 10^8` | |

     So "MIN^pair (`H = infinity`) for width-cheap strata from `n ~ 10^6`" is correct as stated. But the known
     32-class projection-cheap strata need `n >~ 2.2 * 10^6` even if their collision counts are as good as
     possible.

## 2. Own computations

All scripts are mine (`referee_antitheory_checks/`). Each writes or tees an `_out.txt` file.

### 2.1 `tables.py` (`tables_out.txt`)

* Theorem 3.1 table: exact `T_b`, and `p_b >= 0.4/sqrt b` on `[2, 3000]`.
* §3.1 Carbery-Wright table: both normalisations of `sigma`.
* Table 6.1, every column: `-11650.2, -12714.0, -9759.4, +165947, +317403, +6.04e6, +1.20e10` and
  `-72051, -94380, -179223, -516341, -599266, +1.776e6, +1.17e10`.
* Prop 4.2 `q_max`.

All agree with the write-up.

### 2.2 `bits.py` (`bits_t0_out.txt`, `bits_half_out.txt`)

Theorem 4.4 recomputed from scratch, with the Gamma bound evaluated by `ln G(s+1) = sum_{u<s} ln u!` (Barnes `G`).
Everything agrees with the author once `theta_0` is taken as `1 - (1/2)(4r)^{-1/2}`, the code's value, not the
text's.

### 2.3 Extension of Theorem 4.4 (`bits_large_out.txt`)

Same construction, `theta = 1/2`, `W = 3n ln n`, `C = 1`:

| `k` | `m` | failing `q` (point in `{Gamma <= abar}`, `AC_T(.;1/2)` false) |
| --- | --- | --- |
| 26 | 408307 | 10..14 |
| 30 | 1745994 | 13..15 |
| 31 | 2507493 | 14, 15 |
| 32 | 3599458 | none |
| 33 | 5164719 | none |
| 34 | 7407624 | none |

So the author's Theorem 4.4 is true at least for `k in {12..22, 26, 30, 31}`. Here `k = 23..25` and `27..29` were
not run; that they hold is interpolation, i.e. evidence. The construction **stops working at `k = 32`**, earlier than
the necessary condition `K o_q > d`, which lasts to `k = 36` (`necK_out.txt`).

This sharpens §1.6: the density-form failure on the averaged region is a finite-`n` phenomenon of this family, and
it ends at `n = 2^32`.

### 2.4 `necK.py` (`necK_out.txt`)

* The necessary condition holds up to `k = 36` (`1/2`) and `k = 19` (`theta_0`).
* At `theta_S` it holds at `k = 12` (`q = 2, 3, 4`) and `k = 13` (`q = 3`), contradicting "never at `theta_S`".
* `q_max(theta_S) = 4` for `k = 12, 30, 60`.

### 2.5 `ledgerQ.py` (`ledgerQ_out.txt`)

Exponent-`m` crossovers for `Q`-class strata (`r = n`, `W_min`, `C = 1`): `Q = 16, 32, 64, 128, 256` at
`n ~ 5.3e5, 2.2e6, 9.4e6, 4.0e7, 1.7e8`.

## 3. Literature check (WebSearch snippets only)

arxiv.org, theoryofcomputing.org, dokk.org and drops.dagstuhl.de all failed DNS for WebFetch, and curl was refused by
the proxy (403). So, like the author, I saw abstracts and snippets only.

* **Meka-Nguyen-Vu**, Theory Comput. 12 (2016) art. 11, arXiv 1507.00829.
  * Existence, venue, multilinear setting and rank notion confirmed. Search snippets give the baseline as
    `(log r)^{O(1)}/sqrt r` for fixed degree.
  * I could **not** confirm the verbatim quote "non-trivial for degrees up to `c log r/loglog r`", nor the remark
    about degree `log_2 r`. The author labels it [L, quoted from the arXiv text] but also says arXiv was unreachable,
    so the label is inconsistent. It should be [M] or "paraphrase".
  * It does not matter: the rank is at most `d/m` for any degree-`m` polynomial, so MNV is vacuous here on any
    reading.
* **Jin-Kwan-Sauermann-Wang**, arXiv 2505.23335 ("Algebraic aspects of the polynomial Littlewood-Offord problem",
  robust irreducibility). Confirmed.
* **Grebennikov-Kwan**, arXiv 2505.24699 ("Geometric Littlewood-Offord problems via lattice point counting",
  bounded Chow rank, bounded degree). Confirmed.
* **Kwan-Sauermann**, "Resolution of the quadratic Littlewood-Offord problem", Compositio 2025. Confirmed (ISTA
  record).
* **Grebennikov 2610.08708.** I also found only the MathWorld entry. Agree: unverified.
* **Filmus-Mossel**, PTRF 175 (2019). The degree condition is `o(n)` in the journal version but `o(sqrt n)` in the
  CCC 2016 / arXiv v1 version. The author's "degree `o(n)`" is right for the journal version. Under the `o(sqrt n)`
  version our `m ~ 11.3 sqrt(d ln 4r)` is not covered at all.
* **MOO, Carbery-Wright, DSTW, Kane, Nguyen-Vu, dimension-free Remez:** not re-checked beyond the author's
  descriptions. None is used in any proof, only in the "not applicable" column. Every "not applicable" verdict is
  correct on the stated hypotheses, because of degree, domain, scale, or rank `<= d/m`.

**Conclusion on literature.** No quoted theorem is used as a premise of a proved statement, so citation inaccuracies
cannot propagate. The survey's verdict, that nothing applies at degree `m ~ sqrt(d ln r)`, scale `e^{-Theta(n ln n)}`
and the slab domain, is correct.

## 4. Fake-box test (rubric b)

* **Proved statements.** Prop 1.2, Thm 3.1, Lemma 4.1, Props 4.2-4.3, Thm 4.4, Lemma 5.1 and Lemma 5.3 are box or
  cube statements with no arithmetic input. Fakes cannot refute them, and none of them claims anti-concentration.
* **§5.3(3)** uses the column-span relation (fake-violated when profile rank `< d`) only to show that kernel
  certificates vanish at the actual point. That is correct and harmless.
* **ILO_m.** This is the only place a fake box could bite:
  * as stated ("every regime profile"), it includes fake pattern boxes, and is false on any regime (Mar) fake box
    whose `TC_V` (not merely `TC_H`) has a nondegenerate real point;
  * conversely, fakes at points of `TC_H \ TC_V` do not touch ILO_m but do defeat the route (§1.7(1)).

  The author missed this split; it is the precise sense in which ILO_m + `H`-gap, not ILO_m alone, is "as hard as
  the barrier question".
* **Theorem 4.4 and Sylvester boxes.** Theorem 4.4's density-form counterexamples are profile-free (cube) statements.
  Fake Sylvester-type boxes have `mu(A) >= theta_S` as well, so the `theta_S` analysis applies to them verbatim. No
  new fake counterexample arises there.
* I found no fake-box counterexample to any proved statement.

## 5. Corrected statement

At `T = 1`, `w = sqrt(8 d ln 4r)`, `j = floor(w/2)`, `m = ceil(4w)`, regime `d >= 316 ln 4r`, scale
`delta_0 = d^m eps_0`, `C = C' = 1` for all numbers:

1. *(Proved.)* `AC_A(u) <=> kappa_V(u) > eps_0`, with the moment-matching dual (two probability measures on `A`
   agreeing on moments `< m`).
   * `AC_T(u;theta) => AC_A(u)` for `theta <= mu(A)`. The steelman's `AC = AC_T(.;1/2)`.
   * Strictness of the density form is witnessed only outside the regime (uncut cube).
   * `kappa_V = kappa_{H*}` with `ln H* = (E/2) ln E` (circuits). Neither Minkowski-Vaaler counting nor the circuit
     bound brings the height below the rigidity tolerance `~ wW/2`. The `H`-gap is open, not proved real.
2. *(Proved.)* For `b = floor(d/m) >= 2` there are top-chaos, equal-influence (`1/b`), rank-`b` degree-`m`
   polynomials with `mu{Q != 0} <= exp(-0.4 m^{3/2}/sqrt d)` (exact table as in Thm 3.1). At `b = 1`
   (`n = 2900, 4096`) the example is a monomial, and no degree-`m` polynomial is regular at all. So degree,
   influence, chaos or rank hypotheses alone cannot give anti-concentration at degree `m`.
3. *(Proved; finite `n`.)* For Sylvester `(n, r) = (2^k, 2^k - 1)`, `W = 3n ln n`, `AC_T(.;1/2)` fails at explicit
   points of `{|u|_inf = 1, Gamma(u) <= abar}` for `k = 12..22` (author) and `k = 26, 30, 31` (referee). The
   construction fails for `k = 32, 33, 34`.
   * `AC_T(.;theta_0)` with `theta_0 = 1 - (1/2)(4r)^{-1/2}` fails for `k = 12..16` (and `k = 17` at `W = 4n ln n`).
   * At `theta_S` the exceptional bit strata are the `<= 4`-column cheap range for `k <= 60`.
   * Nothing is shown at `W_min`, nothing on `G`, and nothing for `k >= 32` (`theta = 1/2`) or `k >= 37` (necessary
     condition). **It is not shown that an asymptotically true anti-concentration statement must see the slabs.**
4. *(Proved.)* Fibre restriction: `dist(l_u^m|_A, Poly) >= eta^m dist(l_v^m|_{F_t}, Poly)`, height-preserving, with
   `eta = dist_inf(u, L)` and `dist_inf(v, L) = 1`. The grid-model Lemma 5.3 is also proved.
5. *(Conditional.)* ILO_m (`kappa_V` form) plus the separation lemma gives MIN^pair at `H = infinity`, provided every
   projection-cheap `L` has a good collision count. At exponent `m` the room suffices from `n ~ 5.3 * 10^5` for
   16-class spans and from `n ~ 2.2 * 10^6` for the known 32-class spans `L_{2,1}`. It would not give growth even if
   proved, because of the `H`-gap. ILO_m is refuted by nondegenerate real `TC_V` points of regime fake boxes.
6. No growth improvement and no uniform bound: `M <= (2/3 + eps) log R/loglog R` stands. MIN^pair is open.
