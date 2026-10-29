# Papers check (judge): seven full texts against Question B

Inputs: four reader reports (`read_mr.md`, `read_ms_lms.md`, `read_huang_toric.md`, `read_huang_local.md`) and the
adversarial check `check.md` (74 items), all in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/papers_check/`. I also re-extracted the
two statements that carry the most weight (erratum p. 1 and LMS p. 2) to stdout. The repository was not edited, and
no PDF was copied.

The papers, with their abbreviations:

| abbr. | paper (uploaded file) | pages |
|---|---|---|
| MR-I | McKinnon-Roth, *Seshadri constants, diophantine approximation, and Roth's theorem for arbitrary varieties*, arXiv 1306.2976v2 (Invent. Math. 200 (2015) 513-583, per LMS p. 12) | 55 |
| MR-L | McKinnon-Roth, *An analogue of Liouville's theorem and an application to cubic surfaces*, arXiv 1306.2977v1 (Eur. J. Math. 2 (2016) 929-959, per LMS p. 12) | 17 |
| ERR | McKinnon-Satriano, *Erratum for Approximating rational points on toric varieties* (no date or venue printed; file name says 03-29-23) | 2 |
| LMS | Lehmann-McKinnon-Satriano, *Approximating rational points on surfaces*, arXiv 2403.02480v1 | 12 |
| H-ANT | Z. Huang, *Rational approximations on toric varieties*, arXiv 1809.02001v6 (ANT 15 (2021) no. 2, 461-512, per LMS p. 12) | 64 |
| H-I | Z. Huang, *Approximation diophantienne et distribution locale sur une surface torique*, arXiv 1703.01772v2 (Acta Arith. 189 (2019) 1-94, per H-ANT p. 61) | 64 |
| H-II | Z. Huang, part II, arXiv 1805.03920v2 (Bull. SMF 148(2) (2020) 189-235, per H-ANT p. 61) | 48 |

Page numbers are PDF page indices. Where a page number is printed, it equals the index. Labels:
* **[V]**: verified verbatim by a reader and by the checker.
* **[inf]**: an inference, by a reader, the checker or me. Each one is marked.

Throughout, `X = Lbar_5^tw` is smooth, with `-K_X` nef and big, and `X'` is its Gorenstein Fano contraction.

## 0. Bottom line

**None of the seven papers settles Question B, narrows it, or supplies a usable unconditional method for it.** None
contains a theorem stronger than tau/Liouville at a rational point on a Zariski-dense set.

Every unconditional lower bound in them is of one of three kinds:
* Liouville, `alpha_x >= eps_x` (MR-I p. 24);
* Liouville with the asymptotic base locus, `alpha^ess >= gamma` for `gamma <= gamma_eff = tau` (MR-L Thm 3.3);
* Roth or Schmidt type, `alpha^ess >= beta < tau` (MR-I Thm 5.1) and `alpha >= eps/2` (MR-I Thm 6.3).

Huang's lower bounds (H-ANT, H-I, H-II) are verified to be explicit single-section or product-of-sections Liouville
inequalities in Cox coordinates. Every exact value of `alpha_ess` in these papers equals the section value `tau`, and in
all of them `tau >= dim` [inf]. These are Huang's `Y_4`, `Y_3`, `Y_4 x P^n` and his Picard-rank-2 toric varieties.

The only statements reaching `alpha^ess_e(-K) >= 4` are conditional on Vojta's Main Conjecture:
* ERR Prop 0.1: smooth projective X, -K big;
* LMS Prop 1.4: smooth projective X.

Both apply as printed to the smooth X at e. Neither assumes splitness, and both use a Vojta instance on X itself with
a normal-crossings divisor through e, not V_E(5) on Bl_e X.

The papers also make it plausible that the τ < dim regime is open even for surfaces [inf]. LMS make essential
boundedness a *hypothesis* (Thm 1.5) and prove McKinnon's conjecture for split surfaces only conditionally on Vojta
(Cor 1.6).

## 1. Verified key statements (with pages)

### 1.1 Normalisation and definitions (all [V])

* **MR-I Defs 2.7-2.9, p. 13.** "We fix a single place v, archimedean or nonarchimedean."
  * For `x in X(kbar)`, `L` any line bundle, and `{x_i}` distinct points of `X(k)` with `d_v(x,x_i) -> 0`:
    `A({x_i},L) = {gamma in R | d_v(x,x_i)^gamma H_L(x_i) is bounded from above}`.
  * `alpha_x({x_i},L) = inf A`, and `alpha_{x,X}(L)` is the infimum over sequences.
  * "All height functions in this paper are multiplicative and relative to k."
  * This is the record's normalisation, `gamma = 4/(1+2 delta)` with `H = R''^2` and `dist ~ diam/R''`.
* **MR-I Prop 2.12, p. 15.** `alpha_x({x_i},L)_K = (d/m_v) alpha_x({x_i},L)_k`. For `Q(i)/Q` at infinity,
  `d/m_v = 2/2 = 1` [inf, checker item 34]. So lower bounds proved over `Q(i)` on the split form restrict to `T_5(Q)`
  with gamma unchanged.
* **LMS Def 1.1 (p. 2) and Defs 2.1-2.3, 2.7 (pp. 4-5).** The essential approximation constant `alpha_ess(P,A)` is
  "the infimum of all alpha({x_i},P,A) over all sequences of k-rational points {x_i} approximating P such that {x_i}
  is Zariski dense". Def 1.1 is stated for "any algebraic point P". Def 2.1 takes `P in X(k)`, a single place, and
  heights defined up to a bounded positive factor.
* **H-ANT Def 4.7, p. 29 (credited to Pagelot).** `alpha_ess(Q,Y) = sup_{V subset Y} alpha(Q,V)`, where V ranges over
  Zariski open dense subvarieties with `inf_{P in V(K)} d_nu(P,Q) = 0`. The same definition is in H-I Déf 2.3 (p. 7,
  "parties constructibles denses") and H-II Déf 2.2 (p. 10). The sup-over-U form equals the LMS
  inf-over-dense-sequences form [inf, checker item 55: diagonal argument over countably many proper closed sets].
* **LMS Def 1.3, p. 2.** "Let X be a smooth algebraic variety defined over a number field k, and let P in X be any
  k-rational point. … Then P is *essentially bounded* on X if for every ample divisor A and every rational eps > 0, we
  have alpha_ess(P, -K_X + eps A) >= dim X."
  * For X with -K big this is equivalent to `alpha_ess(P,-K_X) >= dim X` [inf, checker item 21: write
    `(1+eps/c)(-K) = (-K+eps A) + (eps/c)(-K - cA)` and use LMS Cor 2.8 concavity].
  * "Canonically bounded" is **not defined** in ERR. Its definition is in MS21, which is unread.

### 1.2 Unconditional lower bounds (all [V] unless marked)

* **Liouville, MR-I p. 24.** For `x in X(k)` and `0 < gamma < eps_x(L)`, `L - gamma E` is ample and has the Northcott
  property. "This proves the inequality alpha_x(L) >= eps_x(L), a stronger inequality than alpha_x(L) >= ½ eps_x(L)".
  For x of degree [K:k] the bound is `alpha_x >= eps_x/[K:k]`. At e this gives 2, which is sharp on the diagonal
  conics.
* **Liouville with base locus, MR-L Thm 3.3 (p. 7).** Hypotheses:
  * `x in X(kbar)` any point, `d = [K:k]`, `L` an **ample** line bundle;
  * `gamma > 0` rational with `L_gamma = pi^*L - gamma E` in the effective cone;
  * `B'` the asymptotic base locus of `L_gamma`, and `B = pi(B')`.

  Conclusions:
  * (a) If infinitely many `x_i` lie outside B, then `alpha({x_i},L) >= gamma/d`.
  * (b) If `alpha_x(L) < gamma/d`, then x lies in B and `alpha_x(L) = alpha_x(L|_B)`.
  * (c) A converse statement.
  * *Remark:* B may be replaced by the Zariski closure of `B(k)`.
  * *Cor 3.4:* `alpha_x(L) >= eps_x(L)/d`.

  At e it must be applied on X' with `-K_{X'}`, because L must be ample [checker item 12]. The values are unchanged.
  It gives `alpha^ess_e >= tau_e(-K) = 10/3` (Cilleruelo-Córdoba), and `15/4` on balanced points. **This is the
  ceiling of the mechanism: gamma <= gamma_eff = tau.**
* **Non-split quadric, MR-L p. 8.** X is a smooth quadric with no k-lines, and `pi^*L - 2E` is effective with base
  locus the two (non-rational) lines through x. So `B(k) = {x}` and `alpha_x(L) = 2` for all k-points. This is the only
  twist-sensitive device in these papers. It shrinks the exceptional set and never raises the exponent above
  `gamma_eff`.
* **Roth/Schmidt type, MR-I.**
  * *Def 3.1 (p. 20) and p. 25.* `eps_x(L) = sup{gamma : pi^*L - gamma E nef}`.
    `gamma_eff = sup{gamma : L_gamma effective}`, which is the record's `tau_x`. For **ample** L,
    `beta_x(L) = int_0^{gamma_eff} Vol(L_gamma)/Vol(L) d gamma`. Since `f <= 1`, f is continuous and
    `f(gamma_eff) = 0`, **beta < gamma_eff strictly**.
  * *Thm 5.1 (pp. 28-29).* X irreducible of dimension n over k, `x_v in X(kbar)`, L ample. If
    `sum_v beta_{x_v}(L) R_v > 1`, "there is a proper subvariety Z of X such that for all infinite sequences {x_i} of
    distinct points of X(k)\Z(k), there is at least one v in S so that alpha_v({x_i},L) >= 1/R_v". The proof uses
    Faltings-Wüstholz. This is **essential-type** and gives `alpha^ess_x >= beta_x(L)` [inf].
  * *Thms 6.1-6.2 (p. 33).* Dichotomies about the infimum `alpha_x`. The essential form of 6.2 is the "alternate
    statement" on p. 4: there is a proper `Z` such that only finitely many `y in X(k)\Z(k)` satisfy
    `d_v(x,y) < H_L(y)^{-((n+1)/(n eps_x(L)) + delta)}`.
  * *Thm 6.3 (pp. 4, 33).* For all X (possibly reducible), all `x in X(kbar)` and all ample L,
    `alpha_x(L) >= eps_x(L)/2`. Equality holds iff both constants are computed on a k-rational curve C with (1) C
    unibranch at x, (2) `kappa(x) != k`, (3) `kappa(x) subset k_v`, and (4) `eps_{x,C} = eps_{x,X}`. At a rational x the
    inequality is strict. Cor 6.4 has exponent `2/eps_x(L) + delta`.
  * *At e, on X'* [inf, readers 1 and checker item 4]: `beta_e(-K) in [2.31, 2.61]`. The lower bound is
    `(4/5)·70^{1/4}`, from Cor 4.2 (p. 27). The upper bound comes from tangent-line concavity of `Vol^{1/4}`, using
    `gamma_eff = 10/3` from round 23. **So the Roth-type essential bound at e is below Liouville's 10/3.**
  * *Thm 8.7 and Cor 8.9 (pp. 47-48).* The beta-hat / eps-hat versions over étale covers give nothing here. X and X'
    are simply connected, because the rays generate N [inf, checker item 9].
  * *Thm 10.1 (pp. 53-54).* Vojta for X with `-K_X` ample and D a finite set of points, assuming
    `eps_x(-K) > (n+1)/n`. At e this is exponent 1. **It is not V_E(5)** (pre-emptive caution: the record does not cite
    it).
* **Huang's lower bounds (section-type, verified by reading the proofs).**
  * *H-ANT, setting.* Thm 1.2/6.1 (pp. 5, 34): X is smooth projective toric of dimension at least 2, **split** over K,
    and satisfies Hypothesis (*), "The cone of pseudo-effective divisors Eff(X) is simplicial" (p. 5). The point is
    `Q_0 = (1,...,1) in T(K)`.
  * *H-ANT Thm 6.1 conclusions.* (1) If L is nef, then for every place nu, `alpha_{L,nu}(Q_0,X) = beta`, where beta is
    the minimal L-degree of a centred primitive relation. (2) If L is ample and `X != P^n`, alpha is properly achieved
    on a finite union of fibres `P^{N_i} x {1}`. The proof (p. 43) gives `alpha(Q_0, X\Y) >= beta + delta`, so
    `alpha_ess > alpha`.
  * *H-ANT Thm 7.5 (pp. 54-55).* For Picard number 2 and nef L, `alpha_ess` is computed exactly and equals the minimal
    L-degree of very free curves.
  * *H-ANT mechanism (pp. 40-48).* The algebraic integer `Y_{i0}(P_0) - X_{i0}(P_0)` is non-zero; the product formula
    bounds it from below, and Salberger's monomial sections compare it with the height. These are toric character
    differences of order 1 at `Q_0`, or products of two of them. Equations (64)-(69) and (92)-(95) were not checked.
  * *H-I.* Prop 5.1 (p. 37): `H(P) d(P)^2 >= 1` for `P != Q` on `Y_4`, a single-section Liouville bound. Thm 5.4
    (p. 38): `alpha_ess(Q) = alpha(Q,Y_4) = 2`, with no locally accumulating subvariety. Prop 7.1 (pp. 50-51):
    `alpha_ess = n+3` on `Y_4 x P^n`, a product of Liouville bounds.
  * *H-II.* Thm 1.4 (p. 4): on `Y_3`, `alpha = 2` and `alpha_ess(Q) = 5/2`. Prop 4.2 (pp. 16-17) proves
    `H(P) d(P)^{5/2} >= D` by a product of non-zero integers. This is Liouville for
    `F = (x-y)^2(t-s)^2(xs-yt)yt in H^0(2(-K_{Y_3}))` with `ord_Q F = 5` [inf, hand-checked by reader 4 and the
    checker], so `alpha_ess = tau = 5/2`.
* **Why Huang's theorems do not apply to X even over Q(i)** [inf, checker item 30].
  * (*) fails for `Lbar_3`, `Lbar_4` and `Lbar_5`, since no maximal cone satisfies H-ANT Lemma 6.3 (**). In the basis
    `f_i = 1_{A_i}`, the ray `1_{A_2\A_1} = f_2 - f_1`.
  * Thm 7.5 needs Picard number 2, but X has geometric rank 26.
  * Thm 6.1(2) needs L ample, but `-K_X` is not ample.
  * Huang works with smooth X, so the singular X' is not covered.
  * Run over `Q(i)`, the mechanism becomes the record's pair bound `|t_ab| >= 1` and its products, which is class L:
    10/3 for arbitrary profiles and 15/4 for balanced ones.
* **Curve formula, MR-I Thm 2.16 (p. 19; primary).** For C a singular k-rational curve, L ample on C, and normalisation
  `phi: P^1 -> C`: `alpha_{x,C}(L) = min_{q in phi^{-1}(x)} d/(r_q m_q)`. Here `r_q` is 0 if `kappa(q)` is not in
  `k_v`, 1 if `kappa(q) = k`, and 2 otherwise. LMS Thm 3.3 (p. 6) and H-I Rem 2.5 (p. 8) quote it with weaker
  hypotheses. This formula gives only **upper** bounds on alpha. It reproduces the record's values:
  * `2(max v - min v)` on `C_v`;
  * `15/(2·2)` for split Thm 2;
  * `30/(2·4) = 15/4` for the round-24 class, and `15/2` if its branches are rational.

  The quadratic field in Huang's Pell equations is the field of these branch slopes. It is a split-surface phenomenon,
  not anisotropy [inf].

### 1.3 Conditional statements (Vojta), all [V]

* **ERR p. 1.** "We correct an error in [MS21, Proposition 1.3]: Vojta's Conjecture does not imply canonical
  boundedness in general. However, the implication does hold for toric varieties. A corrected statement and proof of
  Proposition 1.3 follows."
  * **Proposition 0.1.** "Let X be a smooth projective variety over a number field k, and let P in X(k) be a k-rational
    point. If the anticanonical divisor class −K_X is big, then Vojta's Main Conjecture implies that X is canonically
    bounded at P."
  * *Proof.* `S = {v}`, and D is "the union of any n normal crossings divisors that intersect properly and
    transversely at P". Off a proper closed Z it shows `m_S(D,Q) >= -n log dist_v(P,Q) + C`. The blow-up is used
    only to compare proximity functions. The conclusion is `dist^{n-eps'} H_{-K_X} >= C_4`, and the final line prints
    `alpha_P(K_X)` (sign typo).
  * Despite the toric sentence, the proposition is stated for **every** smooth projective X with -K big. No
    splitness is assumed.
* **LMS Prop 1.4 (p. 2; proof pp. 5-6).** "Let X be a smooth, projective variety defined over a number field k, and
  let P in X(k). Then Vojta's Conjecture for X implies that P is essentially bounded." D is a normal-crossings divisor
  over k of multiplicity dim X at P, for example dim X transverse smooth divisors through P in `|mA|`. No bigness and
  no splitness are needed.
* **LMS Thm 1.5 / Cor 1.6 (p. 3).**
  * Thm 1.5: X a **split** (trivial Galois action on Pic, Def 4.1, p. 10), smooth, projective surface; A ample;
    `P in X(k)`. "If P is essentially bounded on X, then Conjecture 1.2 is true for P on X with respect to A."
  * Cor 1.6: the same, assuming Vojta for X.
  * Prop 3.2, Lemma 3.6, Thm 3.9 and Lemma 3.10 all take essential boundedness as a hypothesis.
  * The only unconditional `alpha_ess` statement in LMS is degenerate: it is infinite when rational points are not
    Zariski dense (p. 11).
* **Consequence for Q-B** [inf].
  * ERR Prop 0.1 and LMS Prop 1.4 give: Vojta for X implies `alpha^ess_e(-K_X) >= 4`, hence `A_5(delta)` for every
    `delta > 0`, hence the uniform bound. This uses round22 README §2 Step 1 and gauss Thm 3.4.
  * The transfer to X' is immediate, because phi is an isomorphism near e and `-K_X = phi^*(-K_{X'})`.
  * The Vojta instance used is (X, D through e, {inf}), not `(Bl_e X, E, {inf})`. The conclusion is the same.

## 2. Corrections to the record

Each item gives the file, the sentence (quoted or located), and the corrected text. Two attributions in the reader
reports must **not** be adopted:
* R2-8. LMS p. 2 credits [MR15] with cubic surfaces. That is a mis-citation in LMS: the cubic-surface theorem is MR-L
  Thm 4.5, i.e. MR16, so the record is right.
* R2-4, the "sup_U vs inf" discrepancy. The two definitions are equal.

### 2.1 `docs/claude_review/round22/lit.md`

1. **§0 item 3**, "they give no lower bound on `alpha_ess` beyond `alpha_e`. All are split or homogeneous; the
   conditional ones take their quantitative input from Vojta on `Bl_P X`, i.e. from `V_E(5)` itself."
   **Replace with:** "They concern `alpha_e`. Exceptions: Huang ANT Thm 6.1(2) (proof p. 43) and Thm 7.5 (pp. 54-55),
   and Huang II Thm 1.4 (`alpha_ess = 5/2 > alpha = 2` on `Y_3`). These give `alpha_ess > alpha`, by section-type
   Liouville bounds equal to `tau`. The unconditional cases are split or homogeneous. The Vojta-conditional
   boundedness statements (McKinnon-Satriano erratum Prop 0.1; LMS Prop 1.4) are **not** split. They use Vojta on X
   itself, with S = {v} and a normal-crossings divisor of multiplicity dim X at P, not Vojta on `Bl_P X` with D = E.
   Their conclusion `alpha^ess >= dim X` is the same as `V_E(5)`'s."
2. **§0 item 4**, "Roth-type theorems … are weaker than Liouville at a rational target."
   **Keep, and add:** "(confirmed: MR-I p. 24, 'a stronger inequality than alpha_x(L) >= ½eps_x(L)'; MR-L p. 1 for
   P^1). MR-I Thm 5.1 (p. 29) is essential-type, `alpha^ess >= beta`, but `beta < gamma_eff = tau` strictly (p. 25), and
   at e on X' `beta_e(-K) in [2.31, 2.61]` [inference]."
3. **§0 item 5 and §5 H2 (lines 270-273)**, "Acta Arith 2019" for both parts, and "studies a toric surface `Y_4`, using
   Pell-Fermat equations".
   **Replace with:** "Part I (arXiv 1703.01772): Acta Arith. 189 (2019) 1-94 (H-ANT bibliography p. 61). It studies
   `Y_4` (`P^1 x P^1` blown up at its 4 torus-fixed points, weak dP4), with `alpha = alpha_ess = 2` (Thm 5.4). Part II
   (arXiv 1805.03920): Bull. Soc. Math. France 148(2) (2020) 189-235 (H-ANT p. 61). It studies `Y_3` (3 of the 4 points
   blown up), with `alpha = 2` and `alpha_ess = 5/2` (Thm 1.4, p. 4)." Also refine "natural scale": Huang zooms at
   radius `eps B^{-1/r}` with `r >= alpha_ess`. For `r < alpha_ess` the only statement is H-I Prop 2.7, which is
   Liouville itself.
4. **§1.1, lines 51-59**, "[S: definitions as quoted in search snippets …]", and "`L` ample".
   **Replace with:** "[V: MR-I Defs 2.7-2.9, p. 13, for any line bundle L, `x in X(kbar)`, and one place v, archimedean
   or not; essential constant: H-ANT Def 4.7, p. 29 (Pagelot), verbatim as below, equal to LMS Def 2.7 (p. 5,
   infimum over Zariski-dense sequences)]."
5. **Lemma 1.1, lines 61-62**, "`L = -K_X`, which is ample (`X` is a smooth toric Fano variety …)". This was already
   killed in round 22 (F1). Note that it matters for citations too. MR-L Thm 3.3 and MR-I's beta require **ample** L,
   so they are applied on X' with `-K_{X'}`. Sections of `m(pi^*L - gamma E)` agree, so the values are unchanged.
6. **§1.3, lines 113-118**, "McKinnon-Roth 2016's Liouville analogue with stable base locus".
   **Replace with:** "MR-L Thm 3.3 (p. 7): L ample, gamma rational with `pi^*L - gamma E` effective, and B the image
   of the **asymptotic** base locus. By the Remark, B may be replaced by the Zariski closure of `B(k)`. The exponent is
   `gamma/[K:k]`." **Lines 156-161**, "`tau_e(-K) ∈ [10/3, 15/4]`, the exact value being open": superseded by round 23
   width Thm 2.5, `tau_e = 10/3`.
7. **§2, lines 172-176**, "Conjecture 2.7 there [S]".
   **Add:** "quoted verbatim as LMS Conj 1.2 (p. 2) and H-ANT Conj 1.1 (p. 3). The McKinnon-Roth form is MR-L Conj 4.2
   (p. 9) with Def 4.1 (p. 8): x a k-rational point, L ample, and a k-rational curve through x."
8. **§2, line 179**, "cubic surfaces, points off lines (McKinnon-Roth 2016)".
   **Replace with:** "split smooth cubic surfaces (MR-L Thm 4.5, pp. 9-10: X isomorphic over k to the blowup of P^2
   at six k-rational points). It covers **all** k-points for -K: `eps in {1, 3/2}`, `alpha in {1, 3/2, 2}`. Thms 4.6-4.7
   cover nef D at points not on a (-1)-curve. The mechanism is MR-L Thm 3.3(b) with gamma = 2 plus Thm 2.8. The
   non-split smooth quadric with no k-lines is also in MR-L (p. 8)." Keep the MR16 attribution. LMS p. 2's "[MR15]"
   is a mis-citation.
9. **§2, lines 180-181**, "San19 … author/journal not verified".
   **Replace with:** "probably [Ca19] D. Castañeda, *Rational approximations on smooth rational surfaces*, PhD thesis,
   Waterloo 2019 (LMS p. 11), credited with 'all split smooth rational surfaces of Picard rank at most five' (LMS
   p. 2). The identification is an inference."
10. **§2, lines 182-183**, "(Z. Huang, ANT 15 (2021) 461-512, §5)".
    **Replace with:** "(ANT 15 (2021) no. 2, 461-512, per LMS p. 12; Thms 1.2/6.1, H-ANT pp. 5, 34). Hypothesis (*):
    'The cone of pseudo-effective divisors Eff(X) is simplicial'. X is smooth projective and split. (1) needs only L
    nef. (2) needs L ample and `X != P^n`."
11. **§2, lines 188-192**, McKinnon-Satriano "(… with a 2023 erratum)", "[S]".
    **Replace with:** "LMS p. 2 summarises MS21 as 'Vojta's Conjecture implies Conjecture 1.2 for all Q-factorial
    split toric varieties' (no 'terminal'; unresolved without MS21). The erratum (undated; '2023' comes from the file
    name; LMS cite it as '[MS21a] … Trans. Amer. Math. Soc., to appear') corrects MS21 Prop 1.3 to Prop 0.1, quoted in
    §1.3 above. 'They also show that Vojta implies every point is essentially bounded' is now [V]: LMS Prop 1.4, p. 2."
12. **§2 verdict item 3, lines 201-205**, "comes **from Vojta's conjecture on `Bl_P X`** with `D = E`; that is exactly
    the record's `V_E(5)`".
    **Replace with:** "comes from Vojta's conjecture on X with S = {v} and D a normal-crossings divisor of multiplicity
    dim X at P (ERR p. 1; LMS pp. 5-6). This is a different instance from `V_E(5)`, with the same conclusion
    `alpha^ess >= dim X`." In the same item, "Huang toric … torsor-integrality" is now **verified** section-type
    (H-ANT pp. 40-48).
13. **§2 verdict item 4, lines 206-209**, "All proved cases are split". **Add:** "MR-L p. 8 treats a non-split
    (homogeneous) quadric, by the B(k) remark of Thm 3.3. Round22 README's 'split or homogeneous' already covers it."
14. **§3, lines 216-217**, "Full texts not read this round".
    **Replace with:** "Full texts read (papers_check, round 27): MR-I = 1306.2976v2, MR-L = 1306.2977v1 (v1, dated
    'Revised on January 18, 2013'; its numbering may differ from the EJM version)."
15. **§3 Statements, lines 220-224**: "a secondary source quotes '`alpha_x(L) >= epsilon_x(L)/2`' (Thm 6.3) [S]; on
    `P^n` the method gives the constant `n/(n+1)` [S]; … No exceptional set appears in the inequality itself".
    **Replace with:** "MR-I Thm 6.3 (pp. 4, 33) [V]. For all X, all `x in X(kbar)` and all ample L,
    `alpha_x(L) >= eps_x(L)/2`. Equality holds iff both constants are computed on a k-rational curve unibranch at x
    with `kappa(x) != k`, `kappa(x) subset k_v` and `eps_{x,C} = eps_{x,X}`, so the inequality is strict at rational x.
    `n/(n+1)` is MR-I Thm 6.2 (Schmidt type) for **any** irreducible n-dimensional X, not only P^n. Thm 6.3 and Cor 6.4
    have no exceptional set. Thm 5.1 (p. 29), the alternate Thm 6.2 (p. 4) and Thm 7.4 **do** have a proper
    exceptional Z, and give `alpha^ess >= beta_x(L)`."
16. **§3 Statements, lines 225-227 (EJM)**: see item 8.
17. **§3 Hypotheses, lines 229-232**: "a single place (archimedean distance)", "Exceptional sets: none in the
    Roth-type bound", "excludes the stable base locus".
    **Replace with:** "any single place, archimedean or not, with multi-place versions (MR-I §§5, 7). Exceptional
    sets: none in Thm 6.3 / Cor 6.4, and a proper Z in Thms 5.1, 6.2-alt and 7.4. The base-locus extension excludes
    `pi(asymptotic base locus of L_gamma)`, or only the Zariski closure of its k-points."
18. **§3 verdict item 3, lines 243-244**, "`beta_x(L) = ∫_0^∞ vol(L - tE)dt / vol(L)` … satisfy `beta <= tau`".
    **Replace with:** "`beta_x(L) = ∫_0^{gamma_eff} Vol(pi^*L - gamma E)/Vol(L) d gamma` for ample L (MR-I p. 25), with
    `beta < gamma_eff = tau` strictly. At e on X', `beta_e in [2.31, 2.61]` [inference]."
19. **§4 table**:
    * *MR16 row:* change "x off lines", "none (computes the infimum)" as in item 8.
    * *McKinnon-Satriano row:* add the erratum's Prop 0.1, which has no splitness and requires smooth X with -K big.
    * *LMS row:* replace "full conjecture for split surfaces | Vojta; split surfaces" with "Thm 1.5: split smooth
      projective surface + P essentially bounded ⇒ Conj 1.2; Cor 1.6: Vojta ⇒ the same; Prop 1.4: Vojta ⇒ every
      k-point of every smooth projective X is essentially bounded (no splitness)".
20. **§5 H3, lines 274-279**: "(simplicial pseudo-effective cone, per a citation in 2308.11847 [S])" is now [V]
    (H-ANT p. 5). Add that Thm 6.1 holds for every `Q in T(K)` (Remark 4.2(2)), every number field K and every place
    nu.
21. **§5 Hypotheses, lines 288-291**: "anticanonical (or ample toric) Weil heights". **Add:** "Huang allows L big
    (H-I p. 7; H-II p. 10) or nef (H-ANT Thms 6.1(1), 7.5)."
22. **§5 verdict 1, lines 294-296**, "*Split only.*" **Replace with:** "Not decisive. Lower bounds over `Q(i)`
    restrict to `T_5(Q)` with gamma unchanged (MR-I Prop 2.12). The decisive failures are:
    * Hypothesis (*) fails for `Lbar_3`, `Lbar_4` and `Lbar_5` (Lemma 6.3 (**) holds at no maximal cone);
    * Thm 7.5 needs Picard number 2 (X has 26);
    * Thm 6.1(2) needs L ample.

    And the mechanism is section-type (item 23)."
23. **§5 verdict 2, lines 297-302**, "[Assessment, not a theorem: I have not seen H3's proofs.]" **Replace with:**
    "[V: H-ANT pp. 40-48; H-I Prop 5.1 p. 37, Prop 7.1 pp. 50-51; H-II Props 4.1-4.2 pp. 15-17]: every lower bound is
    the non-vanishing of one algebraic integer (a toric character difference, order 1 at `Q_0`) or of a product of two
    such integers, compared with monomial sections."
24. **§5 verdict 4, lines 310-318, and §8, lines 401-403** ("In every computed example …"). **Add:** "`Y_4`
    (`alpha_ess = tau = 2 = dim`), `Y_3` (`alpha_ess = tau = 5/2 > 2`) and `Y_4 x P^n` (`n+3`) [inference: tau from
    explicit sections]. All have `tau >= dim`. 'No paper treats a rational point at which `alpha_ess` is shown to
    exceed `tau_x(L)`' is now confirmed for the seven full texts."
25. **§9**: items 1, 2 and 3 are done. Item 3's question ("Pell-Fermat … closest analogue of our anisotropic
    situation") is **answered negatively**: the quadratic field is that of nodal-curve branch slopes on a split
    surface (MR-I Thm 2.16 with `r = 2`). Item 6 is partly done: LMS and the erratum are read, MS21 and Monahan-Satriano
    are not.

### 2.2 `docs/claude_review/round22/README.md`

1. **§0 item 1, lines 29-31**, "In the words of McKinnon-Satriano (TAMS 2021), the strongest form is '`X` is
   canonically bounded at `e`'".
   **Replace with:** "The strongest form is that e is *essentially bounded* on the smooth X (Lehmann-McKinnon-Satriano
   Def 1.3, p. 2: `alpha_ess(e,-K_X+eps A) >= dim X` for every ample A and rational eps > 0). Since -K_X is big, this
   is equivalent to `alpha^ess_e(-K) >= 4` [inference]. McKinnon-Satriano's 'canonically bounded' (MS21; definition
   unread) is the notion in their erratum."
2. **§0 item 2, lines 35-38**, "Roth-type bounds satisfy `beta <= tau`; … the strong statements are split, or
   conditional on Vojta".
   **Replace with:** "Roth-type bounds satisfy `beta < tau` strictly (MR-I p. 25); MR-I Thm 5.1 makes `beta` an
   essential bound, `in [2.31, 2.61]` at e. The Vojta-conditional boundedness statements are not split."
3. **§1.1 item 1, lines 55-59**, "No full text could be read … Every statement about the *mechanisms* … is therefore
   **assessment**". **Replace with:** "Seven full texts were read in papers_check (MR-I, MR-L, ERR, LMS, H-ANT, H-I,
   H-II). The mechanisms of McKinnon-Roth and Huang are now verified, and they are section-type."
4. **§1.1 item 5, lines 77-79**, "concern `alpha_x`, not `alpha^ess`". **Add:** "except Huang ANT Thms 6.1(2), 7.5 and
   Huang II Thm 1.4, which give `alpha_ess > alpha` by section-type bounds (= `tau`)."
5. **§1.1, lines 104-105**, "two citation details are unconfirmed". **Replace with:** "Both are confirmed. ANT 15
   (2021) no. 2, 461-512 (LMS p. 12, secondary). 'alpha >= eps/2 (Thm 6.3)' is in MR-I p. 4 and p. 33 (primary), for
   all X, `x in X(kbar)` and ample L, with the equality clause of §1.2."
6. **§1.2 item 7, lines 153-154**, "Vojta for `(Bl_e X, E, {inf})` gives `alpha^ess >= dim X`. This is gauss Thm 4.2
   and, at snippet level, McKinnon-Satriano's corrected Prop 1.3."
   **Replace with:** "… This is gauss Thm 4.2. McKinnon-Satriano's corrected statement (erratum Prop 0.1, replacing
   MS21 Prop 1.3) and LMS Prop 1.4 reach the same conclusion from a **different** Vojta instance: X itself, S = {v},
   and D = n transverse normal-crossings divisors through P. They require X smooth projective, and the erratum also
   requires -K big. [V]"
7. **§2 table**:
   * *Row "McKinnon-Roth Roth-type (Invent. 2015)":* add "MR-I Thm 5.1 (essential): `alpha^ess >= beta_e(-K_{X'}) in
     [2.31, 2.61] < 10/3`".
   * *Row "McKinnon-Roth Liouville … (EJM 2016)":* add "MR-L Thm 3.3, p. 7; L ample, so on X'".
   * *Row "Huang torsor lower bounds (ANT 2021, Acta Arith 2019)":* replace "no (assessment) | section-type … |
     **assessment**, proofs unread" with "no | section-type, verified (H-ANT pp. 40-48; H-I Prop 5.1; H-II Prop 4.2); for
     X also (*) fails, Picard rank 26, -K not ample | full text". Name Part II as Bull. SMF 148 (2020).
   * *Row "McKinnon conjecture, Coba …":* "snippet-level" becomes "full text for MR-L, LMS and H-ANT".
   * *Row "McKinnon-Satriano (TAMS 2021), Monahan-Satriano, Lehmann-McKinnon-Satriano":* replace "conditional on
     Vojta; split | snippet-level" with "conditional on Vojta (erratum Prop 0.1: smooth X, -K big; LMS Prop 1.4: smooth
     X); **no splitness** in these implications. Splitness enters only LMS Thm 1.5/Cor 1.6 and, per LMS p. 2, MS21's
     main theorem | full text for the erratum and LMS; MS21 and Monahan-Satriano unread".
8. **§2 GAP 1, lines 210-212**, "Their exact definition … is known only from snippets. So is whether they require `X`
   smooth". **Replace with:** "LMS Def 2.7 is verified. The erratum requires X **smooth** projective, `P in X(k)` and
   -K big, so it is applied to `X = Lbar_5^tw`, not to the singular X'. The transfer is trivial near e. Only MS21's
   definition of 'canonically bounded' remains unread."
9. **§2 GAP 2, lines 213-215**: add "and by Vojta for (X, D, {inf}) with D normal crossings of multiplicity 4 at e
   (erratum p. 1; LMS Prop 1.4)".
10. **J1, lines 240-247**, "Both are snippet-level. Nothing in this README depends on the erratum's exact wording."
    **Replace with:** "Verbatim (erratum p. 1): the opening sentence, quoted in §1.3 of papers_check README, and
    Prop 0.1. The toric sentence is motivational. The proposition holds for every smooth projective X with -K big.
    Both killers' paraphrases are accurate fragments."
11. **J3, lines 273-275**, "(`Q`-factorial terminal **split**)": LMS p. 2 says "Q-factorial split toric varieties".
    Mark this as unresolved until MS21 is read.
12. **§4 row 21**, "Huang **assessment** (proofs unread)" becomes "Huang: verified section-type". **§6, lines
    396-399** (owner-side enabler): the seven uploads are done; list the remaining papers of §3 here.

### 2.3 `docs/claude_review/expert_note/expert_note.md` §3

1. **Approximation reading, bullet** "Roth, Schmidt and Ru–Vojta with `S={∞}` do not apply. The target is rational and
   `E` is geometrically irreducible, so `β(L,E)≤τ(L,E)`, which never goes beyond Liouville."
   **Replace with:** "Roth- and Schmidt-type theorems apply but are weaker. McKinnon–Roth (Invent. Math. 200 (2015))
   Thm 5.1 gives, off a proper closed set, `α^ess_x ≥ β_x(L)`, where
   `β_x(L)=∫_0^{τ}\mathrm{Vol}(π^*L−γE)/\mathrm{Vol}(L)\,dγ<τ_x(L)` (ibid. p. 25). Their Thm 6.3 gives `α_x ≥ ε_x/2`,
   strict at rational x. At e (on X'), `β_e(−K)∈[2.31,2.61]` (inference), below Liouville's 10/3. As McKinnon–Roth
   note (ibid. p. 24), at a k-rational point Liouville `α_x ≥ ε_x` is 'a stronger inequality than `α_x ≥ ½ε_x`'."
2. **Same paragraph, end**, "Round 22 checked this at snippet level (two independent searches plus an offline arXiv
   mirror) …". **Append:** "Seven full texts were then read (McKinnon–Roth 2015, 2016; the McKinnon–Satriano erratum;
   Lehmann–McKinnon–Satriano 2024; Huang ANT 2021, Acta Arith. 2019, Bull. SMF 2020). Every unconditional lower bound
   on `α^ess` in them is Liouville with a section, at most `τ_x`. Every exact value of `α^ess` they compute equals `τ`,
   with `τ ≥ \dim`."
3. **"In the approximation-constant literature's own terms"**, "McKinnon–Roth's base-locus Liouville at `e` is
   therefore exactly Cilleruelo–Córdoba." **Add:** "(McKinnon–Roth, Eur. J. Math. 2 (2016), Thm 3.3. It needs L ample,
   so it is applied to `−K_{X'}` on X'. The values are unchanged.)"
4. **Question B, literature form**, the whole box. **Replace** "**canonically bounded** at the identity `e` in the
   sense of McKinnon–Satriano (TAMS 374 (2021)), i.e. does every Zariski-dense sequence of rational points tending to
   `e` have `α(−K)≥\dim X=4`?" **with:** "Is the identity e **essentially bounded** on the smooth weak Fano 4-fold
   `X=\overline L{}^{\rm tw}_5` in the sense of Lehmann–McKinnon–Satriano (arXiv 2403.02480, Def 1.3:
   `α_{\rm ess}(e,−K_X+εA)≥\dim X` for every ample A and rational ε>0)? Since −K_X is big, this is equivalent to
   `α^{\rm ess}_e(−K)≥4` (inference). The same constant is obtained on X', since φ is an isomorphism near e."
5. **Same box, last sentences**, "In McKinnon–Satriano canonical boundedness is a **hypothesis**; by their 2023 erratum
   Vojta implies it when −K is big. The record found no unconditional proof for any variety (snippet level; full texts
   unread)." **Replace with:** "Vojta's Main Conjecture implies it. McKinnon–Satriano's erratum (correcting their
   TAMS 374 (2021) Prop 1.3) has Prop 0.1: 'Let X be a smooth projective variety over a number field k, and let
   P ∈ X(k) … If −K_X is big, then Vojta's Main Conjecture implies that X is canonically bounded at P.'
   Lehmann–McKinnon–Satriano Prop 1.4 gives the same for essential boundedness on any smooth projective X. Neither
   assumes splitness. Both use Vojta on X with a normal-crossings divisor of multiplicity dim X through P, a different
   instance from `V_E(5)` with the same conclusion. In Lehmann–McKinnon–Satriano Thm 1.5 essential boundedness is a
   hypothesis. No unconditional proof was found for any variety where `τ_P(−K) < \dim X`. Full texts were checked for
   seven papers; McKinnon–Satriano (TAMS) and Monahan–Satriano remain unread."
6. **"Vojta predicts it" bullet.** Optional addition: "The literature's instance (Vojta on X, D ⊃ e normal crossings,
   multiplicity 4) gives the same `α^{\rm ess}≥4`."

### 2.4 `docs/claude_review/expert_note/cg_questions.tex`

1. **Lines 639-641**, `\item Roth, Schmidt and Ru--Vojta with \(S=\{\infty\}\) do not apply. … never goes beyond
   Liouville.` Replace as in 2.3 item 1 (LaTeX). Use `\(\alpha^{\rm ess}_x\ge\beta_x(L)\)`,
   `\(\beta_x(L)<\tau_x(L)\)` and `\(\beta_e(-K)\in[2.31,2.61]\)` (inference).
2. **Lines 645-649**, after "…was found for any variety (\rec{round22/README.md} \S2, J3)." Append the full-text
   sentence of 2.3 item 2.
3. **Lines 689-690**, "McKinnon--Roth's base-locus Liouville at \(e\) is therefore exactly Cilleruelo--C\'ordoba."
   Add the ampleness parenthesis of 2.3 item 3.
4. **Lines 697-703 (`named` box)**, "\emph{canonically bounded} at the identity \(e\) in the sense of
   McKinnon--Satriano (TAMS 374 (2021)), i.e.\ does every Zariski-dense sequence … \(\alpha(-K)\ge\dim X=4\)?" Replace
   as in 2.3 item 4.
5. **Lines 707-711**, "In McKinnon--Satriano canonical boundedness is a \emph{hypothesis}; by their 2023 erratum
   Vojta implies it when \(-K\) is big. The record found no unconditional proof for any variety (snippet level; full
   texts unread)." Replace as in 2.3 item 5.
6. **Line 1031**, `\item Roth and Schmidt need an irrational target, and \(e\) is rational.` Replace with
   `\item Roth- and Schmidt-type bounds (McKinnon--Roth) give \(\alpha^{\rm ess}\ge\beta<\tau\) at a rational target:
   weaker than Liouville.`
7. **Bibliography (lines 1125-1146).** Add entries for the cited papers, with the data as printed in other uploaded
   files:
   * McKinnon–Roth, Invent. Math. 200 (2015) 513–583;
   * McKinnon–Roth, Eur. J. Math. 2 (2016) no. 4, 929–959;
   * McKinnon–Satriano, Trans. AMS 374 (2021) no. 5, 3557–3577, and the erratum;
   * Lehmann–McKinnon–Satriano, arXiv 2403.02480;
   * Huang, ANT 15 (2021) no. 2, 461–512.

*Outside the four files (FYI).* round23/README.md line 93, "Vojta, via McKinnon-Satriano's erratum for toric
varieties, predicts it true": the erratum's Prop 0.1 is for all smooth projective X with -K big. It applies to the
smooth split `LM_5`.

## 3. Papers still needed, and what to look for

Bibliographic data in this section that is not printed in an uploaded PDF is from recall, not verified. This
applies to the Philippon, Broustet and Ru-Vojta data.

1. **McKinnon-Satriano, TAMS 374 (2021) 3557-3577, arXiv 2004.05212.** Highest priority. Look for:
   * (a) the definition of "canonically bounded": whether it is `alpha_ess(P,-K) >= dim`, or `-K + eps A`, or a sup
     over open sets, and whether X must be smooth, Q-factorial or terminal;
   * (b) the original Prop 1.3 and the counterexample behind "does not imply … in general";
   * (c) the exact main theorems: smooth split toric surfaces? "Q-factorial (terminal?) split toric", as LMS p. 2
     summarises?
   * (d) **any unconditional case of canonical boundedness**, in particular at a point with `tau_P(-K) < dim X`;
   * (e) how they bound Zariski-dense sequences on toric varieties: whether any argument other than Vojta or sections
     appears.
2. **Monahan-Satriano, arXiv 2308.11847.** Look for:
   * the exact hypotheses (split horospherical, smooth projective; Vojta);
   * whether any statement is unconditional;
   * whether non-split forms or anisotropic tori are mentioned;
   * the exact wording of the citation of Huang's hypothesis (*);
   * any lower bound on `alpha_ess` that is not Vojta or a section bound.
3. **Philippon 1986, "Lemmes de zéros dans les groupes algébriques commutatifs", Bull. SMF 114 (1986) 355-383.** This
   matters for round 26's route-5 zero estimate, not directly for the `alpha^ess` form. Look for:
   * the exact multiplicity estimate on `G_m^n` (and on general commutative groups over a number field), including the
     obstruction subgroup H, the degree/multiplicity inequality and the constants;
   * whether the statement is over an arbitrary field, so that it applies to `T_5` after base change to `Q(i)`;
   * whether multiplicities along a subspace of the Lie algebra are allowed.

   Then check whether the record's "actual point" zero estimate fits that form.
4. **McKinnon, JAG 16 (2007), arXiv math/0504303.** Look for Conj 2.7 verbatim, any example where `alpha_ess` is
   computed, and in particular any example with `tau_P(-K) < dim`.
5. **A del Pezzo surface of degree 2 as a 2-dimensional calibration** [inference, mine; not in any paper read].
   * *The computation.* At a general rational point x of a dP2, take `L = pi^*(-K)` on `Bl_x S`. The member
     `D in |-2K|` with a triple point at x (6 conditions on `h^0 = 7`) has `Dtilde = 2L - 3E` and `Dtilde^2 = -1`.
     Zariski decomposition gives `Vol(L - tE) = 2(3 - 2t)^2` on `[4/3, 3/2]`. So `eps_x(-K) = 4/3` (consistent with the
     known value for dP2 at a general point) and **`tau_x(-K) = 3/2 < 2 = dim`**. This assumes D is irreducible, which
     should be checked.
   * *Along D.* D is rational (`p_a = 3`, with a triple point). MR-I Thm 2.16 gives `alpha_D >= 4/2 = 2`.
   * *The window.* So `alpha^ess_x(-K)` lies in `(3/2, 2)` unless Vojta fails. This is exactly the phenomenon of
     Question B ("alpha^ess > tau at a rational point") on a surface. LMS's conditional treatment of split surfaces
     suggests it is open.
   * *Papers to look for.* Broustet (2006) on Seshadri constants of del Pezzo surfaces. Castañeda's 2019 thesis, to see
     whether it reaches degree <= 4 and with what mechanism. Any paper computing `alpha_ess` on dP2 or dP1.
6. **Ru-Vojta (Amer. J. Math. 2020), Heier-Levin, Evertse-Ferretti.** Check that for D a sum of 4 **Q-rational**
   divisors through e, with S = {inf}, their Schmidt-type bounds collapse to Liouville. With rational linear forms and
   only the archimedean place, the subspace input is trivial integrality [inference]. This closes the remaining "Vojta
   on X with D through e" lane that the erratum and LMS use.
7. Lower priority:
   * Huang IJNT 13 (2017) (split dP6, where (*) fails): the mechanism and value of `alpha_ess` at e;
   * Batyrev-Tschinkel IMRN 1995 (anisotropic tori): local Fourier analysis at infinity in shrinking balls;
   * Pagelot's 2008 thesis (origin of `alpha_ess`).

## 4. Single most promising next step for Question B

**Put the sharpened, literature-native question to McKinnon, Satriano (and Huang), with the MS21 definitions read
first.** The question:

> *Is essential boundedness (LMS Def 1.3) known unconditionally at **any** rational point P of **any** smooth
> projective variety with `tau_P(-K) < dim X`? The simplest case may be a general rational point of a degree-2 del
> Pezzo surface, where (by our computation) `tau = 3/2 < 2`. Our case is the identity of the smooth weak Fano 4-fold
> `Lbar_5^tw`, with `tau_e(-K) = 10/3 < 4`. There `alpha^ess_e(-K) > 15/4` on balanced sequences, or on all
> sequences, would already give the uniform Cilleruelo-Granville bound.*

Why this step:
* The papers settle that no known unconditional method passes `tau`.
* The authors themselves keep essential boundedness as a hypothesis (LMS Thm 1.5) or derive it from Vojta (ERR,
  LMS Prop 1.4).
* So the only remaining value is an expert pointer to an unpublished or unindexed mechanism, or confirmation that the
  `tau < dim` regime is open everywhere. Confirmation would let the record close Question B as "Vojta-hard, open even
  for surfaces".
* Reading MS21 (2004.05212) first is cheap and fixes the definitions the question must use.

Honest probabilities [judgement]:
* MS21 or Monahan-Satriano contains an unconditional result implying or narrowing Question B: **about 1%**.
  LMS (2024), by the same authors, list only conditional results.
* An expert identifies an existing unconditional theorem giving `alpha^ess > tau` at a rational point for some
  variety: **about 3-5%**.
* Such a theorem transfers to X' at e and gives `> 15/4`: **about 1%**.
* Question B proved by anyone within a year: **< 1%**.
* The dP2 calibration confirmed open by experts, and accepted as the right 2-dimensional analogue: **about 60%**. This
  has map value only.

## Files

* This judge report: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/papers_check/README.md`
* Reader reports, in the same directory:
  * `read_mr.md`
  * `read_ms_lms.md`
  * `read_huang_toric.md`
  * `read_huang_local.md`
* Checker report: `check.md`, in the same directory.
* Helper scripts and outputs, in the same directory:
  * `huang_hyp_check.py` -> `.out` (Hypothesis (**) count);
  * `huang_lines_check.py` -> `.out`.
