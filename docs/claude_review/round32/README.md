# Round 32: judge's report on mining the posted Lean "RH/GRH" development for the short-arc problem

Inputs: six analyst reports in this directory (`mine_architecture.md`, `mine_stepanov.md`, `mine_slopes.md`, `mine_cm.md`,
`mine_positivity.md`, `mine_lifts.md`; all six present and read in full), the 4-page remark
(`../dl_remark/remark.pdf`), the Lean tree (`../dl_lean/x/lean/src`, called S2 below; read-only, never built or run),
and the record (`/home/user/Jarnik/docs/claude_review`: expert note sections 1 and 5, `round31/README_XI.md` section 5).
Labels: **[J]** = checked or derived by the judge in this round (scripts in `judge/`); **[A:name]** = an analyst's
claim that I accept after reading; **[Inferred]** = interpretation, not checked against proofs.

Terms used throughout. *Cluster*: M lattice points on an arc of length at most C sqrt(R) of |z| = R, N = R^2,
divided by their Gaussian gcd (normalised norm N', W = log N'). *Rate*: proving M = o(log R / loglog R); the record has
M <= (2/3 + o(1)) log R/loglog R (Theorem 1.1, i.e. 3 M log M <= W + O(M)). *Fake*: an angle model keeping the
genuine primes and exponent profile with real angles; the record has fakes at W = O(M log M) satisfying pair, gcd,
cut, monomial-gap, rational-character and (within scope) point- and character-tested residue inputs, so no such
input gives the rate. *Cut identity*: the exact product formula for prod |z_i - z_j| (expert note section 1).
*Route 5*: auxiliary determinants in monomials of the cluster points with forced divisibility, blocked by the open
zero estimate ZE^all_H (expert note section 5 item 3).

## 0. Answer in five lines

1. No. Nothing in the method gives, or plausibly leads to, M = o(log R/loglog R); my estimate for the whole
   development as a source of the rate is below 1%, and below 0.5% for the general worst case.
2. Its engine is the same scheme as our Theorem 1.1, a product-formula ("divisibility <= size") inequality: p-adic
   factorial divisibility at inert primes against Stirling-size archimedean cost. The new part is amplifying it
   m ~ sqrt(x) times by a chain of theta lifts, and that amplification reduces an *error term* of an averaged prime sum.
3. Every concrete transfer the analysts tried either (a) holds for fakes, (b) is single-place Liouville (exactly
   critical: on one circle, multiplicity just raises the cut identity to a power), (c) is already in the record, or
   (d) needs route 5's zero estimate. The striking "a/p + a/p^3 + ..." exponent is, on the circle, the all-level inert
   collision count; I checked it exactly: it adds a constant 0.2287 * M^2/2, i.e. O(M) to W, and the constant 3 stays 3.
4. The one lesson worth keeping is negative and conceptual: S2 never proves a zero estimate; it makes nonvanishing free
   (positivity, exact period/norm identities) and puts all arithmetic into divisibility. For arcs, free nonvanishing is
   already known to be subcritical, and the supercritical version is exactly ZE^all_H.
5. Recommendation: no new experiment. At most, add three short lane-closing entries to the record (section 4) so these
   ideas are not retried.

## 1. What the Lean development contains

**Status.** S2 is a posted, uncompiled (by us) Lean 4 tree (toolchain `leanprover/lean4:v4.33.1`) of **175,187**
`.lean` files with **no comments or docstrings**; proofs often end in machine-generated `first | exact ...` cascades.
My text scan [J]: no occurrence of `sorry`, no `axiom` declaration (the single regex hit is a structure field named
`axioms`), no `native_decide`, no `implemented_by`/`extern`. This means only that there are no *textual* holes. It does
**not** mean that the tree compiles, that the kernel accepts it, that the definitions mean what their names suggest,
that the many `def X : Prop` "Statement" files are ever proved, or that the main theorem's hypotheses are not
smuggled in through structures. Concretely [J]: `TruncFinal`, `LatticeChainBound`, `GramZ` (SmithDeterminants) and
`LandauStep` (Formalism/Statement.lean) are `def ... : Prop` referenced nowhere outside their own definition files, so
they are stated but neither proved nor used under those names. `AskMainBalance`, by contrast, is proved
(`theorem colour_balance`, `ChainSlopeWindows/StronglyIrreducible.lean` line 695). Nothing below certifies any proof.

**Skeleton** [A:architecture, spot-checked].
* Top: `allChi : EveryCharacter` (L(s,chi) != 0 for Re s > 1/2) := `FinalTheorems.everyCharacter`, which is
  `everyCharacter_of_row` applied to a claimed theorem `galoisLocalGlobalStrong_recParameter`; it then goes through
  `every_character_of_literature_of_vWin` (three "open literature" inputs plus `VWin`, a quantitative weakening of
  Vandiver: total log-size of certain Sha-type groups <= b^(1/2)), then theta-integral local matching, and finally
  `zeroFreeAll_of_topG`.
* Analytic shell: Landau's lemma. A one-sided bound on a Cesaro-weighted Chebyshev error with nonnegative
  coefficients, G(x) <= x^(1/2+O(eps)) on windows, gives no zeros in Re s > 1/2. This part is classical.
* Quantitative heart: one averaged inequality `MeanBoundSecondForm`: `umean ASK <= umean Main + ErrX`. ASK is an
  inert-prime divisibility count sum_p 2 log p sum_i (a/p^i - w)_+; Main is the archimedean size (m copies of
  Stirling a log a - a, normalised by a Chowla-Selberg period). The exact identity `AskMainBalance` says
  ASK - MAIN is a character average of twisted Chebyshev errors, which is how divisibility controls primes.
* Amplification: m theta lifts U(N(0)) -> ... -> U(N(m)), each adding a copy of the gain to one Arakelov degree.
  Parameters: m in [x^(1/2-eps), ...] (`RestSize`, checked [J]), b = exp(x^eps), a = x b, the group H = Z/2^e with
  2^e ~ x^(2+8eps) (a layer of the cyclotomic Z_2-tower, so every prime <= x has a long Frobenius orbit), and Behrend
  sets to select the "top places". A weaker variant (`ZetaQuarterPowerCorollary`, checked [J]) uses m ~ x^(1/4) and
  concludes zeta(s) != 0 for Re s > 3/4.
* Lemma layer: Fock/Gram determinants (positive definite, factorial divisibility, Hadamard), Smith-form conversion of
  "y^m kills mod varpi^S" into det divisibility, Harder-Narasimhan slopes and the unfiltered slope inequality,
  Dwork-type jet bounds for theta functions in a height-2 formal group (exponent gam(l) = l/(l^2-1)), and finite
  theta-lift linear algebra (seesaw = finite Fubini, Rallis = T T* = scalar, denominators multiply along a chain).
* Not found by any analyst: the remark's "2/3" (only gamma < 2/3 occurs), a Stepanov module, the path-quadric toy
  proof, Weil positivity W(1_S) >= 0, Hurwitz/Kummer congruences, Siegel-Weil at degenerate targets.

**Spot-checks of analyst READ claims [J]** (opened the file, compared with the paraphrase).

| claim (analyst) | file | result |
|---|---|---|
| `rho l = 1/(l^2-1)`, `gam l = l * rho l` (cm) | `Formalism/FrobeniusBaseCase/Statement.lean` l.14-16 | matches |
| `TruncFinal`: L(y^m v) in varpi^S O^d => varpi^(dS) divides varpi^(dmt + m n_t S) det L (stepanov) | `Formalism/SmithDeterminants/StatementDVR.lean` l.108 | matches; it is a `def : Prop`, unused elsewhere |
| `LandauStep`, one-sided `-C x^sigma0 <= G chi x` (positivity) | `Formalism/Statement.lean` l.32 | matches; `def : Prop`, unused elsewhere |
| `ASK`, `MAIN`, `AskMainBalance` (architecture) | `Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean` l.84-101 | matches; proved as `colour_balance` |
| `ChainScalar`: T_j T_j* = c_j on S_j => composite = prod c_j (lifts) | `Formalism/ThetaLiftIndependence/Statement.lean` l.244 | matches |
| `RestSize`: x^(1/2-eps) <= m, k <= n^2 x^(1/2-eps) (architecture, lifts) | `Formalism/SizeWindowBrackets/SChain.lean` l.53 | matches |
| `SlopeInequality`: slope L <= muMax M + log C for injective C-Lipschitz phi (slopes) | `Formalism/LatticeSlopes/Statement.lean` l.75 | matches |
| `ZetaQuarterPowerCorollary`: m ~ x^(1/4), k ~ x^(1/2), l ~ x^(1/4) => no zeros for Re s > 3/4 (slopes) | `Formalism/ProductFormula/StatementOCCM.lean` l.137-151 | matches |
| `lpF` = l/(l^2-1) sum Ds - sup of the orbit LP (architecture) | `Formalism/LocalRingCardinality/DataParams₂.lean` l.51-59 | matches |
| `GramZ` (stepanov) | `Formalism/SmithDeterminants/StatementGram.lean` l.37 | matches; the paraphrase omits one extra conclusion (each alpha_a! divides s!) |
| `OrdLargeJ` (architecture) | listed under `BehrendSets/Statement.lean` | it is in `StatementE.lean` l.13; content matches |

No substantive mismatch. Two minor imprecisions (a file location and one omitted conjunct). The real caveat is the one
above: several quoted "lemmas" are unproved, unused `Prop` definitions, so describing them as S2's mechanism is an
inference about intent.

## 2. Idea-by-idea table

Tests: **F** = fake test (do fakes satisfy the input? if yes, no rate); **P** = characteristic-p test (valid verbatim
over F_p(t)? if yes, no uniform bound); **L** = single-place Liouville (exactly critical); **R** = already in the record.

| # | idea | what it is in S2 | analogue for arcs | F | P | L | R | verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | One-sided bound suffices (Landau) | `LandauStep`, `ZetaLandau` | our count is already one-sided; there is no analytic family whose continuation it controls | n/a | n/a | n/a | — | dies: nothing to transfer |
| 2 | Divisibility - size = character average (`AskMainBalance`) | proved identity | collisions mod inert q expanded in characters of T(F_q) | yes (P3-char capped at every modulus) | yes | — | round5/construct, README_IX | dies |
| 3 | m-fold chain amplification | m theta lifts, m ~ sqrt(x) | powers z^k; product clusters z_I; tower z -> z^k | yes (angles add) | yes | homogeneous: gain = cost | expert note 5.3 | dies: k-th power cluster has length kC R^(k-1/2) >> critical R^(k/2); product pairs carry only base pair data |
| 4 | Multiplicity determinant (Stepanov) | Smith form + lift chain | confluent Vandermonde on the 1-dim torus = V^(a^2) x factorials | yes (cut identity) | yes | yes | 5.3 (route 5 in higher dim) | dies |
| 5 | Free nonvanishing by positivity | `GramZ`: det of PD integral Gram > 0 | Fock Gram of evaluation functionals = sum over Schur functions times V; arithmetic = pair data | subcritical before fakes matter | yes | — | 5.3 ("free-nonvanishing grids are subcritical") | dies |
| 6 | Multiplicity as varpi-adic divisibility of a Frobenius-degenerate operator | `TruncFinal`, `LatticeChainBound` | z^((l+1)l^(k-1)) = N^(l^(k-1)) mod l^k at inert l; on a rank-M evaluation it returns the l-adic valuation of V | yes (point residues) | yes | — | Thm 1.1 collisions, Thm 1.2 | dies |
| 7 | Inert exponent gam(l) = l/(l^2-1) = 1/l + 1/l^3 + ... | `GaussJetBounds`, height-2 formal group | all-level inert collisions: sum v_l(z_i - z_j) >= sum_k (M^2/(2(l+1)l^(k-1)) - M/2)^+ | yes | yes | — | two_thirds/proof.md Lemma 3.3(I), line 656 | dies; [J] verified, see below |
| 8 | Averaging instead of pointwise control | `TheoremM`, averaged `SlopeIneq` over twist "words" | average over orders/degrees, twists mod tau, sign flips, lambda^(4k) | yes (grafts, P3-char, profile-only) | yes | — | 5.1, 5.2 | dies |
| 9 | Filtered Bost slope inequality | HN + `SlopeIneq` (unfiltered only) | jets on P^1 / the conic: K^2 x cut identity; Bombieri-Pila conic det = N^(d(d+1)/2) V / prod z^d | yes | yes | yes | round20/arakelov.md L3.4, B.2 | dies; [J] identity verified |
| 10 | Slope along t -> z^t instead of a zero estimate | Bost injectivity | Schneider-Lang/Baker on G_m^d | single forms faked | — | — | round26 ze.md; 5.3, 5.4 | dies |
| 11 | Projector cut / descent (`GoodBadCut`, `DescentSlope`) | Arakelov bookkeeping | exponent-row lattice (Z^M, G), barycentre projector | yes (rows genuine in fakes) | — | — | Thm 1.2 | dies |
| 12 | Nonvanishing from an exact period identity | seesaw/Rallis/Chowla-Selberg (not inspected in S2) | Kronecker-factorised Gram of monomials over E_R + Cauchy-Binet | yes as formulated: moments sum_{E_R} z^a zbar^b = N^b a_(a-b)(N), a product over p in the angles, copied by fakes [J] | yes | — | ZE^all_H, 5.3 | dies as formulated; see section 3 |
| 13 | Exact-norm chain (`ChainScalar`) from the integer cofactors t_ij | Rallis T T* = scalar | an integral operator built from t_ij with T T* scalar | unknown (t_ij is what fakes lack) | — | — | cocycle = "the problem itself" (XI.2), triangle identity uncapped (5.6) | no S2 content; see section 3 |
| 14 | CM theta jets on E = C/Z[i] at torsion points (cm "Programme C") | `ThetaZeroJetsCurve`, `ThmAGauss` | jets of th0 at 1/z_i = zbar_i/N mod Z[i] | not tested | supersingular analogue exists | [J] yes in a field of degree ~N^2: the torsion values live in ray class fields of conductor ~N, so a product-formula argument is Liouville in huge degree | not in record | dies on L [J, heuristic] |
| 15 | Sign-definite moments over Hecke characters ("iterated Parseval") | footnote 3 of the remark; not in S2 | sum_k F_K-weighted |a_(4k)(N)|^(2j) | yes (identity holds for angles) | yes | — | overview 5.4 rows 5-6 | dies; exact diagonal C(2j,j)^s caps it at M <= 2^(s/2) |
| 16 | Positivity of G on minors / Sym powers | `SymGramPosDef` | 2x2, 3x3, hypermetric conditions on exponent distances | automatic for every profile | — | — | ref8_extend, ref9_restate | dies |
| 17 | Behrend thick sets, long Iwasawa orbits | `BehrendSets` | choose S or sub-parents avoiding coincidences | constants only | — | — | random-parent rate; K_C | dies (our Galois action (Z/2)^s has orbits of length 2) |
| 18 | Siegel mass / Chowla-Selberg / rank-1 Siegel-Weil | mass ratios, CM periods | r(L,T) for T = (z_i zbar_j) | vacuous (genus of Z[i] has one class, r = 4) | — | — | 5.2 for averages over N | dies |
| 19 | Frame-to-span lemma (`FrameLA`) | near-tight frames => small corank | averaged slope over evaluation maps twisted by units | family is fake-visible | — | — | — | dies (tool without a family) |

**Checks behind rows 7 and 9 [J].**
* `judge/t1_inert_levels.py` (exact integers): the class count of {z mod l^k : z zbar = N} is (l+1) l^(k-1) for
  l = 3, 7, 11 and k = 1, 2. The all-level inequality of row 7 held in 9,600 tests on four norms (consecutive-in-angle
  and random subsets, M = 6..32, eight inert l): 0 violations, minimum slack 0 (tight). The constant
  sum_{l = 3 mod 4} log l/(l^2-1) = 0.2287, against a level-1 mass of 0.52, 1.57, 2.72 at M = 10, 100, 1000. So levels
  k >= 2 add a bounded multiple of M^2 to sum log|c_ij|, which is O(M) in W: the record's line 656 is confirmed.
* `judge/t2_conic_det.py` (sympy, exact Gaussian integers): det[z_j^a (0 <= a <= d), zbar_j^b (1 <= b <= d)] =
  +- N^(d(d+1)/2) V(z) / prod z_j^d held in 54 of 54 cases (M = 3, 5, 7; three norms). Hence |V|^2 >= N^(d^2), and the
  minimum observed log10 margin was 2.16 >= 0. Algebraically [J] this is implied by the cut identity (for odd M every
  (2S - M)^2 >= 1, giving prod |z_i - z_j|^8 >= N^((M-1)^2)), so it is the parity part of the cut identity and adds
  nothing.

## 3. Survivors

None survives all four tests. Two items are not killed outright and are recorded with what would have to be proved.

**(a) Nonvanishing from an exact global identity instead of a zero estimate (row 12; architecture B3, slopes #2).**
*What would have to be proved:* an auxiliary determinant Delta(cluster), with route-5 forced divisibility exceeding
its archimedean size once W = o(M log M), together with an identity that forces Delta(cluster) != 0 *at this cluster*.
That second half is exactly ZE^all_H ("the base locus of the forced space contains no nondegenerate point"). *Why fakes
might not satisfy it:* fakes have no global points, so an identity that uses global integrality or modularity in N is
not copied. *What kills the natural version:* the only identities in sight (Kronecker/Cauchy-Binet over E_R) depend
on the angles through products over p, which fakes copy [J], and they give nonvanishing of *some* minor, not of the
cluster's. Using modularity in N turns it into an averaged analytic input, which grafts cap. *First step (one session,
optional):* write the Gram matrix of monomials z^a zbar^b (a + b <= D) over E_R as a Kronecker product over p and
confirm, on one Paley/Legendre fake profile (`round8/lw.md`, `round31/README_VIII.md`), that the same factorisation
holds. The expected outcome is yes, which closes it. *Probability:* about 0.3% for o(log R/loglog R) in some class,
below 0.1% for moving the general rate. It is a renaming of route 5's open lemma, not a new route.

**(b) An exact-norm operator built from the integer cofactors t_ij (row 13; lifts #2).** *What would have to be
proved:* integral maps T_j built from the t_ij (z_j - z_i = 2 i eta H_ij t_ij) with T_j T_j* = c_j exactly and forced
divisibility g_j, such that sum_j (g_j - log D_j - (dim/2) log c_j) grows faster than the main term. *Why fakes might
not satisfy it:* the exact integer cocycle is the one input fakes lack (expert note 5.6, round 31 XI.2). *Problem:*
round 31 XI.2 proves that this cocycle *is* the genuine problem on a longer arc, so any such operator is an attack on
the problem itself; S2 contributes only the trivial algebra (`ChainScalar`). The obvious candidates (Fourier on angles,
Gauss sums on residues) satisfy Parseval for fakes. *First step:* at M = 3, check by hand whether any 2 x 2 or 3 x 3
integral matrix in t_12, t_23, t_13, A_ij satisfying t12 A12 + t23 A23 = t13 A13 has T T* scalar with a factor not
already in the pair data. *Probability:* below 0.3% (general rate below 0.1%).

## 4. Dead ends worth recording (so they are not retried)

1. **One-dimensional multiplicity is the cut identity.** On the norm-one torus (or P^1, or the conic), any
   interpolation determinant with multiplicity a at the M points is V^(a^2) times factorials, so its product-formula
   inequality is the a^2-th power of the cut identity. The norm-one torus at inert l has a height-1 formal group, so the
   S2 mechanism (D^p = 0 in a height-2 formal group) has no counterpart. [A:cm Prop B, A:slopes Prop S; standard]
2. **The Bombieri-Pila conic determinant is the parity part of the cut identity**: det = +- N^(d(d+1)/2) V / prod z^d,
   giving |V|^2 >= N^((M-1)^2/4). [J verified, 54/54]
3. **All inert levels l^k add O(M) to W** (constant 0.2287 * M^2/2 in sum log|c_ij|); S2's l/(l^2-1) is this.
   [J verified]
4. **Powers and products of the cluster are not critical** (length k C R^(k-1/2) against R^(k/2)), and the tower
   z -> z^k is homogeneous. [A:positivity 5.1, A:lifts 2.2]
5. **Sign-definite Hecke-character moments** are capped by the cube energy C(2j,j)^s: at best M <= 2^(s/2). [A:positivity 2.2]
6. **Averaging over twists mod tau, sign flips, lambda^(4k), or orders/degrees** is P3-char, profile-only, or the
   record's large-sieve problem. [A:lifts, A:stepanov]
7. **Even granting S2's conclusion (GRH)**, nothing follows for q_P: the residual problem concerns signed subset sums
   of s fixed Gaussian-prime angles at precision about exp(-s log P/4), far below the P^(-1/2) equidistribution
   that GRH gives. [A:lifts; J agrees]
8. The remark's "2/3" is not in S2, and its coincidence with our 2/3 (= 2/c, c = 3) is numerical only.

## 5. Recommendation

None. No experiment or lemma from this development is worth a session aimed at the rate. If anything is done, write
items 1-3 of section 4 into the expert note's section 5 as barrier entries ("Lean RH development, round 32: its
inert-multiplicity and slope devices reduce, on one circle, to the cut identity and the all-level collision count").
The single pointer of `round31/README_XI.md` section 5(iv) (the explicit parents q_P, for a specialist in linear forms
in logarithms) is unchanged.
