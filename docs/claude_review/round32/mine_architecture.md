# Round 32, block "architecture": the top-level skeleton of the Lean GRH development, and what transfers to arcs

Labels: **[Read]** = read in S2 (path given); **[Inferred]** = my reading of definitions whose names/fields I saw but
whose proofs I did not read; **[Heuristic]**, **[Speculative]** as usual. The Lean has **no comments or docstrings
at all** (grep for `/--` over Formalism/ and Development/ returns 0 files), so every interpretation of a field name is
at best Inferred. Nobody has compiled it; nothing below certifies that any proof is correct. Root of S2 below:
`S2 = .../scratchpad/dl_lean/x/lean/src`.

## 0. Results in five lines

1. The final theorem is GRH for Dirichlet L-functions in the form "L(s,chi) != 0 for Re s > 1/2"
   (`UpperTheory.EveryCharacter`), and it is derived **only** through a classical Landau-type lemma: a one-sided
   bound `G_Lambda(y) <= y^(1/2 + O(eps))` for a nonnegative-coefficient Cesaro-weighted Chebyshev error implies
   zero-freeness. [Read]
2. The one-sided Chebyshev bound comes from a single averaged inequality `umean ASK <= umean Main + ErrX`
   (`MeanBoundSecondForm`), where ASK is a sum `sum_p log p * sum_i (a/p^i - w)_+` over inert primes (a p-adic
   divisibility count) and Main is the archimedean size computed by Stirling plus Chowla-Selberg
   (`log Omega = gamma + log(pi/v) - L'/L(1)`). It is a product-formula (Liouville) inequality for Hermitian lattices
   attached to iterated theta lifts, phrased through slopes. [Read for the statements; Inferred for the meaning]
3. An exact identity (`AskMainBalance`) says ASK - MAIN equals an average over characters psi of the Galois group of
   the twisted Chebyshev errors G_psi evaluated at points a/w. So "divisibility minus size" IS the explicit formula. [Read]
4. The parameter m (number of lifts) multiplies the signal (each lift contributes a full copy of ASK and Main) while
   the dominant error has one m-independent part (about x*b per point) and one part quadratic in m (about m^2*b);
   signal per point is m*b*G(x). Hence G(x) <~ x/m + m, optimal at m ~ sqrt(x), exponent 1/2. The concrete
   parameter file confirms m ~ x^(1/2 - eps/16), group order n ~ x^(2+8 eps), weights b = exp(x^eps), a = x*b. [Read
   for the parameters; Inferred for the balance]
5. The "2/3" of the timeline does not occur in the Lean (the only literal 2/3's are `gamma < 2/3` used in
   bookkeeping). Its natural reading is a zero-free half-plane Re s > 2/3 from a cost x/m + m^2. For circles the
   coincidence with our 2/3 (a constant in front of log R/loglog R) is **only numerical**. Nothing in the block
   transfers to a proof of M = o(log R/loglog R); one meta-idea (replace a zero estimate by a computed nonzero
   period, and gain multiplicity by an orbit-glued family) is recorded with a first step and a probability below 1%.

## 1. The skeleton, top down

### 1.1 Top level [Read]

* `S2/Main.lean`: `riemannHypothesis := rhOfAll`, `LFunction_ne_zero_of_half_lt_re := allChi`, etc.
* `S2/Statements.lean`: `DirichletLNonvanishingRightOfHalf : forall N chi s, 1/2 < Re s -> (chi != 1 or s != 1) ->
  L(s,chi) != 0`; also zeros-on-the-line forms.
* `S2/Supplements/Kernels/MainTheorems.lean`: `theorem allChi : UpperTheory.EveryCharacter :=
  UpperTheory.FinalTheorems.everyCharacter`; `rhOfAll := Formalism.VW.rh_of_every allChi`;
  `critAll := Formalism.VW.critLine allChi ...` (strip version via the functional equation, `ClaimStrip.lean`:
  `strip (h : HalfPlane) : Strip`, `zeta (h : HalfPlane) : Zeta`, using `prim_reflect`, `gamma_ne_zero`).
* `UpperTheory.EveryCharacter` (`Development/AutomorphicArithmetic/CharacterBoxes/EveryCharacter.lean`) is
  literally the same statement as `Claim.HalfPlane` (`HalfPlane.lean`).

### 1.2 From EveryCharacter to the "literature" and to VWin [Read]

`Development/AutomorphicArithmetic/CharacterBoxes/MainConclusion.lean`:
```
theorem UpperTheory.FinalTheorems.everyCharacter_of_row (h : UpperTheory.GalStar.Row0b ...) : EveryCharacter :=
  every_character_of_literature_of_vWin ⟨h, recTExists, h⟩ (UpperTheory.One.vwin_of_row h)
```
* `every_character_of_literature_of_vWin (O : OpenLiteratureStatements) (hV : VWin)`
  (`.../ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean`, 814 lines). `O` packages three
  open "literature" inputs (a Caraiani-type local-global statement `M.Row`, existence of local Langlands
  reciprocity `RecTExists`, Galois representations of polarized pairs). Dozens of further big theorems (CM main
  theorem, Hecke-Tate FE, Deligne-Ribet, Siegel-Klingen, Bhatt-Scholze, arithmetic Hilbert-Samuel, principalisation,
  Kawamata-Viehweg, ...) are fed in as `theorem ... := <name>_holds`, i.e. claimed formalised in Development/.
  I could not check any of them.
* It calls `UpperTheory.Lit.every_character_of_theta_integrals_local_matching`
  (`.../ThirdPlaceSeries/SiegelEisensteinSecondTermHypotheses.lean`), whose hypotheses are
  `UnitaryGroupThirdPlaceHypotheses`, `ThetaLocalInputs`, `ThetaIntegralsAgree`, `LocalMatching`,
  `CovolumeBoundOfX`, `RecTExists`, `VWin`; it produces `Nonempty (Hyps ...)` (test vectors, Fourier-coefficient fit,
  evaluation grids: the local theta-correspondence bookkeeping) and calls
  `UpperTheory.Top.every_character_of_green_v_win`, which goes to `plugAll`/`plugWin`
  (`CharacterBoxes/EveryCharacter.lean`) and then to `zeroFreeAll_of_topG` (`CharacterBoxes/ZeroFreeAllOfTopG.lean`)
  with `la = 3/2`, and finally `UpperTheory.Top.imprimitive` (reduction to primitive characters).

**VWin = the "weakening of Vandiver" (timeline item 17).** `CharacterBoxes/ExistsPageN.lean` [Read]:
`VVle E t B` says: for every prime p and 0 < e < t the Selmer/Sha-type length `shaLen p e` of
H^1(o_E[1/p], Q_p/Z_p(-e))-type data is finite, and `sum_{l in S} log l * sum_{1<=e<t} shaLen_l(e) <= B` for every
finite set S of primes. `VWin` asks this, with B = sqrt(exp(x^eps)) = b^{1/2} and t <= x^{401}, for every field
E in an explicit family (`IsChiField psi e q d E`: subfields of Q(zeta_N) cut out by `a = +-1 mod 2^{e+2}`,
`psi(a) = 1`, `(a/q) = (a/d) = 1`), with q = 3 mod 8, d = 7 mod 16 primes, q | d+1, every odd prime <= Y^{eps/5}
dividing d+1, log q ~ Y^{eps/5}, log d ~ Y^{eps/5} log^2 Y (`WinPair`), and 2^e ~ x^{2+8 eps} (`PageN`). So the
total log-size of the cohomology classes that cannot be killed by twisting is at most b^{1/2}, negligible against
the signal (of order m*n*b*x^{1/2}). In `HypV` (`DivinNOfHypV.lean`) the same bound is imposed on each member of
the family. [Read; the role "negligible against the signal" is Inferred]

### 1.3 The analytic layer: Landau [Read]

`Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/TopWin.lean`:
```
theorem zeroFreeOn_of_tiedWin_cyc {Q : PChar → Prop} {la : ℝ} (hla : 1 < la)
  (h112 : ∀ ψ, Q ψ → ChebyshevInequalityOnWindows @core @Tied (Level.cyc ψ)) (hB : ∀ ψ, Q ψ → (Level.cyc ψ).FormB)
  (h : ∀ ψ, Q ψ → ChebyshevInequalityAtWindowScale @core @Tied la (Level.cyc ψ)) : ZeroFreeOn Q
```
* `Level.cyc psi` is the cyclic group generated by a primitive chi; `Level.Gm y = Gz y + sum_{psi in Nt} wt psi *
  Re gChi(psi, y)` (`ExplicitFormulaErrorOfHostConditions.lean`), with `Gz y = sum_{n<=y} Lambda(n)(y/n - 1) -
  y log y + (1+gamma) y` (`Formalism/StatementOC.lean`). So Gm is the Cesaro-weighted Chebyshev error of the
  Dedekind zeta of the cyclic field cut out by chi (coefficients >= 0).
* `ChebyshevInequalityOnWindows`: whenever a family `Fm` satisfies `Tied` on its windows,
  `Gm x <= x^(1/2 + (eD + eL + 7) eps)` there. `FormB`/`zero_free_of_Chebyshev_growth` (`Landau.lean`) is Landau's
  lemma: an upper bound `Gm y <= y^theta` on windows [Y, Y^la] (la > 1) implies no zeros in Re s > theta.
  Letting eps -> 0 gives Re s > 1/2.
* `Formalism/Statement*.lean` are all forms of this one classical step: `LandauStep`, `RHofG` (Statement.lean:
  `-C x^(1/2+eps) <= G_chi(x)` for G built from `Lambda(n)(1 + Re chi(n))`), `OCStep`/`OCAvg` (StatementOC: a
  one-sided bound on `Gz - Gchi` or a weighted average over a family forces equality of zero orders of zeta and
  L(chi)), `FamilyLandauC`/`FamilyLandauCW3` (StatementFL: one-sided bound on `Gz - Re sum c_i G_{chi_i}` gives
  zero-freeness of zeta), `ZetaLandau*`, `SampleAverageApprox`, `RiemannSumApprox` (StatementU: smoothing in x),
  and `ForcedSumCM`, `TwoTypes`, `HalfTypes` (StatementCM: for a CM type Phi in (Z/l)^* and its dual type, the
  class-restricted sums `sum Lambda(n)(a/n - b)` over n in the complement of Phi plus the dual add up to
  `b * Sz(a/b)` with error `<= 3a/l log(a/b)(1+log(a/b)) + a log l/(l-1)`; this is how "inert primes only" is
  completed to all primes). None of the Statement files contains the hard step; they are the (standard) analytic
  shell. [Read]

### 1.4 The quantitative heart [Read for statements, Inferred for meaning]

`Tied = ChainTied` (`Formalism/SharpRestriction/GlobalSchurBounds/TieGlobal.lean`) only ties the concrete family
to an abstract `ChainCore` (`Formalism/SizeWindowBrackets/SChain.lean`, line 60), which bundles: a finite abelian
group H, a `Data H W Top m f` record, constants, and the hypotheses `SizeData₂ Dd C Xpp` (=
`HidaSetting ∧ DIV ∧ EComp5 ∧ HComp ∧ Count5 ∧ Slopes ∧ QuotientSlopeBoundsWith ∧ LossBounds ∧ SlopeIneq ∧ HHBound ∧
FParams`), `ErrWindows`, `BNBounds`, `GrowthInequalityRemainingHyps`, `CostWindows`, `RestSize`. Then
`ChainAt : ... ChainInputs Cd eD eL eps x -> Gz x <= x^(1/2 + (eD+eL+6) eps)`.

The record `Data` (`Formalism/RankGapBounds/MeanBoundSecondForm.lean`) has, for each u0 in H and each "word" u:
a rank `d`, slopes `muN, muMaxN` of a lattice N_u, slopes `muM, muMaxM` of a comparison lattice M_u, a height `h`,
a dual slope `muCC` with minimum `muMinCC`, prime sizes `T`, local error sums `EeJ` per lift j, and the two
character-sum quantities `ASK u0` and `Main u0`. The single inequality from which everything follows is
```
def MeanBoundSecondForm : Prop :=
  Nonempty U0 → D.HidaSetting C → D.DIV C → D.EComp C → D.HComp C → D.Count C → D.Slopes C → D.SlopesN C →
    D.Losses C → D.SlopeIneq → D.HHBound C → D.FParams →  umean D.ASK ≤ umean D.Main + D.ErrX C
```
with (all Read, same file):
* `Main u0 = sum_{j<m} sum_{t in Top} a_t (log a_t - 1 - gamma + LpL - (1/f) sum_rho log w(u0,j,rho))`, where
  `LpL` = L'/L(1, chi_K) and `logOm = gamma + log(pi/v) - LpL` (Chowla-Selberg / Lerch: the log of the CM period).
* `DIV`: `muCC u0 + h(u0,u) - Fo u0 <= sum_t [ m (a_t log a_t - a_t - a_t logOm) - a_t log T(u0,u) ] + cDiv nT a La
  + cDiv m (B0 + Nc) x^(eDiv eps)`. The bracket is m times log of `a!/Omega^a` (Stirling): the archimedean size of
  the m lifts of a weight-a CM derivative value. [Inferred meaning]
* `SlopesN`/`QuotientSlopeBoundsWith`: `muM + Ga + Fo - Ee - Lo <= muN` and `muMaxN <= muMaxM + Ga + Fo - Ee + ...`:
  the lattice N_u of the iterated lift is the comparison lattice M_u shifted by the p-adic gain `Ga` (a quotient /
  "division" of slopes). `LossBounds`: `ASK u0 - Ga u0 <= LoP` (the gain realises the divisibility count).
* `SlopeIneq`: `s(u0) - log|W| - (1-del)^{-1}(4 uF log(1 + HH/uF) + (eta + 4 sqrt e) HH) <= muCC(u0)`, with
  `s = weighted mean of (muN - h)`: a slope (Bost/Harder-Narasimhan type) inequality for a nonzero map. [Inferred]
* Hence, schematically, `Ga <= (archimedean size) + errors`, i.e. **divisibility <= size**: a nonzero integral
  object cannot be more divisible than it is large. This is the product formula, made quantitative by Arakelov
  slopes, and averaged over u0 in H (the inequality is for `umean`, not for each u0).

The bridge from ASK to the Chebyshev function (`Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean`) [Read]:
```
def ASK (D : ArithData H) (N : ℕ) (a : ℝ) (t : H) (w : E0 × H → ℝ) : ℝ :=
  (∑ p ∈ (range (N + 1)).filter (IsInert D), 2 * Real.log p * askP D N a t w p) + Real.log D.d * askD D N a t w
-- askP = ∑_{i ≤ N, i•Frob_p ∉ Φ} max (a / p^i - w (partner t (i • Frob_p))) 0
def AskMainBalance : Prop := ∀ D, HostFacts D → ∀ ell a t w, 0 < a → ... →
      (ASK D N a t w : ℂ) - MAIN ell a w = De D ell a t w
-- De = (1/|Φ|) ∑_{g ∉ Φ} ∑_ψ conj(ψ g) * w(partner t g) * GPsi D ell ψ (a / w(partner t g))
```
`GPsi psi y` is the twisted Chebyshev error `sum_{n<=y} Lambda(n) primVal(psi,n)(y/n - 1) + y*ell(psi)` (and `Gz` for
psi = 0). So: the divisibility count minus the main term is exactly a Fourier average over the characters psi of
E0 x H of Chebyshev errors at the points a/w. Bounding the left side above therefore bounds a positive combination
of G_psi above, which is what Landau needs.

### 1.5 The parameters, and the roles of m, X, H [Read for numbers]

`SizesTextAt eps x` (`SChain.lean`) gives a model choice: H = Z/n with `x^(2+8eps) <= n <= 2 x^(2+8eps)`;
`x^(1/2 - eps/16) <= m <= 2 x^(1/2 - eps/16)`; `b = exp(x^eps)`, `a = x b`; `k = x^(9/2 + 10 eps)`;
`|TT| ~ x^(-eps/8) n/2`; `B0 = 2 n m b x^eps`. `RestSize` requires `x^(1/2 - eps) <= m` and
`k <= n^2 x^(1/2-eps)`. `Ranks` (`Formalism/LocalRingCardinality/Params.lean`) has k, m even, t odd, `m+3 <= k`,
`k+m+4 <= t`, `3m < t`, `m^2 <= k + m(t+1)`. The Hecke parameters of the lifts (`SizeData₂` in
`LocalRingCardinality/DataParams₂.lean`) have the initial weight `u_0 = a_t + (k-1)/2` (about x*b) and lift weights
`w_j = e_{j-1}` in `[b, b + etaSz + 1]`.

* **X and m.** X is the length x of the von Mangoldt sum; m ~ sqrt(x) is the number of iterated lifts. In
  `Main` and in `Dd.ASK u0 = sum_t sum_{j<m} ASK(a_t, t, w_{u0,j})` (`Ident`) every lift contributes one copy, so
  the signal per point of H is about `m * b * G(x)` (because `b Sz(a/b) = a(log(a/b) - 1 - gamma) + b Gz(a/b)`,
  `ForcedSumIdentity`, with a/b = x). The error `ErrX` contains `cDiv * nT * a * La` (no factor m: about
  n*x*b*x^eps) and `cDiv * m * (B0 + Nc) x^(eDiv eps)` (about n*m^2*b*x^(2 eps)). Dividing by `m * nT * b` gives
  `G(x) <~ x^{1+O(eps)}/m + m x^{O(eps)}`, minimised at m = x^{1/2}: this is the "sqrt X" of the remark and the
  "each lift wins one unit of multiplicity". [Inferred from the displayed formulas]
* **The Iwasawa tower.** H = Z/2^e, 2^e ~ x^{2+8eps}, is `OrdGrp s = (Z/2^{s+2})^* / {+-1}`
  (`Formalism/BehrendSets/Statement.lean`), the Galois group of the layer Q(zeta_{2^{s+2}})^+ of the cyclotomic
  Z_2-tower; the fields of `IsChiField` are its composita with the fields of psi, q, d. Lemmas `OrdIff`, `OrdLargeJ`
  (`orderOf (frob s l) <= J -> 2^{s+2} <= J(l+1)`), `ShortOrdSq`, `ShortOrdLin` show every prime l <= x has
  Frobenius of order >= 2^s/(x+1) in H, i.e. long Frobenius orbits: the analogue of passing from F_q to F_{q^n} in
  Stepanov-Bombieri. [Read for the lemmas, Inferred for the role]
* **Behrend.** `BehrendFree`, `BehrendBound`, `ThickSet`, `SpaceFree`, `Spacing` (`BehrendSets/Statement.lean`,
  `StatementE.lean`) [Read]: a Behrend-type set S of density `exp(-C sqrt(log n log L))` with no small-coefficient
  three-term relation, so that any segment `{t + i h : j0 <= i <= j0 + Lw}` of a Frobenius orbit meets the "thick"
  set in at most 2 points. It is used to choose the set TT of top places (|TT| >= x^{-eps/8} n/2, `RestSize`) so that
  the partner points `partner(t, i * Frob_p)` in ASK do not collide too often. [Inferred]
* **"Division"/"gcd glued across Hecke characters" (items 14-15).** The closest Lean object is the LP deficit
  `lpF` (`DataParams₂.lean`) [Read]: for a prime l, an orbit length d and budgets `Ds : Z/d -> R`,
  `lpF l d Ds = l/(l^2-1) * sum_s Ds(s) - sup { sum_{s,i} x(s,i) : x >= 0, x(s,i) = 0 unless i odd and i < d,
  sum_i l^i x(s,i) + sum_i x(s-i, i) <= Ds(s) for all s }`. Here `l/(l^2-1) = sum_{i odd} l^{-i}` is the ideal
  inert divisibility `a/p + a/p^3 + ...` of timeline item 2, and the shift `s - i` couples neighbouring characters
  along a Frobenius orbit. `FOrb j p z = lpF p (ordP p) (s -> (k-1) w_j(pos p z s))`
  (`LocalToGlobal/Divisibility/Numbers.lean`) is the loss of lift j at p, averaged over the orbit (`FInert`).
  So the divisibility actually realised by the glued integral structure is the LP optimum over an orbit of Hecke
  characters, a gcd-type quantity, and the averaging over u0 in H in `MeanBoundSecondForm` is the "extra averaging".
  [Inferred]
* **2/3.** Not found. Literal 2/3's in the relevant files are `Real.eulerMascheroniConstant_lt_two_thirds`
  (e.g. `PChain.lean`, `DIVfromGLf.lean`, `BudgetOfSizeAndWeight.lean`). My reading of item 14: if the per-lift
  cost were of order m^2 b per point instead of m b (no division by the common divisibility of earlier lifts), the
  balance would be `x/m + m^2`, i.e. m = x^{1/3} and a zero-free region Re s > 2/3. [Heuristic; not located in S2]
* **Not inspected** (budget): the seesaw/Bessel-period step (files exist, e.g.
  `Development/AutomorphicSpectrum/Proofs/MembershipProofs/SeesawThetaDualPair.lean`), GGP determinants, the
  relative trace formula, Chowla-Selberg/Siegel-mass norm computations. I make no claim about them.

### 1.6 Assessment of the skeleton [Inferred/Heuristic]

The architecture is: Landau shell + explicit identity (divisibility - size = character average of G_psi) +
one product-formula inequality for lattices of iterated theta lifts, with multiplicity m ~ sqrt(x). Structurally
it is Chebyshev's 1852 method (log a! = sum Lambda(n) floor(a/n), compared with Stirling), upgraded: the integer a!
is replaced by CM-period-normalised derivative values (Hurwitz-type numbers, divisible by p^{a/p+a/p^3+...} for p
inert), and Chebyshev's constant-factor loss is beaten by gaining m-fold divisibility at a cost growing slower than
the gain. All the difficulty is inside the hypotheses DIV, SlopesN, SlopeIneq, LossBounds (i.e. that the iterated
lifts exist, are nonzero, have the stated integral structure and divisibility), which I have not checked and which
are where an error would most plausibly hide. A sceptical note: the claimed method gives an *upper* bound for Gz
from a *lower* bound on divisibility; a Chebyshev-type method with such huge weights (b = exp(x^eps)) and an error
x/m + m needs the divisibility to be exact to relative precision x^{-1/2} per lift, which is an extraordinarily
strong claim about p-adic properties of the lifts at all p <= x simultaneously.

## 2. Translation to lattice points on short arcs

### 2.1 Is there a structural analogue of the 2/3 and of "division"? 

Our 2/3: `3 M log M <= W + O(M)` (Theorem 1.1) from the cut identity
`prod_{i<j}|z_i - z_j|^8 = N^{M(M-2)} prod_p p^{sum_tau (2S_{p tau} - M)^2} prod_{i<j}|c_ij|^8`, an archimedean upper
bound for chords and pigeonhole collisions at primes q <~ M (mass `sum_{q<=M} log q/q ~ log M`). This is exactly the
same *type* of inequality as `ASK <= Main` (divisibility <= size for a nonzero integral product), at multiplicity
one. [Inferred]

* The two 2/3's are different kinds of number: theirs (if my reading is right) is an exponent in a zero-free
  half-plane from a balance x/m + m^2; ours is the constant 2/(3) = 2/c in M <= (2/c) log R/loglog R with c = 3
  coming from how the collision mass log M is shared between the (2S - M)^2 cut terms and the reduced chords.
  The coincidence is **numerical only**. [Heuristic, confidence high]
* "Division by a gcd" in our problem already exists and is used: normalisation by the Gaussian gcd of the cluster
  and the reduced chords `c_ij = (z_i - z_j)/gcd(z_i, z_j)` are built into the cut identity. Expert note §5 item 7:
  gcd normalisation is invisible to forced-divisibility, auxiliary-form, Liouville, pair, cut, character and
  single-form inequalities. Their division is of a different kind: it removes divisibility already credited to
  earlier lifts in an m-step chain. We have no chain whose cost grows with its length, so there is nothing to divide.
* "Gluing across Hecke characters along a Frobenius orbit": the closest object for circles is the family of the
  2^s lattice points of the same norm, indexed by sign vectors sigma (z_sigma = unit * prod_p pi_p^{(1+sigma_p)/2 e_p}
  conj(pi_p)^{...}), together with the characters lambda^{4k}. The Frobenius-orbit coupling `x(s - i, i)` has no
  counterpart: for circles the "Galois group" acting on E_R is (Z/2)^s (complex conjugation prime by prime),
  whose orbits have length 2, not long cyclic orbits. The long-orbit trick (Iwasawa layer so that all primes <= x
  have long Frobenius orbits) uses that the primes p <= x are *fixed* and the tower level n -> infinity; in our
  problem the primes are those of N and there is no auxiliary tower in which they become inert with long orbits
  that preserves the arc condition. [Inferred]

### 2.2 Candidate borrowable ideas, tested

**B1. Explicit identity "divisibility count - archimedean main term = character average of twisted sums"
(`AskMainBalance`).** Analogue: write the collision mass in the chord product,
`sum_{i<j} log|c_ij|` restricted to collisions modulo inert q (Theorem 1.2's mechanism), as
`sum_q log q * #{(i,j): z_i = z_j mod q-level}` and expand the collision indicator in characters of the norm-one
torus T(F_q) (order q+1) and of Z[i]/(q^r): `#collisions = (1/(q+1)) sum_chi |sum_i chi(z_i/zbar_i)|^2`.
Main term = trivial character; the rest are twisted Weyl sums over the cluster. Test: these are rational-character
and residue inputs at moduli q <= M^A: capped by fakes (§5 item 1; P3-char capped at every modulus, round 31 IX.3).
Already in the record (`round5/construct.md`, `round31/README_IX.md`). **Fails the fake test.** Probability of rate:
< 0.2%.

**B2. Multiplicity by a chain whose gain is linear in m and whose cost is quadratic in m (balance x/m + m).**
Analogue: amplify the chord-product inequality m-fold. The obvious realisations: (a) powers z_i^k, k <= m: the arc
of the powered points has length about k C R^{k-1/2}, far above the critical scale R^{k/2}, so each power costs
(k-1)/2 * log N per pair while the extra collisions (fibres of size gcd(k, q+1) in T(F_q)) gain O(log M) per
pair; dead at once [Inferred]. (b) Confluent/generalised Vandermonde determinants in monomials of the points
(Stepanov with multiplicity): this is route 5, whose missing ingredient is the zero estimate at the trivial
subgroup (§5 item 3; round 26 `ze.md`); single-place versions are capped at tau^* < 2. (c) Iterated theta lifts
for U(1) -> U(2) -> ... over Q(i) with the CM character lambda^{4k}: coefficients of theta series of Hermitian
Z[i]-lattices of rank r count r-tuples with prescribed Gram matrix; but r points on one circle have a rank-one Gram
matrix (z_i zbar_j z_j zbar_i = N^2), so the r-tuples of a cluster sit at singular (degenerate) Fourier
coefficients, where theta lifts carry Eisenstein/degenerate information only [Inferred]. **(a) dies on cost,
(b) is in the record, (c) has no nondegenerate place to land.** Probability of rate: < 0.3%.

**B3. Replace the zero estimate by a nonvanishing that is an exact period/inner-product identity.** In the
development, the nonvanishing of the auxiliary object is not proved by a zero lemma but obtained from a computed
norm (seesaw/Rallis inner product, Chowla-Selberg, Siegel mass) [timeline items 9, 12; not inspected in S2]. This
is the one idea that addresses exactly route 5's missing lemma. Analogue programme: find an auxiliary determinant
Delta(cluster), built from monomials in the points, such that (i) Delta is a Gaussian integer with forced
divisibility exceeding the archimedean bound when M >> log R/loglog R (route 5 already supplies this for suitable
designs), and (ii) Delta != 0 follows from an identity expressing a positive combination of |Delta(S)|^2 over
configurations S as an explicit nonzero product. The natural identity: E_R is a product set
(E_R = unit * prod_p {pi_p^a conj(pi_p)^{e_p - a}}), so moment (Gram) matrices of the whole circle in a monomial
basis factor as Kronecker products over p and have explicit nonzero determinants; Cauchy-Binet writes such a
determinant as a sum of |minor(S)|^2 over subsets S. Test: (1) this gives `sum_S |Delta(S)|^2 = explicit > 0`, i.e.
an upper bound on each term and nonvanishing of *some* minor, not of the minor at the cluster; we need the cluster's
minor nonzero, so the identity must be localised, and the only localisation available is the zero estimate again.
(2) Fakes: the product structure in p and the Kronecker factorisation hold verbatim for angle models (they keep
genuine primes and exponent profile), so any conclusion from this identity alone holds for fakes. (3) Char p: the
same factorisation holds over F_p[t]. **Fails all three tests as formulated.** It survives only if the period
identity uses the *modularity* in N (theta series of Z[i] as a modular form of weight 4k+1), which fakes do not
have; but modularity speaks about all N at once, and the record's averaged-analytic inputs built from it are capped
by grafts (§5 item 2). [Speculative] Probability of rate: about 0.5%.

**B4. Behrend-thick sets to control orbit collisions.** Analogue: choose a large subset of sign vectors / of
characters k avoiding short three-term relations. Our obstruction is not repeated collisions of a few orbits but
the size of signed subset sums at the Liouville threshold; nothing here. **Not transferable.**

### 2.3 Filters, for the record

* Characteristic p: the development is intrinsically characteristic 0 and archimedean (Stirling, Chowla-Selberg,
  slopes at infinity); in that sense its *mechanism* passes the filter. But none of the transfers above keeps that
  feature except B3 with modularity, which then becomes an averaged analytic input.
* Single-place Liouville: the core inequality `divisibility <= size` is single-place Liouville at multiplicity
  one; the whole content of the development is the m-fold amplification. Our attempts to amplify (B2) are either
  costlier than the gain or are route 5.

## 3. First steps (one session each) and probabilities

* B3 (only candidate worth a session). Lemma to prove or refute: for E_R with normalised N = prod p_j^{e_j} and the
  monomial basis {z^a zbar^b : a + b <= D}, compute the Gram matrix G_D = (sum_{z in E_R} z^{a} zbar^{b}
  conj(z^{a'} zbar^{b'}))_{...} as a Kronecker-type product over p, and decide whether det of the restriction to a
  route-5 design space factors with an explicit nonzero value; then check on the Paley/Legendre fakes of
  `round8/lw.md`, `round31/README_VIII.md` whether the same factorisation holds (prediction: yes, hence capped).
  Probability that this leads to M = o(log R/loglog R) even in a restricted class: about 0.5%.
* B1, B2, B4: no session recommended; each is either in the record or fails on cost.
* Overall: the architecture block yields no idea that gives the rate. Probability that anything in this block moves
  the worst-case 2/3: below 1%.
