# Block "stepanov": the Stepanov determinant and the role of multiplicity

Labels: **Read** = statement read in S2 (path given); **Inferred** = my reading of how statements fit
together; **Heuristic**; **Speculative**. S2 is uncompiled and unverified; nothing below treats its
statements as established mathematics. About 30 files were read (statements only).

## 0. Results in five lines

1. S2 contains no module, definition or comment named "Stepanov" (a grep over all `.lean` files returns
   nothing). `Development/WeilBound` does not hold a Weil bound for curves proved by an auxiliary
   function. It holds GL_n representation theory, Siegel–Weil and theta-lift material. The timeline's
   "Stepanov determinant" has no single counterpart in S2. Its ingredients are spread over three
   kinds of lemma. (Read: the grep and the directory listing.)
2. In S2, **nonvanishing of the determinant costs nothing.** It is the Gram determinant of linearly
   independent vectors under a positive-definite integral Hermitian (Fock) pairing, so `det > 0`
   (Read, `SmithDeterminants/StatementGram.lean`). **Multiplicity becomes ϖ-adic divisibility** through
   a Smith-form lemma: if `L∘y^m ≡ 0 mod ϖ^S`, then `ϖ^{dS}` divides `ϖ^{dmt + m·n_t·S}·det L` (Read,
   `StatementDVR.lean`, `TruncFinal`). The divisibility of `y^m` is gained one power per theta lift,
   along a chain `θ_0, …, θ_{m-1}` (Read, `LatticeChainBound`). There is no zero estimate anywhere. What
   takes its place is positivity, or the Bost slope inequality with a filtration by vanishing order
   (Read, `MonomialMeanResidues`, `DivisorInequality`).
3. `LogVandermondeBounds` does not bound Vandermondes of integer or Gaussian points. Its Vandermondes are
   **Weyl dimension formulas** (`dimVN x N = ∏_{i<j}(x_i-x_j)/∏_{i<N} i!`). The module is the bookkeeping
   ledger of how much dimension and mass each theta lift gains or costs (Read).
4. For our problem, nothing in the block yields an inequality sharper than Theorem 1.1. The sharp
   divisibility of a Vandermonde of points that S2 implicitly uses (superfactorial, or Bhargava,
   divisibility) is, on the circle, exactly multi-level pigeonhole collisions. Levels `l^k` with `k ≥ 2`
   change `W` by `O(M)` only, so the constant 3 stays 3. These are point-residue inputs, so fakes satisfy
   them.
5. The one structural idea worth stating is: "make nonvanishing free by positivity, and put all the
   arithmetic into divisibility gained by operators." The record already caps its simplest form: grid
   and Gram determinants with free nonvanishing are subcritical (expert note §5 item 3). It needs a
   high-rank family carrying operators with few large elementary divisors, and one circle does not
   supply such a family. Probability that anything here yields `M = o(log R/loglog R)`: ≤ 1%.

## 1. What S2 actually contains (statements)

### 1.1 The determinant: a Fock Gram matrix (Read)

`Formalism/SmithDeterminants/DefsGram.lean`, `StatementGram.lean`, namespace `Formalism.A`.

* `mfact g = ∏_v g(v)!` for a multi-index `g`. The Fock pairing on `Z[x_1..x_k]` is
  `fockZ f g = Σ_m f_m g_m m!`. On `Z[i][x]` it is `fockG f g = Σ_m f_m conj(g_m) m!`. This is the
  Bargmann–Fock inner product `∫ f conj(g) e^{-|x|^2}` up to normalisation.
* `famZ E α a = E_a · x^{α_a}`, where `E_a` is homogeneous of degree `e_a` and
  `deg α_a + e_a = s`.
* **`GramZ`** (the Gaussian version **`GramGauss`** is the same). Hypotheses: `E_a` homogeneous of degree
  `e_a`; `deg α_a + e_a = s`; `‖E_a‖_1 ≤ k^{e_a}`; the family `E_a x^{α_a}` linearly independent over `Z`.
  Conclusions: (i) `α_a!` divides every entry of row `a` of the Gram matrix `G`; (ii) `G` is positive
  definite and `det G > 0`; (iii) `det G ≤ ∏_a G_aa` (Hadamard); (iv) `G_aa ≤ α_a!·(k^2 s)^{e_a}`;
  (v) `∏_a α_a!` divides `det G`, and `1 ≤ det G/∏_a α_a! ≤ ∏_a (k^2 s)^{e_a}`; (vi) for every prime
  `p`, `v_p(det G) ≤ Σ_a (v_p(α_a!) + e_a log(k^2 s)/log p)`.

So the matrix is the Gram matrix, under the Fock pairing, of monomial multiples of a fixed family of
forms. The "multiplicity" is the exponent `α_a` of the monomial factor. The factorial `α_a!` is the
divisibility it buys. The upper bound is Hadamard. The ratio `det/∏α!` lies in `[1, ∏(k^2 s)^{e_a}]`, so
the lemma is a two-sided control of an integral structure. It produces no contradiction by itself.

### 1.2 Multiplicity becomes divisibility of a determinant (Read)

`Formalism/SmithDeterminants/StatementDVR.lean`, namespace `Formalism.A`. `O` is a DVR with uniformiser
`ϖ`. `SmithData y d` is a pair of bases with `y(b_i) = a_i b'_i`. `n_t := nAbove ϖ t` counts the `i` with
`a_i ∤ ϖ^t`.

* **`TruncFinal`**: for `y : M → M` with Smith data, `L : M → O^d`, and `L(y^m v) ∈ ϖ^S O^d` for all
  `v`, one has `ϖ^{dS} | ϖ^{dmt + m·n_t·S} · det(L)`. That is,
  `v(det L) ≥ dS − dmt − m n_t S`.
  In words: if the evaluation map `L` kills `y^m` modulo `ϖ^S`, and `y` has only `n_t` elementary
  divisors above `ϖ^t`, then `det L` is very divisible.
* `TruncDet`, `LengthQuotientBound`, `TwoTier`, `GradedHorn`, `SatChain` are the length estimates
  behind it, for example `length(M/(im y^m + ϖ^S M)) ≤ rank·m·t + m·n·S`.
* **`LatticeChainBound`**. Data: spaces `V_i`, lattices `M_i`, subspaces `E_i`, operators `y_i, T_i` on
  `V_i`, and maps `θ_i : V_i → V_{i+1}` that commute with `y` on `E`, with `θ_i M_i ⊂ M_{i+1}`. Hypotheses:
  `T_i = ϖ^{e_i} y_i + ϖ^{f_i} c_i` on `E_i`; `θ_i T_i M_i ⊂ ϖ^{s'_i} M_{i+1}`; `s'_i ≤ f_i`. Conclusion:
  `θ_{m-1}∘…∘θ_0 (y_0^m φ) ∈ ϖ^{Σ_{i<m}(s'_i − e_i)} M_m` for `φ ∈ M_0 ∩ E_0`.
  This is the formal content of "each lift wins one unit of multiplicity". It needs a Hecke-type
  operator `T_i`, congruent to a scalar modulo `ϖ^{f_i}`, whose difference from that scalar is
  `ϖ^{e_i}` times the derivation-type operator `y_i`, and lifts `θ_i` that make `T_i`-images
  ϖ-divisible. The algebra is two lines (`θ y v = ϖ^{-e}(θ T v − ϖ^f c θ v)`). All the content is in the
  hypothesis `θ_i T_i M_i ⊂ ϖ^{s'_i} M_{i+1}`.
* **Inferred**: a grep for `LatticeChainBound`, `TruncFinal` and `GramGauss` outside
  `SmithDeterminants` finds no importer. Either these statements are not used, or they are used under
  other names (for example `MahlerMeasure`'s `SmithZp`).

### 1.3 Averaged control of the bad elementary divisors (Read)

`Formalism/MahlerMeasure/StatementTM.lean`, **`TheoremM`** (re-exported in
`PolynomialDimension/FinalM.lean`). It is stated over `Q_p` with `p ≠ 2`, letters `LL` and tie data `T`
(degree parameters `n_i` running over arithmetic progressions `n_0 + dd·j`, `j < Y`, with
`gcd(dd, p) = 1`). Conclusion: the average over the progression of the number of elementary divisors
of valuation above `t(E + C_0κ) + C_0κ` of a p-adic Fock-type Gram pairing, on a space `Y` of
polynomials, is at most `Σ_i (2S_i/Y_i + (C_0κ + 2S_i lgV2/Y_i)/E)` times the average dimension.
**Inferred**: multiplicity is not available pointwise, but it is on average over a progression of
degrees. This is the role of the timeline's "pigeonhole on F_{q^n}/F_q" / Iwasawa tower: average
over many levels and find most of them good.

`HasOrder v e f` (Read, `MahlerMeasure/HasOrder.lean`) is the order of `f` in the variable `v`.
`GamDeg` and `GamLead` say the normalised pairing `tieGam`, as a function of the degree `n`, is a
polynomial of degree `≤ total + deg − max(e, e')`, with leading coefficient `e!·C(z+e, e)` times the Fock
pairing of the initial parts. So **vanishing order lowers the degree in `n`**. This is the
Stepanov-type mechanism in S2 that is most concrete: order of vanishing in one variable becomes
polynomiality of low degree in the level parameter, and averaging over a progression of levels then
bounds the bad set.

### 1.4 "Multiplicity at the identity" (Read, outside the listed modules)

`Formalism/SharpRestriction/GrowthBounds/CoarseMultiplicity/CoarseMultiplicityAtIdentity.lean`,
`CoarseMultiplicityAtIdentity`: there exists `F` in a polynomial space such that `w^{Σ blockOrd}·F` is
integral, `F` agrees with a prescribed symmetric family modulo `w^c`, and for every GL-translate `g'`
whose rational plane differs from that of the identity, `F` lies in the `w`-adic range `InWjRg … g' F c`.
`KeyOne.lean`, theorem `loc4One`, computes the substitution of a local polynomial as
`∏_j X_j^{p^{g_j} − 1} − w^{Σ_j(p^{g_j}−1)} ∏_j X_{j,c}^{p^{g_j}−1}`. **Inferred**: here an auxiliary
polynomial with Frobenius-type exponents `p^g − 1` vanishes `w`-adically to order `c` at all
non-identity points of a finite family. This is the nearest thing in S2 to "a zero estimate at the
trivial subgroup", and it is p-adic and Frobenius-driven, not archimedean.

### 1.5 Lower-bound protection: small denominators of projections (Read)

`Formalism/PolynomialDimension/StatementHar.lean` and `StatementNum.lean` (namespace
`Formalism.CountingSumsFrames`), on polynomials in `x_1..x_k, y_1..y_k` with the `sl_2` operators
built from `h = Σ x_i y_i` and the contraction:
* `BinomVal`: `v_p(C(n, m)) ≤ log_p n`, and the variants with `log_p(A+B+k)` and `log_p(8a+k−1)+1`.
* `FockProj`: `p^{⌊log_p(A+B+k)⌋}` times the harmonic projection of a Fock-integral `f` is
  Fock-integral, and `p^{2⌊log_p(A+B+k)⌋}` times each isotypic (level) projection is.
* `GaussProj`: the same with loss `p^{(k+1)⌊log_p(A+B+k)⌋}` for the Gauss-norm lattice.
* `FockGaussXY`: Fock and Gauss lattices agree up to the factor `A!` and `p^{(k−1)log_p A}`.
* `HarQuot`: `p^{(k+1)log}·f = w + h·G` with `w` harmonic and both integral.
* `FrameMultinomial`: `v_p(d!/((d−ΣM)!∏M_i!)) ≤ (r+1)(log_p(2a+k)+1)`, summed `≤ Nk(log_p(2a+k)+1)`.
* `WeylFlat`: a comparison of Weyl dimensions under adding one row.
**Inferred**: these lemmas show that decomposing into harmonic or isotypic pieces loses only
logarithmically many powers of `p`. The divisibility gained in §1.2 therefore survives the
representation-theoretic projections.

### 1.6 "LogVandermondeBounds": Weyl dimensions, not point Vandermondes (Read)

`Formalism/LogVandermondeBounds/LogRatioLowerBound.lean` defines
`vandN x N = ∏_{i<j<N}(x_i − x_j)` and `dimVN x N = vandN x N / ∏_{i<N} i!` (the Weyl dimension of the
GL_N representation with shifted highest weight `x`). `rowPar k v b` concatenates the `k−1` weights `v`
with the `m` new weights `−b_j` (the weights added by a theta lift). `PlaceData` encodes interlacing
(`v_i < u_i`, `u_{i+1} < v_i`) and the size windows.
* `VandPlace`: `log vandN(rowPar k v b, k−1+m) − log vandN(u, k) ≥ (k−1)Σ_j log b_j − 2k log(8x e^{x^ε}) − 2k(1+log k)`.
* `UpperPlace`, `LogRatioUpperBound`: `(1/(fk)) log(D*/d*) ≤ C m L`, with `L = x^ε`.
* `LogRatioLowerBound`: `(1/(fk)) log(D*/d*) ≥ avg_ρ Σ_j log(w_{ρj} vE/π) + m log(g_1 g_2) − C(L + m^2L/k)`.
  Here `D*` and `d*` are mass × index × Weyl dimension of the lifted and the original representations,
  with masses carrying `∏_M ((M−1)! (vE/π)^M)^f · L(M)`.
**Inferred**: this is the "cost of each lift" ledger (timeline item 14), comparing the size of a lifted
determinant against the original one. It contains **no** inequality for Vandermondes of integers or
algebraic integers under divisibility constraints. Its only arithmetic fact about Vandermondes is the
classical one: `∏_{i<j}(x_i − x_j)/∏ i!` is a positive integer for decreasing integer `x`.

### 1.7 Slope inequalities in place of zero estimates (Read)

* `Formalism/MonomialMeanResidues/MonomialMean.lean`: `HSClaim`, `HSJoint`, `GradedSum`, `ShareZ`,
  `SlopeLimit`, `SlopeCor`, `MinkowskiCount`. Graded pieces indexed by vanishing order `a < T−D` have
  rank `≤ C(a+D−1, D−1)δ`. The degree of each graded piece is bounded below by
  `C(a+D−1,D−1)δ(−J + (T−a−D−3)η − aλ − l_A)`. The total degree `deg M_T` dominates the sum. The output
  is a lower bound of the form `(T−1 choose D)δ(T h − … ) ≤ deg M_T`, followed by a limit `τ → ∞`
  (`SlopeLimit`) giving `δ'(sl − c + t) ≤ δ(h + t)`. **Inferred**: this is Bost's slope method
  (arithmetic Hilbert–Samuel against a filtration by order of vanishing) in Arakelov bookkeeping.
  `ZerosPerVariable` is Schwartz–Zippel: a nonzero `G` with `deg_i G ≤ e` vanishes on at most
  `q e |C|^{q−1}` points of `C^q`.
* `Formalism/DivisorInequality/StatementObj.lean`, `BoundsMu`: an upper bound for the maximal slope
  `μ_max` of the dual of a hermitian lattice over a CM field by `meanLog κ + Σ_j lcv_j/D` plus local index
  terms. `DivIneqTopAnchorU.lean` holds numeric constant checks (`DivisorIneqSharp`,
  `DivisorIneqFailsAtTwo`: the constant 184·Cmax² works for Cmax ≥ 11/5 and fails at 2). This is the
  "Bost slope inequality" of the timeline.

### 1.8 Remaining modules of the block (Read)

* `ResultantValuation`: `WindowKappaEq`: for `l` odd and `v_l(u − 1) = 1`,
  `v_l(∏_{i=b}^{b+S−1}(1 − u^i)) = S + Σ_{i=b}^{b+S−1} v_l(i)` (LTE). `WindowValuation`:
  `Σ_{b≤i<b+S} v_l(i) ≤ S/(l−1) + log_l(b+S)`. `ResultantConjugate`: for root sets closed under
  `α ↦ conj(α)^{-1}` up to permutation, `Res = (−1)^{ab}(∏α)^b(∏β)^a conj(Res)`. Also `VandermondeSplit`,
  `LaplaceRows`, `PermDetOnePlus`, and `SubspaceCount` (Gaussian binomials).
* `TruncatedSeriesCounting/StatementX.lean`: `TruncRem1`:
  `(π_0(π_0+ℓ)^{ℓ−1})^y ∈ (ℓ^c, π_0(π_0+ℓ)^K)` when `c ≤ ℓy − K` (Lubin–Tate / ramified tower algebra).
  `UnitCoefficientPowerMod`. `LemmaDetPoly`: coefficients of `det(Σ X_i a_i)` lie in the subalgebra
  generated by characteristic-polynomial coefficients. `PropEta`: a congruence-module length bound
  `≤ Σ_{x∈Zg} dim V_x·(vcc(x) + v(A))` for a Hecke-type operator acting by `A` on one eigen-part and by
  0 on the rest (Hida's η-invariant).
* `RankGapBounds`: parameter arithmetic for the ranks `(k, m, t, r, N)` of the unitary tower
  (stable-range and Siegel–Weil range inequalities). It also contains `UniformMeanBound`, a
  slope-and-mean bookkeeping. It has no Diophantine content.
* `LogFactorialBounds` and `LogFactorialIndices`: Stirling sandwiches for `ph(u) = u log u − u`,
  `Interlace` (a sum of `ph` over interlacing sequences changes by `≤ (w+2S)log(w+S+2)`), `HookFact`,
  `NormRel` (`det(YGY*) = σ(det Y · c det Y)·det G`), `LogProductFormula`, and `KWCancel*` (the cancelled
  exponent of the period `logΩ` after a lift is `(2j−1−m)·w`; **Inferred**: this is the "cancellation of
  order-m terms" of timeline item 9).

### 1.9 How the pieces fit (Inferred, low confidence)

The determinant is a Gram determinant of automorphic or theta vectors in a Fock or Gauss lattice.
* Upper bound: Hadamard together with masses and Weyl dimensions (§1.1, §1.6).
* Lower bound: divisibility. It comes from factorials (`α!`) and from the operator chain
  (§1.2, `LatticeChainBound`), where each of `m` theta lifts turns one power of the derivation `y`
  into `ϖ^{s'−e}`. It is protected by the logarithmic-loss projection lemmas (§1.5).
* Nonvanishing: positivity (`PosDef`).
* Missing pointwise control: replaced by averaging over levels (§1.3) and by slope inequalities (§1.7).

"Multiplicity m from m lifts" is therefore **ϖ-adic divisibility of `y^m`**, not archimedean vanishing
at a point. I found no lemma of the literal form "a section vanishing to order m at the relevant points,
with small norm, is zero". The nearest are `TruncFinal` combined with `GramZ` (positive determinant,
large ϖ-adic valuation, Hadamard upper bound), and the slope inequality (§1.7).

## 2. Route 5 against the S2 determinant, entry by entry

Route 5 (`round10/determinant.md` §1–2; `overview/overview.md` §6.2):
* Entries: `F_c(z) = Σ_{k∈A} c_k G_k(z)`, where `G_k(z) = N^{(D−|k|_1)/2} ∏_x z_x^{k_x^+} conj(z_x)^{k_x^-}`
  and `k` runs over a slab box `A ⊂ Z^n_0`. This is one number, not a determinant. The interpolation
  version uses `det(G_k(z^{(j)}))`.
* Multiplicity: `ord(c) ≥ m` means `Σ c_k P(k) = 0` for `deg P < m`. This is archimedean vanishing to
  order `m` at the trivial point `δ = 0`, so `|F_c(z)| ≤ N^{D/2}|c|_1(d|δ|)^m`.
* Arithmetic: forced Gaussian divisibility (`v_π(G_k(z)) = e_p D/2 + <k, a_p>`). This gives
  `Norm F_c(z) ≥ N^D ∏ p^{−h_supp(a_p)}` **if `F_c(z) ≠ 0`**.
* Missing: `F_c(z) ≠ 0` (ZE, MIN^pair, ZE^all_H).

| | S2 (Read/Inferred) | Route 5 |
|---|---|---|
| object | Gram det of lin. indep. Fock vectors | value (or interpolation det) of an exponential polynomial at the cluster |
| nonvanishing | free (positive definite) | open (ZE) |
| multiplicity | ϖ-adic: `L∘y^m ≡ 0 mod ϖ^S`, gained by m theta lifts | archimedean: order m at `δ = 0` |
| arithmetic lower bound | `∏ α_a!` and `ϖ^{Σ(s'−e)}` divide det | forced π-adic valuations at split primes |
| upper bound | Hadamard with masses and Weyl dimensions | Taylor: `|δ|^m` |
| dimension count | lifts: dimension gain against cost per lift (`VandPlace`, `KWCancel`) | slab box: `ord/width ≈ 0.12 sqrt(n/log r)` (supercritical) |

The decisive difference is **where the multiplicity lives**. In S2 multiplicity and arithmetic
live at the same non-archimedean place: the derivation `y` is ϖ-adically degenerate (timeline item 2:
`D^p ≡ 0` on a CM curve for inert p), so vanishing and divisibility are the same thing, and
nonvanishing is supplied separately by positivity. In route 5 multiplicity lives at ∞ at the point
`δ = 0`, while the arithmetic lives at the split primes at the actual point. Joining the two requires
the actual point to avoid the base locus, and that is the zero estimate.

## 3. Answers to the block's questions

### (a) What plays the role of the missing zero estimate in S2?

Nothing plays its role literally; the design avoids the need (Inferred from §1.1, §1.2, §1.7).
1. **Positivity.** The object is a Gram determinant of a linearly independent family, so it is positive
   (`GramZ` (ii)). The independence comes from the construction (distinct monomial supports, harmonic
   decomposition, `HarQuot`), not from a zero estimate.
2. **Arithmetic source of multiplicity.** `T = ϖ^e y + ϖ^f c` together with
   `θ T M ⊂ ϖ^{s'} M_{next}` gives the conclusion of `LatticeChainBound`. The order `m` equals the number
   of lifts. `TruncFinal` converts it into `v(det) ≥ dS − dmt − m n_t S`.
3. **Averaging over levels instead of pointwise control** (`TheoremM`), and the **slope inequality**
   with filtration by vanishing order (`HSClaim`, `SlopeLimit`).
So the answer to "is there a lemma 'section vanishing to order m with small norm is zero' with an
arithmetic source of the order of vanishing?" is: **only in aggregate form**, as `TruncFinal` combined
with `GramZ`, or as the slope inequality. The order of vanishing is ϖ-adic and comes from a Hecke
operator congruent to a scalar plus ϖ^e times a derivation.

### (b) Is there an analogue for M Gaussian integers on a short arc?

I tested four candidate sources of high-order vanishing.

**B1. Frobenius at inert primes (the honest analogue of `f(x^q)' = 0`).** Let `l ≡ 3 (4)` and `l ∤ N`.
Frobenius on `Z[i]/l` is conjugation. Hence every lattice point satisfies `z^{l+1} ≡ z z̄ = N (mod l)`,
and by lifting, `z^{(l+1)l^{k−1}} ≡ N^{l^{k−1}} (mod l^k)`. Equivalently, the ratios `u_i = z_i/z_1` lie
in the norm-one torus `T(Z_l) ≅ μ_{l+1} × Z_l`. So the polynomial `z^{(l+1)l^{k−1}} − N^{l^{k−1}}`
vanishes to l-adic order `k` at **every** point of the circle. This is genuine arithmetic "multiplicity"
(fakes have no l-adic structure). But on a one-dimensional point set the Smith-form lemma
`TruncFinal`, applied to `L` = evaluation at the cluster, returns exactly the l-adic valuation of the
Vandermonde `∏(z_i − z_j)`, that is, pair data. (Bombieri–Pila, Heath-Brown and Salberger collapse to
pair data on one circle: `overview/overview.md` row 16, `growth/auxiliary.md`, unrefereed.) The forced
l-adic valuation is the multi-level pigeonhole count of §3(c). It improves Theorem 1.1 by `O(M)` in
`W` only.
*Fake test*: fakes satisfying point residues at prime powers `≤ 2M` satisfy it (round 31 VIII.3).
*Char p*: the same Frobenius structure exists over `F_q(t)`, where the bound fails. *Verdict*: local
and dead for the rate.

**B2. Harmonic degree (CM theta series).** `Σ_{|z|^2=N} (z/|z|)^{4k}` is the N-th coefficient of the
theta series with character `λ^{4k}`. Within one N, raising the harmonic degree multiplies by `u_i^4`,
where `u_i = z_i/|z_i|`. The determinant `det(u_i^{4k})_{i, k<M}` is a Vandermonde in `u_i^4`, and
`u_i^4 − u_j^4` contains `Im(z_i z̄_j)`, so again we get pair data (Inferred). The S2 derivative `D`
(timeline item 2: `D^a θ(0)/Ω^a` divisible by `p^{a/p + a/p^2 + … − O(log a)}`) concerns Taylor
coefficients at a CM point. These are L-values (Hurwitz-type numbers), sums over all N, not one
Fourier coefficient. I see no way to make them see one circle (Heuristic).

**B3. Gram/positivity replacement for route 5's single value.** Replace `F_c(z)` by the Gram
determinant of the family `{F_{c^{(j)}}}`, or of point-evaluation functionals on polynomials of
degree `< D` under a Fock pairing. By Cauchy–Binet this Gram determinant equals
`Σ_{S, |S|=M} |det(z_i^s)_{s∈S}|^2 w_S` = `|V(z)|^2 Σ_λ |s_λ(z)|^2 w_λ` (Schur functions). Nonvanishing is
free, but the archimedean smallness is only `|V(z)|^2`, and the arithmetic is the cut identity (pair
data). The record caps the general form: "grid determinants with free nonvanishing are all
subcritical; a supercritical design needs an irreducible factor of degree ≥ 2, whose nonvanishing is
route 5's open zero estimate" (expert note §5 item 3, `round31/README_II.md` §1.3). *Verdict*: already
in the record and capped.

**B4. Bost slope method in place of ZE (Speculative).** One could try to replace route 5's ZE by
Bost's injectivity-by-Zariski-density: filter forms on `X_M` by order of vanishing along the
one-parameter subgroup `t ↦ z^t` through `e`. The germ's direction is `(log z_x)`, which is not
algebraic, so the jet maps have no integral structure. The method then becomes Baker/Schneider–Lang on
`G_m^d` with points `z^j`. The record treats this as round 26 `ze.md`: it is met by the trivial
subgroup for `|j| ≲ sqrt(n)/λ`, and single linear forms are faked (§5 item 4). *Verdict*: no new
input. The slope inequality reorganises the same Hilbert-function count that the zero estimate uses.

**What would have to exist for the S2 mechanism to transfer (Speculative programme).** A family of
rank `d ≫ M` of integral objects attached to the cluster, carrying an operator `y` with:
(i) `L∘y^m ≡ 0 mod ϖ^S` for the cluster evaluation `L`, with `S/m` large;
(ii) few elementary divisors above `ϖ^t` (`n_t ≪ d`);
(iii) `ϖ` running over many primes, so the gains add up across places.
On one circle the natural `L` has rank `M`, and its determinant is the Vandermonde. I could not
produce such a family. The record's theorem that the Gaussian-integer cocycle is the genuine problem
(round 31 XI.2) suggests there is no hidden high-rank structure (Heuristic).

### (c) Sharper Vandermonde divisibility than Theorem 1.1?

S2 has no such inequality for points. `LogVandermondeBounds` is about Weyl dimensions; `WindowKappa`
concerns geometric progressions `u^i` with `u ≡ 1 (mod l)` exactly, and cluster points are not a
geometric progression. The sharp classical fact behind `dimVN ∈ Z` is that a Vandermonde of `M`
integers is divisible by the superfactorial `∏_{k<M} k!`, which is optimal by consecutive integers.
For a subset of `Z[i]` the analogue is Bhargava's generalised factorial.

**Proposition (B-factorial on the circle; proof sketch, Inferred).** Let `l ≡ 3 (4)` with `l ∤ N`.
The circle `{|z|^2 = N}` meets each class of `T(Z/l^k)` (of size `(l+1)l^{k−1}`) in at most
equidistributed fashion in the worst case. Then for any `M` distinct points,
`v_l(∏_{i<j}(z_i − z_j)) ≥ Σ_{k≥1} (M^2/(2(l+1)l^{k−1}) − M/2)^+`.
The sum over `k ≥ 2` is `≤ M^2 · 1/(2(l+1)(l−1)) + O(M log_l M)`.
Weighted by `log l` and summed over `l ≤ M`, the `k ≥ 2` levels add `O(M^2)` to `log|Δ|`. Through the
cut identity this is `O(M)` in `W`. So the inequality `W ≥ 3M log M − O(M)` keeps its constant `3`. For
split `p ≤ M` the same holds with fibres of size `p−1` or `p` (Inferred, same computation).

*Fake test*: these are point residues at prime powers, which the window fakes satisfy. So even if the
`O(M)` term were improved it could not move the rate. **No inequality from this block raises the
constant 3, let alone the rate.**

## 4. Borrowable ideas, tests and first steps

| idea | analogue | fake test | char p | in record? | first step | P(rate) |
|---|---|---|---|---|---|---|
| Free nonvanishing via positivity (Gram) | Gram det of evaluation functionals / route-5 forms | fails only through integrality, but the design is subcritical before fakes matter | works over F_q(t) as well | yes, capped (§5 item 3) | Cauchy–Binet check that the Gram det equals `|V|^2 Σ|s_λ|^2 w_λ` (a 15-minute lemma; settles B3) | <0.5% |
| Multiplicity = ϖ-adic divisibility via an operator chain (`LatticeChainBound`, `TruncFinal`) | Frobenius `z^{l+1} ≡ N` at inert l | passes the residue axiom at prime powers, so faked | Frobenius exists in char p, so the uniform bound cannot follow | yes, as point residues / collisions | prove the B-factorial proposition of §3(c) with explicit constants and confirm the change is `O(M)` | ~0 |
| Averaging over levels (`TheoremM`) | average route 5's MIN^pair over a progression of orders `m'` or degrees | averaged inputs are grafted-faked (§5 item 2) | n/a | close to round 10 analytic / round 20 | none recommended | <0.5% |
| Slope inequality instead of ZE | Bost slope on `X_M` along `t ↦ z^t` | single linear forms in logs are faked (§5 item 4) | n/a | round 26 ze.md covers it | none recommended | <0.5% |

Overall: the block yields **no transferable mechanism for the growth rate**. The one conceptual lesson is
that S2 never proves a zero estimate. It makes the determinant positive by construction and puts all
the arithmetic into divisibility gained by operators that are ϖ-adically degenerate (Frobenius or
Hecke congruent to a derivation). Our problem has such operators (Frobenius at inert primes), but on a
rank-M evaluation they reproduce the Vandermonde, and hence pair or collision data, which fakes meet.
Something new would need a high-rank family attached to one circle with operators of small Smith
defect. I found none, and I give it ≤ 1%.
