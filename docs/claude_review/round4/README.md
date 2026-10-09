# Round 4: flipped Hadamard cores, (L_7), hardness, and methods outside the barred classes

Each write-up has an adversarial referee report; the statements below are the referees'
corrected versions. Nothing here proves the uniform theorem or a general growth improvement.

* `flipped.md` (fully flipped Hadamard phase lemma). Exact reduction: for any normalised Hadamard
  core with arbitrary flips, each label pinning is a nonzero linear form in two logarithms,
  `M Log(G_a/conj G_a) - 2 Log(P_a/conj P_a) + b Log i`, with error at most `M Delta`, where `G_a` is
  the label's block and `P_a` the Hadamard-pattern product of the flipped primes. **Conditional
  Theorem L:** a Lang-Waldschmidt-type two-logarithm bound over `Q(i)` with naive heights implies
  the lemma whenever the flipped mass is at most `(1/2 - delta) W` (this covers the Prop-P barrier
  class). Proved Baker-type bounds (product of heights) do not suffice. The referee shows the
  heavy-flip regime is vacuous for the canonical Paley assignment (`W_F <= W/2 + O(M log C)`), so
  the residual case is the critical band `W_F/W -> 1/2`. New integer prime-field uncertainty lemma
  for Paley matrices: `|supp Hn| >= (q+1)/2 - |supp n|`.
* `L7.md` ((L_7) undecided). With the anchor and two labels at `inf, 0, 1`, residues are `+-1` and
  determinants are exact products of contact integers, so (L_7) is a statement about boundary
  heights alone (residue-free formulation; the normalisation is from `growth/ptolemy.md`). A
  Zariski-dense set of `Q`-rational curves of class approximately proportional to the balanced class
  would refute (L_7). The referee found numerically that the class of twisted cubics meeting 13 of
  the 20 Kapranov planes is carried by a covering 3-dimensional family, which combined with
  smoothing of free curves would likely give such a refutation; this is not proved.
* `hardness.md` (search result: no implication from the uniform theorem to a named open problem).
  Proved: the uniform theorem holds iff for every `kappa` some short-arc tuple set `A_M(kappa)` is
  finite (equivalently not Zariski dense); the absolute form is Roth-type finiteness at Dirichlet
  exponent `1/2` for the diagonal line in a fixed `X_M` (known families force `M >= 9`); additive
  energy of squares in two intervals of length `c sqrt X` is elementary
  (`E <= (2 + c^2)|A||B|`), so Rudin-type consequences through energy are already theorems; the
  circle's torus is anisotropic, so there is no equivariant transfer to the Erdos-Rosenfeld or Ruzsa
  divisor problems.
* `outside.md` (GAP/Fourier, energy, incidences, equidistribution). **Theorem A.1:** for large
  primes, all monomial gaps carry exactly the information of the rational characters imposed
  simultaneously. **Theorem A.3 (fake circles):** angle models with genuine primes, a genuine flipped
  Hadamard exponent profile and real angles carry `M`-point clusters on arcs `C sqrt R` with
  `log R <= (15+o(1)) M^2 log M` while satisfying every monomial gap, residue-strengthened character
  gap and exponent identity. So these inputs cannot prove `M = o((log R/loglog R)^(1/2))`. The
  referee stresses that nothing proved excludes an `o(log R/loglog R)` growth improvement from
  simultaneous rational-character admissibility: the decisive open question is whether fake circles
  with `log R = O(M log M)` exist (e.g. fully flipped Paley cores with boundedly many copies).
