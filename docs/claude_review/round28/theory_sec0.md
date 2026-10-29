# Round 28 — theory: multi-shift reach versus the zonotope order

Status: complete (unrefereed). Lane: theory. Date: 2026-10-29. Companion computation lane: `compute.md`.

## 0. Verdict

**No growth improvement and no uniform bound.** `M <= (2/3+eps) log R/loglog R` stands. The comparison of the
multi-shift reach with the zonotope order is **not decided at regime size**; what is proved is below.

1. **(a) Upper bounds.**
   * *No unconditional bound (proved, Thm 2.1, Cor 2.2).* Reach, `rho`, `zeta_pair`, `zeta_all` are all additive under
     products of boxes, so `(M8')^t` has `reach = 11t`, `zeta = 10t` (M8' input mod p, two primes). No inequality
     `reach <= zeta + c` holds on slab boxes. Product (block) profiles violate the character gap, so this does not help
     the route, but it shows that any such bound must use the gap, and that the favourable branch of the decision rule
     can be met vacuously.
   * *A conditional Liouville-capped bound (proved, Thm 2.3, Cor 2.4, Prop 2.5).* **Implied-slab obstruction:** if an
     integer form `phi` takes `R` values on `A` and at most `R - 1` on every `A_l`, then `phi^{R-1} in cap_l J(A_l)`; if
     also `phi^{R-1} notin J(A)`, the criterion stops at `R - 1`. For the **base coordinate** `sigma = -k_n` of a
     balanced sign-column profile, `|sigma| <= 2j` on the `T = 1` box (the averaged slab has the zonotope budget
     `Lambda(A) <= 2w`), so if its top level is rigid under pair moves (and `sigma^{2h} notin J(A)`), `reach(A, L_pair) <= 4j <= 2w < m' - 1`.
   * *But the hypothesis is a design coincidence.* Rigidity holds in M8 (the top level is the 7 Fano lines) and, in
     Sylvester boxes up to `n = 128` (IP), exactly when `h = 2L - 1`, i.e. for `L` a power of 2 (subspace configurations);
     it fails in every decided non-design case, including all decided random profiles. So (a) is **open** in the regime; the exact obstruction to a proof is
     that no other low-degree obstruction family is known (slab-power products cost `>= 2j log_2 n`; non-univariate
     residues exist in models but have no known structure).
2. **(b) Structure.**
   * *The syl8 `e = 6` obstruction is `sigma^6`, the sixth power of the base-point coordinate `k_8` (exact, Thm 4.1).*
     `|sigma| <= 3` because `sum_p s_p = 4 sigma` (Fano incidences); the top level `{1_l : l a Fano line}` is rigid under
     weight-2 moves and not under weight-4 moves, which is why all cube shifts rescue `(SHIFT_6)`.
   * *Why M8' reaches 11:* its base-coordinate top level is not rigid (`h = 4 < floor(7L/4) = 5`; 252 legal pair moves),
     because `L = 3` admits no design and its cuts are cube-bound vertex cuts; no implied-slab obstruction arises below
     `rho`. M8 and M8' differ in exactly the two features (design value of `L`; slab- versus cube-bound extremes) that
     pull in opposite directions at regime size, so neither is evidence about the regime.
   * *Most model failures are implied-slab powers* (6 of 8 defects fully spanned by `phi^e`, `phi` ternary).
   * *Chains and subsets (proved, Thm 3.1):* any chain or subset `A' subset A` gives the same conclusion, so
     `reach_chain >= max(zeta_pair, reach)`; on M8, `A cap {|sigma| <= 2}` lifts the reach from 6 to 7 `= zeta_pair`.
   * *First-order and linkage forms (proved; verified mod p):* `J(A_l) supset J(A) + R d_l J(A)` (equal in models up to
     one dimension), and `(SHIFT_e) <=> J(B)_{2d-e} = sum_l mu_l J(B_l)_{2d-e-deg mu_l}` for the tail `B = C \ A`.
3. **(c) Regime.** `(SHIFT)_A` (or its chain version) at regime size requires the absence of rigid implied slabs of
   range `<= 8j` **and** the absence of non-univariate obstructions in degrees `(4j, 8j)`. Hilbert-function tools
   (Macaulay, Gotzmann, Green, gin, linkage) control dimensions of single leading-form ideals, not the position of
   `~d^2` subspaces, and representation counting gives no obstruction for actual (symmetry-free) profiles. **No known
   asymptotic tool decides it.**
4. **Decision rule.** The theory does not decide it. It advises reading failures that are design effects (rigid
   implied slabs at Sylvester `L = 2^m`) and successes on product or vertex-cut boxes as uninformative about the
   regime, and using `zeta_pair` as the same-hypothesis competitor (§6).
