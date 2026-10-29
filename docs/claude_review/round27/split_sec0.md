## 0. Summary and verdict

**Verdict: `obstruction` (new, partly proved, partly conjectural with exact evidence), plus a new piece of the
exceptional set. No growth improvement, no uniform bound, no proof of split Question B or of Erdős #887.**
`M ≤ (2/3+ε) log R/loglog R` stands.

1. **New exceptional hypersurfaces (Theorem 3.1, proved, exact checks).** The second 2026 five-divisor construction
   (mnn31, Lean, `C = 7/2`: `n = u(u+2)(u+3)(u+5)`, `(u+2)(u+3) = 3(w+1)(w+2)`) is a rational curve with two conjugate
   real-quadratic branches of order 1, `B·C = 7`, hence **`α = 7/2 < 15/4`** (unbalanced: 7 of 15 cuts). It belongs to a
   3-parameter family of conics whose union is dense in the **quadric cone**
   `(u_2−u_1)(u_3−u_1) = (u_4−u_1)(u_5−u_1)`. So every admissible `Z` in split `A_5(δ)` (any `δ < 1/14`) contains these
   15 cones besides round 23's `W`. The cone is not defined over `Q` in the twisted `T_5`.
2. **Split profiles are 11-dimensional, and sections are profile-dependent (§1.3, Lemma 3.0).** The record transferred
   the balanced cap `15/4` to the split form only at the symmetric profile `x ≡ 1/2`. **Theorem 3.2 (proved):** for every
   auxiliary form with centrally symmetric support (all binomial and Vandermonde products, all known equality cases)
   the forced-divisibility cost is *independent* of the near-tight profile, so the cap `15/4` holds on the whole
   near-tight polytope `Π_nt` (exact identity + the refereed Theorem B). **Prop 3.3 (proved):** at the `W`-vertex (a
   "code" profile `1[|B ∩ R| even]`, `|R| = 3`) the cap holds for every form not divisible by the defining polynomial
   `g_W`, and holds in general iff `cost_{x_R}(g_W) ≥ 4·ord_e(g_W)` — one inequality, not decided (`deg W > 14`, exact).
   **Conjecture WTB** (cap at every near-tight split profile): 28 712 random sets, a local search, one-sided polytopes
   and a `k = 4` control, all with maximum ratio exactly 1.
3. **Bounded pair cofactors (Prop 3.4, Lemma 3.5, Theorem 3.6; proved).** Five divisors with `|t_ab| ≤ K` are exactly
   `S_K`-integral points of the affine variety `LM_5 ∖ ⋃Δ̄_ab` (fixed `S`, so Siegel/Corvaja–Zannier-type arguments
   are legitimate). But poles along the root arrangement, paid for by integrality, multiply only what plain sections
   already give: the Vandermonde has ratio exactly `15/4` on all of `Π_nt` (mediant lemma). So integrality is worthless
   wherever WTB holds, and decisive (no distance condition at all) wherever it fails. Each triple satisfies the exact
   identity `x^τ × y^τ = (t_bc, t_ca, t_ab)`, so bounded-residue triples are Möbius (`SL_2(Z)`/Farey) images of fixed
   rationals; this explains `A_3` false and gives no constraint at `k = 5`.
4. **Phase 0 (#887 proper; §3.6).** The forum's reduction (`z_j² − z_1² = 2k_j s + k_j²`, `k_j ≤ C² + 1`) is a family of
   complete intersections of quadrics: K3 for four divisors, general type from five on. The two known families are
   Pell conics on these surfaces (`k = (0,224,225,324,900)` and `(0,3,4,9,12)`, exact). Under Lang + Siegel, absolute
   #887 reduces to the absence of such curves for `K+1` divisors. Phase-0 only; not used for the torus problem.
5. **Transfer to `T_5` (§5).** Nothing ports as a proof. The twisted problem sits exactly at the symmetric profile,
   where the section cap is Theorem B; the split-only freedoms are profile asymmetry (worthless to symmetric forms,
   Theorem 3.2) and weight-1 curves (`W`, the cones), which the twist removes (parity lemma; Galois moves the cone).

**Status labels.** Theorem 3.1, Lemma 3.0, the cost identity, Theorem 3.2, Prop 3.3(a,c), Prop 3.4, Lemma 3.5,
Theorem 3.6: **proved** (unrefereed). Prop 3.3(b) decision, Conjecture WTB: **open** (evidence). §3.6: **derived**
(conditional reading only). §3.7: **heuristic**.
