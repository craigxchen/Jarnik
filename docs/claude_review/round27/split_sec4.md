## 4. Pre-registered tests T1–T5

**T1 (fake test; growth claims).** No growth claim is made: `M ≤ (2/3+ε) log R/loglog R` stands. The barrier
statements (Thm 3.2, Prop 3.3, Lemma 3.5, WTB) are statements about which auxiliary forms *fail*; they are consistent
with fakes, which satisfy every section inequality at the cube profile (record Cor C3). Theorem 3.1 concerns curves
over `Q` and is invisible to fakes. Nothing here is an input that fakes could violate, so T1 is vacuous and passed.

**T2 (`k = 4` negative control).** 
* *Barrier statements.* Their `k = 4` analogues hold and are consistent with `A_4(δ < 1/14)` being false: Theorem B
  with `n = 4` caps sections at `7/2`; the cost identity and Theorem 3.2 hold verbatim on the `k = 4` near-tight
  polytope (same proof); and the `k = 4` WTB was tested on 18 875 random sets with maximum ratio exactly 1
  (`k4_wtb.py`). A capped section method cannot exclude the dense `k = 4` families, as it must not.
* *The only exclusion mechanism in this write-up* is conditional: Lemma 3.5(i)/Prop 3.3(b) ("if WTB fails at `x`,
  then points with profile near `x` are degenerate"). At `k = 4` the dense families of gauss Thm 4.1 are twisted, hence
  at the symmetric profile, where WTB is Theorem B; so the conditional mechanism does not fire on them. Correct
  behaviour.
* *Phase 0.* The surface `S_k` is K3 for four divisors (density allowed: Erdős–Rosenfeld's infinite family) and of
  general type only from five divisors on; the conditional Lang reading (§3.6) therefore says nothing at `k = 4`.
  Correct behaviour.

**T3 (split test).** Nothing here proves split Question B, so absolute #887 is not claimed. Every statement is
compatible with the 2026 constructions: the m-czech curve (`α = 15/4`, profile `x_R`) lies on `W`, and Prop 3.3(a)
shows that no auxiliary form outside `g_W·Z[χ]` can exclude its points; the mnn31 curve (`α = 7/2`) lies on the cone
`𝒬` (exact), which Theorem 3.1 adds to every admissible exceptional set. Both land in the exceptional set; neither is
excluded by anything claimed.

**T4 (criticality test).** My results live inside the section class and say where it stops:
* *Lemma K1 / `τ_e = 10/3` / balanced `15/4`.* The record's balanced cap was proved at the symmetric profile only
  (Theorem B). For split points the profile is 11-dimensional and the section inequality depends on it (Lemma 3.0).
  So a section mechanism at an asymmetric split profile would **not** be one of the mapped caps. Theorem 3.2 shows
  that for all forms with centrally symmetric support (all binomial and Vandermonde products) the cap is nevertheless
  `15/4` at every near-tight profile; Prop 3.3 shows the same at the `W`-vertex for every form not divisible by `g_W`;
  WTB conjectures it in general. So this lane does not escape the cap; it locates the one unexamined corner (WTB at
  realisable one-sided vertices, decided by `g_W`).
* *Class `L` (round 24).* Lemma 3.5 adds integrality at the root arrangement (poles along `Δ̄_ab` paid by bounded
  residues) and shows it multiplies only what single sections already give.
* *Multi-point / Thue / Dyson.* Not used. The 4+1 decomposition (§3.7) is a heuristic count, not an argument.
* *Fixed-`S` Subspace.* Used only with `S = {∞} ∪ {p ≤ K}` genuinely fixed (Prop 3.4), where it reduces to
  Lemma 3.5; no exceptional set is taken uniform over varying prime sets.

**T5 (calibration).** Nothing excludes existing points: round-24's real twisted family is at the symmetric profile
(Theorem B, consistent); split Pell curves at `15/4` are protected by Prop 3.3(a); the cone points at `7/2` are put
into `Z`, not excluded; the record's 8- and 10-point families are untouched (no growth claim).

**Handoff constraints.** No uniform exceptional set from a fixed-`S` Subspace theorem (only `S = S_K` fixed with `K`
fixed); no Gaussian gcd or denominator cost discarded (tuples primitive, exact toric heights, `t_ab` exact); no drifting
phase set to zero (the phase-0 material of §3.6 is stated for #887 at `X = √n` only, where the phase is genuinely 0,
and is explicitly not used for the torus problem); no template-dependent threshold treated as uniform (in §3.6 the
per-template finiteness is legitimate for #887 because `n_0(C)` may depend on `C`, and nothing is claimed for the
torus problem from it).

## 5. Which steps use splitness; Gaussian analogue for `T_5`

| result | uses splitness? | twisted (`T_5`) analogue |
|---|---|---|
| Thm 3.1 (cones at `7/2`) | yes: the construction lives on the split torus; the cone is a split `Q`-subvariety | the real structure `u ↦ 1/ū` moves `𝒬` to a cubic; only `𝒬 ∩ σ𝒬` (a surface) could be twisted-rational. Second instance (after parity) of the twisted exceptional set being smaller |
| Lemma 3.0, Thm 3.2 | uses the 11-dim split profile space | twisted near-tight profiles are only `x ≡ 1/2`; there the statement is Theorem B (record). The split profile freedom is exactly what the twist removes, and Thm 3.2 shows it gives symmetric forms nothing |
| Prop 3.3 | yes: the `W`-vertex has weight-1 cuts | none: parity lemma (twisted curves have even weight on every cut) |
| Prop 3.4, Lemma 3.5 | no | identical over `Z[i]`: bounded residues = integrality at the Gaussian root arrangement; the Vandermonde has ratio `15/4` at `x = 1/2`, so Δ-poles add nothing |
| Thm 3.6 (triples) | integer factorisation of `α_ab` | over `Z[i]` the same identity holds with `x^τ, y^τ ∈ Z[i]^3`, and conjugation relates `y^τ` to `x̄^τ` (cut sides swap), giving `GL_2(Z[i])`-structure with a reality condition; not pursued |
| §3.6 phase 0 | #887 only | the near-axis bound (record lit18 §2.3) |

**Decision.** Nothing in this lane ports to a proof of Question B for `T_5`. What ports is negative and precise: the
two split-only freedoms (profile asymmetry; integrality of residues) are worthless to section methods for symmetric
supports, and the twisted problem sits exactly at the symmetric profile where the cap is already a theorem. The only
split-only positive possibility left by this analysis — WTB failing at a one-sided realisable vertex via `g_W` —
cannot exist in `T_5` (no weight-1 curves) and would at best prove a restricted *split* statement.
