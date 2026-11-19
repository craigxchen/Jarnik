# Referee/analyst report, block "arith": can any proved arithmetic tool replace BMY in Part IV?

Scripts and outputs are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref4_arith/`. All of the code is my own.
* `trace_indep.py`: my own GF(2^d) arithmetic. It rebuilds the trace construction for q = 8, ..., 128 and evaluates (IV.2) on it.
* `kappa.py`, `kappa_opt.py`: Astra's blow-up accounting redone with a general constant κ in place of 3, plus a search over the weights.
* `transfer.py`: the record's pigeonhole/chain-lemma bound transferred to integral fakes over Z, F_2[t], F_3[t] and F_{m+1}[t].

## 0. Verdict

**No proved arithmetic tool replaces BMY, not even in the growth regime or with bad constants.** I found no proof of anything new for integer parents. Astra's sentence "I have no arithmetic substitute" is correct, and so far it is unavoidable. What this block does establish:

1. **(IV.2) is false in characteristic 2.** This is checked exactly. The trace construction satisfies every hypothesis of (IV.2) with δ = E_ij = 0, yet the left side is positive for q ≥ 32: 341/2 at q = 32, 3717/2 at q = 64, 19685/2 at q = 128. So the proof must use characteristic 0. It does: Langer's inequality is a char-0 theorem. This is a consistency check that the proof passes.

2. **The strength of BMY is irrelevant; only its existence matters** (new, computed exactly). Replace 3 by any κ, i.e. assume (K_Y+B)^2 ≤ κ e_orb(Y,B) on this class of lc rational pairs. Then Astra's accounting gives **m ≤ 2κ(κ+1)**: 24, 40, 60, 144 for κ = 3, 4, 5, 8. For κ > 3 the tangency-chain weight must be changed from b = a to some b in (a, 1). With reduced boundary (a = 1, the classical log Miyaoka–Yau), κ = 3 gives m ≤ 32.
   * Consequence: over \bar F_2(t) no inequality c_1^2 ≤ κ c_2 holds for any κ on these pairs.
   * Consequence: an arithmetic orbifold BMY with *any* constants, and an error O(W), would bound m uniformly. That is the full conjecture for integral fakes. Parshin and Moret-Bailly showed that arithmetic BMY implies Szpiro and effective Mordell.

3. **No Bogomolov-only proof exists in the form asked** (a no-go argument; I proved no bound).
   * Y is a rational surface.
   * By Shepherd-Barron (Invent. Math. 106 (1991)), Bogomolov's inequality for rank-2 semistable sheaves holds in characteristic p on surfaces that are neither of general type nor quasi-elliptic with κ = 1. I confirmed this through two secondary sources only.
   * Hence any argument of the following shape is refuted by the trace configuration over \bar F_2: build a rank-2 bundle from the configuration, prove its semistability from the incidence and numerical data, compute Δ < 0.
   * The only escape is a characteristic-0 proof of semistability. In practice that means Miyaoka's generic semipositivity of the (log) cotangent sheaf, i.e. a derivative along the base.
   * The arithmetic Bogomolov inequality (Miyaoka; Moriwaki 1993–95; Soulé 1994) is killed for the same reason. Its geometric shadow on a P^1-fibration only needs generic-fibre semistability. On P^1_K that is strong semistability, which does not depend on the characteristic.

4. **Hodge index is killed by F3**, and the record already capped it (round 20): Faltings–Hriljac, arithmetic Riemann–Roch, Minkowski/Siegel, Zhang minima and Arakelov adjunction all have char-p analogues that the trace surface satisfies.

5. **F3, stated precisely and tested.**
   * The record's proved bound transfers to F_q[t]. Its gain lives entirely on residue fields of size < m, so it is **vacuous at q = m+1 and consistent with the trace construction** (no contradiction).
   * The trace construction sits at W = m log(m+1), which is below the transferred integer bound (≈ 2 m log m for no-norm fakes with any parent, 3 m log m for split parents). So the trace table could never live at its own weight over Z.
   * F3 therefore kills exactly the arguments that do not exploit small residue fields or characteristic 0.
   * The sharp open question is **F3′**: are there growth-scale integral fakes over F_p[t] with p fixed? I found no construction. The trace construction does not descend to F_2[t]: the Frobenius-stable subspaces form a divisor lattice, not a projective geometry. That remark is heuristic.

**Significance for the goal: about 1/10.** Item 2 sharpens what is missing: a *qualitative* arithmetic BMY. Items 3–5 close candidate routes cheaply. Nothing changes the growth rate, the uniform bound, or Questions A and B.

## 1. What the C(t) proof really uses

### 1.1 Characteristic dependence (exact)

`trace_indep.py` builds f_a = (t^{q/2} − L_a)/t and P_a = (f_a, 1) over F_q, q = 2^d, d = 3..7. It checks, for every pair, the three facts below; there are no failures over all 8001 pairs at q = 128.
* D_ab vanishes on (H_a ∩ H_b) \ {0}.
* deg D_ab equals the number of common marks, so the cofactor is constant.
* deg D_ab − critical = −1/4, so δ = 0.

Layers: L = q − 1, every s_h = q/2 − 1, and E_ij = 0. The (IV.2) left side (3/4)(m−24)L + Σ(s_h − m/2)^2 is:

| q | m | (IV.2) left side |
|---|---|---|
| 8 | 7 | −175/2 |
| 16 | 15 | −195/2 |
| 32 | 31 | +341/2 |
| 64 | 63 | +3717/2 |
| 128 | 127 | +19685/2 |

The right side is 0. So (IV.2) holds for q ≤ 16 and fails for every q ≥ 32 in characteristic 2.

This is **not** a counterexample to Part IV, which is over C[t]. It shows that the inequality has char-0 content, so any "arithmetic substitute" must have it too. That is the starting point of filter F3.

### 1.2 Any finite κ suffices (exact; new)

I re-derived Astra's local accounting with κ e_orb − (K+B)^2 ≥ 0 in place of κ = 3 (`kappa.py`, sympy). The four half-fibres cancel exactly for every κ: the remainder is (2a − a^2)S, independent of κ. A centre of k ≥ 3 sections gives F(k,0,0) = κ + 1 − 2ak + a^2 k. A transverse unit of intersection costs (κ−2)a^2. Writing I = (m−1)S/2 and S ≤ Σ(s_h − m/4) (IV.3 with δ = 0), the per-layer contribution is identically

  κ + 1 + a(κa − 4) m/8 − (κ−2) a^2 (s_h − m/2)^2 / 2 .

At a = 2/κ the m-coefficient is −1/(2κ), so δ = E = 0 forces **m ≤ 2κ(κ+1)**.

Tangency chains need a separate check:
* For κ = 3 and b = a = 2/3 the chain inequality holds; this is Astra's case.
* For κ > 3, b = a fails: F(2,a,a) exceeds the unit charge.
* A chain weight b ∈ (a, 1) restores it. `kappa_opt.py` finds (a, b) = (1/2, 143/200) for κ = 4, (2/5, 31/40) for κ = 5, and (1/4, 41/50) for κ = 8. The bounds are 40, 60 and 144, exactly 2κ(κ+1).
* Effectivity of K+B is unaffected, because the extra vertical weights are ≥ 0.
* Reduced boundary (a = b = 1, κ = 3): chains hold and m ≤ 32.

Two consequences:
* **Interpretation.** The weight a = 2/3 is the orbifold weight 1 − 1/3 of a cyclic triple cover branched along the sections, and the half-fibres are double covers. The C(t) proof is Hirzebruch's Kummer-cover use of Miyaoka–Yau. Over \bar F_2, (i) and (ii) below are incompatible with the trace family for every κ. That is the classical failure of Hirzebruch-type inequalities for F_q-rational arrangements (all lines of P^2(F_q) have only (q+1)-fold points).
  * (i) Any linear Chern inequality, with any κ, gives a uniform bound.
  * (ii) The trace family has unbounded m.
* **Arithmetic reading.** An arithmetic orbifold BMY for P^1_Z with m sections and an error A·W would give m ≤ 2κ(κ+1+A) for integral fakes. That is a *uniform* bound, much more than growth. Literature status:
  * Arithmetic BMY (Parshin 1988; Moret-Bailly) is a conjecture, known to imply Szpiro's discriminant conjecture and effective Mordell. I checked this through secondary sources (mathnet survey, Moret-Bailly's Astérisque 183 exposition); I did not read the originals.
  * I found no formulation for P^1 with sections; it would presumably be stated through the triple covers above.
  * So **the arithmetic replacement for BMY is at least Szpiro-strength**, even with terrible constants.
  * Caution: the record's barrier 6 ("abc on single relations is vacuous on squarefree profiles") does *not* apply here. BMY uses all m sections at once, which is why it is stronger than abc on Ptolemy relations.

## 2. Item (1): Bogomolov, arithmetic Bogomolov, Reider

**Statement attempted.** Find an explicit rank-2 bundle E on Y (or on P^1_Z blown up), not Ω_Y(log B), that is (a) semistable and (b) has Δ(E) = 4c_2 − c_1^2 < 0 whenever m > M_0. Bogomolov would then bound m.

**No-go argument (proved, modulo the cited statement).**
* Y is a blown-up Hirzebruch surface, so it is rational.
* Shepherd-Barron (1991): on surfaces that are neither of general type nor quasi-elliptic with κ = 1, Bogomolov's inequality holds in characteristic p for rank-2 semistable sheaves. Confirmed via arXiv 1907.08378 and 2008.09800, both quoting it; I did not read the original.
* Suppose E is built from the configuration (elementary transformations along the C_i, a Serre construction on the contact cycle Z, kernels of evaluation maps), and (a) and (b) are proved from the hypotheses of (IV.2) by intersection and incidence reasoning.
* Then the same construction over \bar F_2 applied to the trace configuration (§1.1: all hypotheses hold, m = q − 1 unbounded) gives a semistable rank-2 sheaf with Δ < 0 on a rational surface in characteristic 2. That contradicts Shepherd-Barron.
* Hence **no such proof exists**.

**Escape routes.**
* The semistability proof uses characteristic 0. On rational surfaces the only standard char-0 semistability sources are Miyaoka's generic semipositivity (applied to Ω^1(log)) and Kähler–Einstein/uniformization. Both concern the (log) cotangent sheaf, or bundles built from derivatives such as principal-parts bundles. These are the bundles the task excludes.
* An h^0 jump in characteristic 2 on the specific trace configuration could invalidate a semistability check done there. But a proof must deduce semistability from the hypotheses, and those are the same over \bar F_2.

**Reider-type arguments.** A Serre bundle 0 → O → E → I_Z(D) → 0 on the contact cycle has Δ = 4 deg Z − D^2. When Δ < 0, Bogomolov instability produces a curve through much of Z with small D-degree. In our setting the sections themselves are such curves, so the conclusion is satisfied by the configuration and gives no contradiction. Reider also holds on rational surfaces in characteristic p (same Shepherd-Barron paper; secondary sources only), so it is F3-killed in any case.

**Arithmetic Bogomolov.**
* *Statement.* The arithmetic Bogomolov–Gieseker inequality: if X is a regular arithmetic surface and (E,h) a Hermitian bundle whose generic fibre is semistable on the geometric generic fibre, then ĉ_2 − (r−1)/(2r) ĉ_1^2 ≥ 0. Moriwaki proved it (arXiv alg-geom/9305005; higher-dimensional version alg-geom/9307004; Amer. J. Math. 1995 venue not confirmed). It is credited also to Miyaoka. Soulé, "A vanishing theorem on arithmetic surfaces", Invent. Math. 116 (1994) 577–600, is the related vanishing input. All of this is from abstracts and secondary sources (Naumann math/0508641; Soulé 2005).
* *Why it fails here.*
  * For X_Q = P^1_Q, semistable means O(a)^2. That is strongly semistable, so the geometric shadow is the *relative* Bogomolov inequality for P^1-fibrations, which has no characteristic-0 content. In char 0 it follows from Bogomolov instability, because the fibre class is nef. In char p, on our rational Y, I expect it from Shepherd-Barron's instability form or from Langer 2004; I did not verify the exact form.
  * Any bound obtained this way would therefore also hold over F_q(t), and the trace construction refutes it. So it is **F3-killed** (with the stated caveat).
* *Where the derivative was supposed to enter.* Soulé (2005) recalls Miyaoka's observation that a general extension of ω_{X/S} by π^*Ω_S is semistable. Over Spec Z there is no Ω_S. That extension class is the Kodaira–Spencer class of the marked family, i.e. the derivative. In our setting it is the extension 0 → π^*Ω_B(log) → Ω_Y(log C) → Ω_{Y/B}(log C) → 0, whose class records the velocities d/dt of the sections. This is exactly the datum the integer configuration lacks.
* **Verdict: GAP, not repairable with known tools.**
* **F1:** not reached (no inequality is produced).

## 3. Item (2): Hodge index and Noether-type identities

**Statement.** Use only the arithmetic Hodge index theorem (Faltings 1984, Hriljac 1985), the arithmetic Noether formula, adjunction and Riemann–Roch on (blow-ups of) P^1_Z with the m sections.

**Filter F3.** Each of these has a char-p geometric analogue: Hodge index holds on every smooth projective surface, and Noether and RR are identities. The trace surface over F_q satisfies all of them. So any inequality derived from them alone holds for configurations with m = q − 1, and none bounds m. The same applies to Minkowski/Siegel lemmas, Zhang's successive minima, arithmetic Hilbert–Samuel and Bost slopes: each has an F_q(t) analogue that the trace data obey. Noether *inequalities* (K^2 ≥ 2p_g − 4 and the like) point the wrong way.

**Record.** `round20/arakelov.md` §§2–3 already shows that on genus-0 arithmetic surfaces Hodge index is an *equality*, namely the cut identity, so the class is pair data and capped. That was for genuine clusters on the conic model; the same holds for integral fakes on P^1_Z. F3 is an independent, cheaper second kill.

**Verdict: killed (CORRECT that it cannot work).**

## 4. Item (3): substitutes for the derivative

* **p-derivations / Fermat quotients.** The record already found these fake-blind (round 29 §2.1 item 7). They also act one prime at a time. BMY's derivative acts at every layer simultaneously through one global section dt of Ω_B(log). A p-derivation at p sees only the layers at p. I see no way to assemble them into a global sheaf with a Chern-class inequality. **Not a substitute.**
* **Archimedean derivative.** This is the only char-0 structure Z has. It acts at one place, the cone at infinity, while all contacts are finite. Any link to the finite places goes through the product formula. Round 20 Theorem C says every height relation obtainable from local data and product formulas is a combination of pair product formulas: pair data, satisfied by fakes. It passes F3 formally, but **F1/pair-capped**.
* **Kodaira–Spencer over Z.** It would require a "derivative in the p-direction", the classical missing ingredient of abc. **Open, at least abc-hard** (§1.2).

## 5. Item (4): filter F3 and the record's proved bound

### 5.1 F3, stated precisely

Let 𝒜 be an argument that is valid for integral fakes (IV.1) over every A ∈ {Z} ∪ {F_q[t]}. Norms are |ℓ| = q^{deg ℓ}, W = (deg N) log q, and A's residue fields and units are used only through their sizes. If 𝒜 yields W ≥ Φ_A(m), then Φ_{F_q[t]}(q−1) ≤ (q−1) log q for all q = 2^d, d ≥ 3.

In particular, no argument of this kind can give W ≥ ω(m) m log m − ε(m) with an error ε that is o(m log q) at q = m+1. A growth improvement over Z must therefore use one of:
* (i) residue fields of size < m (small primes), or
* (ii) the characteristic 0 of Q: archimedean place, derivative, ordering, or global addition.

### 5.2 Testing the record's proved bound

**Transfer (derived here; numbers are evidence).** For integral fakes the exact upper side is

  Σ_{i<j} log|c_ij| + Q/2 ≤ mW/8 + K log C1, with c_ij = D_ij/g_ij.

I re-derived it from (IV.1) and the layer identity: (m−1)s/2 − s(s−1)/2 − m(m−1)/8 = m/8 − (s − m/2)^2/2. It holds verbatim over F_q[t].

The lower side transfers too:
* Distinct prime divisors multiply in any PID.
* At ℓ ∤ N, collisions are counted in P^1(A/ℓ^a), with (|ℓ|+1)|ℓ|^{a−1} classes.
* At p | N, the chain lemma applies with λ = 1/(|p|+1). Rows at the same level t ≥ 1 are forced into one class mod p^t and then spread over |p|^a classes.

`transfer.py` gives the coefficient c in W ≥ c·m log m:

| m | Z, any N | Z, split N | F_2[t] | F_3[t] | F_{m+1}[t] |
|---|---|---|---|---|---|
| 256 | 0.78 | 1.89 | 0.75 | 0.81 | 0 |
| 1024 | 0.98 | 2.10 | 0.96 | 1.01 | 0 |
| 4096 | 1.14 | 2.25 | 1.12 | 1.13 | 0 |

The values increase slowly toward 2 (any N, including F_p[t] with p fixed) and 3 (split N; Astra's "coefficient three"). The leading constants are by Mertens; the error terms are O(m log q) on the F_q side.

**Against the trace construction.** At q = m+1 every residue field has size ≥ m, so the gain is exactly 0. The transferred inequality on the trace data reads Q/2 = L/8 ≤ mL/8 (checked for q up to 2^20). **Consistent, no contradiction**: the proof is not uniform in q, and its only content is small residue fields.

**Integer comparison.** The trace weight is W = m log(m+1), below the integer bound (≥ 2 m log m − O(m) in the no-norm model with any parent). So no integer configuration can mimic it at that weight; this is independent of Astra's symmetry argument.

**What F3 actually filters.** F3 does *not* test the pigeonhole approach. It tests the "geometric" tools: Hodge index, Bogomolov on rational surfaces, Arakelov positivity, geometry of numbers. The relevant open test for small-prime methods is:

**F3′ (open, proposed).** Is there, for fixed p (say p = 2), a family of integral fakes over F_p[t] with deg N = O(m log m / log p) and m → ∞?
* If yes: Z's small primes cannot be the whole story, and growth needs char-0 input, i.e. an arithmetic BMY-type statement.
* If no: the proof of "no" over F_p[t] is a candidate template for Z.

My attempt to descend the trace construction to F_2[t] fails structurally. Columns with F_2[t] entries would need Frobenius-stable F_2-subspaces of F_{2^d}. These correspond to divisors of x^d − 1 (kernels of g(σ)), so they form a distributive lattice with no family of "hyperplanes" pairwise meeting in codimension one. This is heuristic and does not settle F3′.

## 6. Item (5): other ideas, with filter status

| idea | precise statement attempted | F1 | F3 | status |
|---|---|---|---|---|
| any-κ BMY (§1.2) | (K+B)^2 ≤ κ e_orb ⇒ m ≤ 2κ(κ+1) | n/a (C(t) theorem) | correctly char-0 | **proved (C(t))**; arithmetic version ⇒ Szpiro-type; open |
| triple-cover Szpiro | Szpiro-type inequality for y^3 = ∏_i [P_i,(x,z)] (genus ≈ m; discriminant ∏D_ij^2 ≈ (m^2/4)W against conductor ≈ W) | unknown | function-field Szpiro (Szpiro, Esnault–Viehweg, Tan) is char-0 or separable-only, so it passes | **speculative**; whether the genus-dependence of the constant (≥ m^2?) leaves anything is unchecked; conditional either way |
| design-matrix rank bounds (Barak–Dvir–Wigderson–Yehudayoff 2011) | quantitative Sylvester–Gallai on the layer incidences | — | valid in char 0 or large char, so passes | **not set up**: over Z the incidences live in different residue fields F_p, with no common ambient space |
| Mason/Wronskian on quartets | Plücker 3-term relations among D_ij | vacuous on squarefree profiles (record barrier 6) | Mason fails in char p, passes | **vacuous**: with f_i = (f_i,1) normalisation, coprime parts have degree ≈ L/8 against radical ≈ 3L/8 |
| S-unit counts for the split unit equations λ_iA_i − λ_jA_j = λ_{ij}A_{ij} | — | count depends on S (record barrier 5) | — | **killed** |

## 7. Novelty and significance

**New against the record.** These have no prior entry for Parshin, Shepherd-Barron, Soulé or Raynaud in the record:
* the explicit char-2 failure of (IV.2);
* the "any κ" generalisation with m ≤ 2κ(κ+1);
* the Shepherd-Barron no-go for Bogomolov proofs;
* the precise F3 and F3′ and the transfer table;
* the identification "the arithmetic replacement for BMY is at least Szpiro-strength".

`round4/L7.md` only noted that log BMY has slack at seven marks. The Hodge-index kill duplicates round 20.

**Significance.** Low. The block shows that the gap between (IV.2) and the integers is the classical gap between function-field and arithmetic BMY (Szpiro/abc territory), and that no proved Arakelov tool crosses it. The one useful sharpening for planning is §1.2: the missing input is qualitative (any κ, any O(W) error), so a future route needs *some* arithmetic Miyaoka–Yau-type inequality for P^1_Z with sections, not a sharp one. That is still a major open conjecture.
