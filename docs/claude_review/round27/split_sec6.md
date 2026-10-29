## 6. Scripts and certificates

All under `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/`; runs under
`prlimit --as=3500000000`. The repository was not edited.

| file | checks | label |
|---|---|---|
| `dl/Erdos887.lean` | mnn31 Lean file (downloaded as data, read, not executed) | source |
| `czech72.py` | `C = 7/2` family along 40 Pell points: 5 divisors, `C_act → 2√3`, `γ → 7/2`, `log(dist^{7/2}H)` periodic | exact integers |
| `fam72.py` | second factorisation identity; exact Jacobian rank 3 of the 4-parameter family | symbolic / exact |
| `cone_checks.py` | (a) cone identity at 6 Pell points; (b) m-czech points on none of the 15 cones; (c) triple identity `x×y = t` on 201 tuples; (d) Galois image of the cone | exact |
| `phase0_k.py` | phase-0 `k`-tuples and the identity `z_j² − z_1² = 2k_js + k_j²` on both families | exact |
| `profiles.py` | `dim Π_nt = 11`, 787 vertices (278 one-sided), binomial cost `≥ 4` at all vertices (`|m_a| ≤ 3`) | exact rank; LP vertices (float) |
| `reglib.py` | regularity index mod `2^31−1` (checked: `reg(nR_5) = ⌊10n/3⌋`, `n ≤ 4`) | exact mod `p` |
| `search_beat.py` → `sb_rand.out`, `sb_lop.out` | 28 712 sets: `max 4reg/min_{Π_nt} cost = 1` | mod-`p` reg, float LP |
| `local_search.py` → `ls7.out`, `ls8.out` | local search from zonotopes: `max(4reg − mincost) = 0` | same |
| `onesided_vol.py`, `qf_reg.py`, `qf_reg_big.py` → `qf_reg_big.out` | one-sided free polytopes: volumes, `reg(nQ_F)/n ≤ 3.5` (`n ≤ 5`) | exact mod `p`; float volumes |
| `cone_section.py` | quadric-cone form: min cost over `Π_nt` = `25/3 > 8` | mod-`p`, LP |
| `k4_wtb.py` → `k4_wtb.out` | `k = 4` control: 18 875 sets, max ratio 1 | mod-`p`, LP |
| `w_equation.py` | no form of degree `≤ 14` vanishes on `W` (3100 points over `F_{1000003}`) | exact mod `p` |
| `tangent_W.py` | branch directions of `W` at `e` span the irreducible octic `Σ_W` | symbolic |
| `../split_notes.md` | working notes N1–N6 | — |

## 7. Open questions produced by this lane

* **Q-split-1 (decisive for the section method on split points).** For the round-23 hypersurface `W` with defining
  polynomial `g_W` and code profile `x_R`: is `Σ_{B : |B∩R| even} π_B(g_W) ≥ 4·ord_e(g_W)`? (`ord_e(g_W) ≥ 8`; if the
  order-2 branches give two sheets, `16`.) Equivalently: is `(W̄·C) ≥ 4·mult_e(W)` for a torus translate `C` of a
  `W`-curve? *No* gives a restricted theorem (near-tight split points with profile near `x_R` lie on `W` plus a fixed
  closed set); *yes* (with WTB) closes the section method on split points. It is a finite computation (tropical
  implicitisation of `W`).
* **Q-split-2 (WTB).** Prove `4·reg(S) ≤ Σ_B x_B π_B(S)` on `Π_nt`; it needs an asymmetric slicing lemma (§3.4).
* **Q-split-3.** Is `K_{M̄_{0,7}} + D_int` big (log general type of `LM_5 ∖ ⋃Δ̄_ab`)? If not, Lang–Vojta does not even
  predict degeneracy of bounded-residue points without the A-point condition.
* **Q-split-4 (phase 0).** Does any complete intersection `S_k` (six entries) contain a rational curve with two conjugate
  real-quadratic points at infinity (a six-divisor Pell family)?
