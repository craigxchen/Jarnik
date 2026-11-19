# Referee report, round 31, Part III, block "proj" (astra_free_III §§III.7, III.8, III.10)

Scope: `astra_free_III/answer.md` §III.7 (global projective reconstruction, robustness to deletions, the
"almost every parent in a pool" strengthening, the planted fourteen-prime example), §III.8 (five-mark
telescoping (III.6), five-mark consistency, gluing on regular quartets), §III.10 (the residual
"one scale per mark" problem, integral realization, the primitive-column criterion, the quadratic form).
I wrote my own scripts. They are in `round31/ref3_proj/`:

- `gauss.py`: Gaussian arithmetic. From a list of lattice points it builds the star matrix (gcd removal,
  factorization, exponent array `a_ip`, absorption of units by `1, 1+i, i, 1-i`, then `g_ij`, `d_ij`, `t_ij`). It asserts
  `g_ij | [P_i,P_j]` and `t_ij != 0`.
- `recon.py`: the III.1/III.7 pipeline. For each quartet it checks Plücker and the lcm formula for `(A,B,C)` from
  the `d_ij` alone. It certifies local uniqueness (some modulus `> 2 U U'`) and recovers the small kernel ray from
  `(A,B,C)` and the bounds only, by exhaustive enumeration of a reduced 2-D lattice. It then compares with the true
  ray and reconstructs the directions in a chart (anchors to `∞,0,1`), with iterative propagation.
  `out_apriori.txt` uses bounds computed from `d_ij` and the arc width only. `out_oracle.txt` uses the exact `|t_ij|`
  as bounds, which is a best case.
- `algebra_checks.py`: random prime-power arrays. It checks the telescoping and the unit/valuation claims of
  III.8, the pattern-indicator expansions used in III.7/III.8, and the III.10 saturation, realization and
  primitive-column criterion (against brute force over multipliers).
- `frame_check.py`: the frame bound `λ_max(G) ≤ 2` for sign Gram matrices with negative off-diagonal.

## Verdicts

| Claim | Verdict |
|---|---|
| III.7 local uniqueness (one projective `[x:y:z]` when `B > 2U^2`), cross ratio determined by the `d_ij` | **CORRECT** (proof checked; numerically confirmed, no mismatch on 0+4+8·5 a-priori-certified and 152 oracle-certified quartets) |
| III.7 main theorem: for `m` large, parent + labelled `d_ij` determine all star directions up to one `PGL_2(Q)` | **CORRECT**. "Sufficiently large" is explicit and absolute: `m ≥ 16,390` with Astra's constants (`τ=1/32`), about `10^4` with sharp counting. |
| III.7 robustness (≤ `m/10` deleted partners per mark; `o(m^2)` deletions → `m−o(m)` marks) | **CORRECT**, needs `m(1−6η) ≳ 1.64·10^4`, i.e. `m ≳ 4.1·10^4` at `η=1/10` |
| III.7 "almost every parent in a pool of `K ≥ m^{1+η}` primes excludes every assignment on every Paley quartet" | **CORRECT** (union bound checked). Peripheral; Paley-profile-specific. |
| III.7 planted fourteen-prime parent (one passing quartet, no five marks pass) | Not recomputed by me. It was **reproduced independently** in `ref3_quartet.md`; I rely on that report. |
| III.8 (III.6): three independently solved quartets on anchor pair `ij` telescope when `Q > 2U_1U_2U_3` | **CORRECT** (proof checked line by line; 0 failures of the algebraic ingredients in 20,000 random prime-power arrays) |
| III.8 five-mark consequence on regular five-sets (`τ=1/64`) | **CORRECT** |
| III.8 gluing on all regular quartets of a large profile | **CORRECT, sketchy write-up**. The finite-set step is better done by the induction given below. Explicit size needed is about `m ≳ 5·10^5`. |
| III.10 `t_ij = c s_i s_j h_ij` parametrization | **CORRECT** (trivial) |
| III.10 integer two-row realization via saturation; all realizations are `M·B`, `det M = d` | **CORRECT** (checked on random cases) |
| III.10 primitive-column criterion (columns of `B` primitive, and at each `p | d` the reductions miss a point of `P^1(F_p)`) | **CORRECT** (proof checked; 47 random cases agree with brute force) |
| III.10 exact norms ⇒ `S = M^tM`, determined by three directions | **CORRECT** |
| Does any of this give a bound? | **No.** These are rigidity/uniqueness statements conditional on existence. They change neither the growth rate, the constant 2/3, the uniform bound, nor the status of Questions A, B or route 5. |

## III.7 — global projective reconstruction

**Local uniqueness.** Take two solutions of `Ax − By + Cz = 0` with `|x|,|z| ≤ U`. Since `A, B, C` are pairwise coprime,
`A(xz' − x'z) ≡ 0 (mod B)`, so `B | xz' − x'z`. The left side has absolute value at most `2U^2 < B`, so `x:z = x':z'`, and `y` follows.
This needs `x, z ≠ 0`, which holds because minors are nonzero for distinct directions. It also needs pairwise coprimality,
which comes from the `h_3 − h_2` rule of III.1. That rule is correct, and the lcm formula `A_k = sqrt(L/D_k)` holds at each prime because at most one
`A_k` has positive valuation there. So the raw products of the six minors are fixed projectively. Hence every cross
ratio of the quartet is fixed, and it depends only on the six `d_ij` and the bounds. The bounds depend on `d_ij` and `N`.
A second realization with another orientation array but the same `d_ij` has the same `(A,B,C)` and the same bounds.
This is the point that makes "prime orientations chosen separately" legitimate. The step is correct.

**Availability of certified quartets.** I checked each ingredient.
- `G_ij = <f_i,f_j> = 1 − 2 w(a_i−a_j)/W`, with total layer mass 1 when `W = log N`.
- Pair separation makes all `G_ij < 0`. `G` is PSD, so `λ_max(I−G) ≤ 1`, and `I−G` is entrywise nonnegative. Perron–Frobenius
  then gives `‖I−G‖ ≤ 1`, hence `λ_max(G) ≤ 2`, which is the frame bound. `frame_check.py` found a worst value of 1.874 over 317 random admissible Grams.
- The count of bad partners is correct. Astra's `2/τ^2` is loose: the self-term uses 1 of the budget 2, so fewer than `1/τ^2` partners are bad.
- The 4-mark indicator `(1/8)(1 + Σ_{pairs} + quartic)` was checked exhaustively. So each 2|2 matching has log weight
  at least `(1−7τ)W/8`. The surviving exponent `h_3 − h_2` is exactly the mass of the 2|2 layers of that matching, so this is `log B`.
- `|t_ij| ≤ sqrt(d_ij)/(4N^{1/4}) = e^{−W G_ij/4}/4` (identity checked), so `U = e^{τW/2}` bounds each product of two cofactors.
- At `τ = 1/32`, `(1−7τ)/8 − τ = 0.0664 > 0`.

**Gluing.** First take anchors `i, j, k`, pairwise good, placed at `∞, 0, 1`. Every extension `x` passing the four tests gets its chart value from the one
certified quartet `ijkx`. A missed `x` is placed through `y, z, w` taken from the recovered set, with `xyzw` certified.
The two realizations agree at `y, z, w`, so they agree at `x`. The argument is correct. Explicit constant: anchors exist once `m > 4097`.
The recovered set has at least `m − 8195` marks. The choice of `w` needs more than `8194` candidates. So **`m ≥ 16,390` suffices with Astra's constants**.
With sharp counts (pair failures `< 1/τ^2`, quartic failures `< 2/τ^2`) about `10^4` suffices.

**Deletions.** At each choice, also avoid the deleted partners of `x` (for `y`), of `x, y` (for `z`) and of `x, y, z` (for `w`). Then all six constraints of
`xyzw` are retained. The recovered set loses at most `3ηm`. The condition `6η < 1` is right, and at `η = 1/10` it needs `m ≳ 4.1·10^4`. For the
`o(m^2)` version, the principal subGram of the remaining rows still has negative off-diagonal and the frame bound, so it is fine.

**Pool strengthening.** The union bound runs over at most `4^r` allocations and `(2U+1)^3` coefficient triples, with each triple forcing
the middle product. By unique factorization, conditional on the exposed positions the middle `k ≥ r/16` positions are a uniformly
random `k`-subset of the remaining `K − r + k` primes. Since `binom(K−r+k, k) ≥ (m^η/2)^{r/16}`, the bound is
`exp(−Ω_η(m log m))`. Correct. It is peripheral, because it concerns the Paley profile only, and that profile is not forced on clusters.

**Numerical test on the record's explicit clusters.** Two remarks first.
- `growth/families.md` contains **no explicit 10-point cluster**. It states (§§0, 9) that none was found. In its place I used the
  shortest 10-point window on `R^2 = 5^2 13^2 17·29·37·41·53·61`, which has `C ≈ 25.7`.
- All the record's clusters have `C ≈ 5–6` in the `C = Φ N^{1/4}` normalization. They are therefore outside the theorem's normalization `Δ ≤ ½ N^{−1/4}`: the minimum `log d_ij / W` is 0.36–0.47, below 1/2. Only the local criterion and the algorithm can be tested on them, not the theorem itself.

Results (all quartets, exact arithmetic):

| cluster | m | quartets | Plücker & lcm formula | certified (a-priori bounds) | recovered = true | marks reconstructed (a-priori) | certified (oracle bounds) | marks reconstructed (oracle) |
|---|---:|---:|---|---:|---|---:|---:|---:|
| golden n=4 (`R^2=537606725`) | 8 | 70 | 70/70 | 0 | – | 3/8 (anchors only) | 10 | 8/8 |
| golden n=7,10,…,22 (six members) | 8 | 70 each | all | 4, 8, 8, 8, 8, 8 | all | 4/8 | 10–14 | 8/8 |
| record 8-point, C=5.673925 | 8 | 70 | 70/70 | 0 | – | 3/8 | 31 | 8/8 |
| record 7-point, C=5.269927 | 7 | 35 | 35/35 | 0 | – | 3/7 | 22 | 7/7 |
| 8-window, C=14.8 | 8 | 70 | 70/70 | 0 | – | 3/8 | 12 | 7/8 |
| 10-window, C=25.7 | 10 | 210 | 210/210 | 0 | – | 3/10 | 11 | 9/10 |

Every certified quartet gave exactly the true primitive ray, with no ambiguity and no mismatch. Every chart value
equalled the true cross ratio. So the recovery mechanism is sound. But **at m = 8–10 the reconstruction
from `d_ij` alone determines at most one mark beyond the three anchors**. The theorem's multiplicity
hypothesis (`m ≳ 10^4`) is essential. Even with oracle knowledge of the exact cofactor sizes, the
windows are not fully reconstructed.

**Novelty.** `growth/ptolemy.md` Theorem B identifies the content of the quartet system as the rational point of
`M_{0,m}` plus its boundary contacts. III.7 adds a converse-type rigidity: for `m ≳ 10^4` under pair separation, the boundary
contacts (`d_ij`) together with the cofactor height bound determine that point. I found nothing equivalent in the record by grep
("determin", "unique", "cross ratio", "PGL"). It is **new**, and it is distinct from Astra's Part I scalar reconstruction, which
uses the exact `t_ij`.

## III.8 — five-mark compatibility

**(III.6).** At a layer separating `{i,j}` from `{k,l,h}`, the uncrossed matching `ij|·` of each of the three quartets has
reduced valuation at least `v_p(Q)`. The two crossed matchings both have valuation `h_1 + h_2`, equal to that of the gcd, so their reduced quotients
`C_s, D_s` are units mod `Q`. The reduced congruences `C_s a_s ≡ D_s b_s (mod Q)` follow from `[ik][jl] − [jk][il] = [ij][kl]`
after division by the common gcd. The rational identity `(g_ik g_jl/(g_jk g_il))(g_il g_jh/(g_jl g_ih)) = g_ik g_jh/(g_jk g_ih)` gives
`C_1C_2D_3 = D_1D_2C_3` exactly. Multiplying the congruences and cancelling the unit gives `Q | a_1a_2b_3 − b_1b_2a_3`. This is bounded by
`2U_1U_2U_3 < Q`, so the identity holds. **Correct.** The random check found 0 failures of the telescoping, of coprimality and of
`Q | uncrossed` on 20,000 prime-power arrays.

**Five-mark consequence.** The indicator `(1/16)(1 + 10 pair + 5 quartic)` was checked exhaustively. The margin
`(1−15τ)/16 = 0.0479 > 3τ/2 = 0.0234` holds at `τ = 1/64`. Place marks 4 and 5 using `1234` and `1235`. Then apply (III.6) with anchor pairs `12`, `13`, `23` and the third
mark as auxiliary. This gives `1245`, `1345`, `2345`, so all five quartets are consistent. **Correct.**

**Gluing on regular quartets.** The finite-set step is stated loosely ("successive single-member replacements"). A clean
proof, which I checked, goes as follows. In a regular set `S` with chart triple `T = {1,2,3}`, place every `x` by `123x`.
- A quartet meeting `T` in two marks, say `12xy`, lies in the regular five-set `123xy`, so it is consistent.
- A quartet `1xyz` lies in the five-set `12xyz`. Its frame `1, 2, x` places `y` and `z` by `12xy` and `12xz`, which are already consistent. So `1xyz` is consistent.
- A quartet `wxyz` lies in the five-set `1wxyz`, which is handled the same way.

The infinite gluing (transitions through auxiliary triples `J`, cofinite agreement, uniqueness of the projectivity from
three distinct values) is correct. Explicit size: an extension of a `d`-set excludes `2(d + C(d,3))/τ^2 = 8192(d + C(d,3))` marks. For the
7-sets used this is about `3.4·10^5`, so `m ≳ 5·10^5`.

**Significance.** This is a **negative result about a method**. It says that cross-ratio compatibility across quartets carries no information beyond the single
quartet tests. Large multiplicity cannot be exploited through consistency of independently solved quartets. Astra
draws this conclusion correctly.

## III.10 — what remains

Everything in III.10 is correct, elementary linear algebra.
- **Parametrization.** `V_i = λ_i M u_i` gives `t_ij = det(M) λ_iλ_j h_ij`.
- **Realization.** The saturation `Π ∩ Z^m` has a basis with primitive Plücker vector, and every integral realization is `M B` with
  `det M = d`.
- **Primitivity.** Necessity: mod `p | d`, a rank-1 `M` has a kernel line that must miss every column. Sufficiency: use a CRT-lifted primitive row, completed to a unimodular matrix, with the second row
  multiplied by `d`. Checked against brute force on 47 random cases.
- **Quadratic form.** `S = M^tM` is determined by three directions, because a binary quadratic form has at most two projective zeros.

Combined with III.6, where exact norms force isotropy and hence genuineness, III.10 restates the residual existence problem for a
given labelled profile as follows: find rational scales `c, s_i` making all `c s_i s_j h_ij` nonzero integers of the
prescribed small size, then satisfy the primitivity condition, and then fit one positive form. No theorem, bound or obstruction is
proved about this problem. It is a precise problem statement.

## Is "sufficiently large" explicit? Bound or rigidity?

- **Explicit:** yes. It is absolute and independent of `N`, `W` and the radius, so no template-dependent threshold arises. III.7
  needs `m ≥ 16,390` (Astra's constants), or `≈ 10^4` with sharp counting, and `≈ 4.1·10^4` with `η = 1/10` deletions. III.8's gluing needs
  `m ≈ 5·10^5`. Astra never states these numbers; they should be stated.
- **No bound:** none of III.7, III.8 or III.10 bounds `m`. III.7 is a uniqueness assertion about hypothetical
  clusters of size about `10^4`. It does not count clusters, and it does not exclude an exceptional parent or assignment. III.8 shows that
  the quartet tests glue, which removes a possible source of contradictions. III.10 reformulates.
- **Standing cautions:** none is violated.
  - No fixed-`S` Subspace input is used.
  - Gaussian gcds are carried in `g_ij`, and units are absorbed explicitly.
  - No drifting phase is set to zero; only rational directions are used.
  - The threshold in `m` is uniform.

## Significance

Correct and new, but **peripheral to the goal**. III.7 is a clean rigidity theorem: in a hypothetical
cluster of size about `10^4` or more, the labelled separation divisors determine the rational point of `M_{0,m}`. III.8 is a no-go for
exploiting multiplicity through cross-ratio consistency. III.10 is a problem restatement. Nothing here moves the
growth rate `log R / log log R`, the constant 2/3, or the uniform bound `M(C)`. Nothing changes the status of Question A (`L_8`; `m = 8` is far
below the threshold, and my computations show that reconstruction from `d_ij` fails at `m = 8`), of Question B, or of route 5.
