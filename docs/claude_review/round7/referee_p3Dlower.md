# Referee report on `round7/p3Dlower.md` (P3-all\*: can designed prime residues be forced to cost `>> M log M`?)

Referee: adversarial, independent. Scripts are in `round7/referee_p3Dlower_checks/`; none imports the
author's code (the author's scripts were only re-run, separately, to confirm their stated outputs).
Each script runs from that directory with `python3 <script> [args]`; outputs are stored next to it as
`*_out.txt`. The repository and the round-5/6/8 folders were not modified.

Labels: **proved** = re-derived here with all quantifiers; **exact** = exact finite computation;
**evidence** = a computation no proof uses. `§n` is a section of `p3Dlower.md`, `R§n` of this report.

---------------------------------------------------------------------------------------------

## 0. Verdict

**The numbered mathematics survives. Several headline sentences do not, and the evidence table is wrong.**

* *Correct (re-derived, and exact-checked where finite).* Theorem 1.1 (both items, as formally stated),
  Theorem 2.1 (all four items), Proposition 3.1, Theorem 3.2 (under GRH), Definition 4.1,
  Lemma 4.2 (now checked at **every** level `a <= A_q`, not only level 1), Theorem 4.3, Corollary 4.4
  and Lemma 5.1. No mathematical error was found in any of their proofs. The author's scripts reproduce
  (`check_theorems.py`: `ALL OK`; `designs.py` reproduces `designs_out.txt`).
* *Not correct as stated.*
  1. **Evidence table (§6) is computed with a unit-bookkeeping bug** (major). `designs.py` classifies a
     collision by `<v,ell>/(m_q/4)` without normalising by `iota = log i`, which is `m_q/4` at some moduli
     and `3m_q/4` at others (21 vs 13 of the odd `q < 150` with the author's generators). So its
     "per-unit" demand mixes the units `i` and `-i` across moduli. Proof that the numbers are wrong: the
     author's `GEN` entry at `k = 1` is `2.51 log M` at `M = 50`, but genuine labels obey
     `m_(e_j,s) <= |D_s(pi_j)| <= sqrt(2 p_j)`, i.e. at most `1.68 log M` there (likewise at all five `M`:
     author `2.51-3.03`, hard cap `1.60-1.68`). With correct bookkeeping
     (R§3) genuine labels sit **at their own cap, `~1.6 k log M`, far below the random ones**, coherent
     labels at `~0.8 k log M` (not `1.5 k`), and twisted labels at 2-3.4 times the coherent value at
     `k = 1` (above genuine). All three "readings" of §6 are artefacts.
  2. **Theorem 1.1(2) heading "only as a veto", and the summary "at `q = +-1 mod 8` it can only veto a
     collision", are false** (major as a claim, harmless for growth). At every `q = +-1 mod 8` with
     `(p_j/q) = 1` the parity *forces* the monomial `(m_q/8) e_j` to collide, for every
     parity-respecting label, while an unconstrained label (`ell_j = 1`) avoids it. At `q = 7` this is the
     norm-1 monomial `e_j`: all 553 primes `p = 1 mod 4` below `2*10^4` with `(p/7) = 1` have
     `7 | D_s(pi)` for some `s` (exact). The probability statement of 1.1(2) and Lemma 5.1 are correct
     because they assume `gcd(v, m_q/2) = 1`, which the summary and the §0 table drop.
  3. **"(M2) fails, proved"** is overstated (major). Nothing proved here, nor in the cited parts of
     `p3lower.md`, shows that "pigeonhole/entropy cannot force a height violation". Lemma 5.1 covers
     only sets of fewer than `m_q/8 <= X/8` monomials at **one** modulus, while the relevant sets have
     `(2r)^k` elements. §5.1 also relies on two statements that the round-7 referee of `p3lower.md`
     rejected: "parity gains at most a factor 2 per modulus" (false as phrased) and "counting over all
     moduli at once gives only B.5" (B.5 is single-modulus and two-term).
  4. **"(M1) is false on twin profiles"** is overstated (major, headline). Corollary 2.2 shows that pairs
     do not constrain the `r - 2M` non-twin labels and the `M` sums `U_x`. They do constrain all `2M`
     twin-column prime residues, through `Lambda = kappa - F`. That is a localisation of (M1), not a
     refutation. The author's own heuristic (last row of the §5.2 table) says the twin channel alone
     produces the random rate.
  5. **§0 item 2, "Certified designs cannot reach `O(M log M)`" and "any fake at `W = O(M log M)` would
     have to use uncertified labels"** (major). Under the §3 definition of "certified", genuine labels
     (actual clusters) are certified, with zero slack. Read literally, the sentence would therefore prove
     growth for actual circles, which the author rightly does not claim. What is proved is narrower: the
     crude certificates `m <= |D|` cannot *prove* membership of these design families at `O(M log M)`.
  6. **Class definition** (major, validity): "P3-all\* ... with ... the completion requirement for `v`
     outside `L` (`sharp.md` §6.1). Actual clusters are members." In `sharp.md` §6.1 the completion
     requirement is the union-bound sum `sum_v min(1, 4 m_v e^(-w/2))` being small. That is a
     *sufficient* condition for constructing angles, not a necessary one. Actual clusters violate it
     (`p3G.md`: the genuine sums are 4052 and 15190 on verified sets). No proved result here uses it as a
     necessary condition, so nothing breaks. But the class must be stated as `F^all_X` (slab conditions
     (1.2) for all `v`), with the union bound only as a construction tool.
* *Further minor issues* are listed in the table below: novelty of Theorem 2.1, the vacuous
  `q = +-1 mod 8` part of Proposition 3.1's consequence, the levels and non-degeneracy in
  Corollary 3.3, the random-fibre constant, the `p3G` citation, the roadmap (I2), and wording.

**No growth claim is made, and none is refuted.** The bracket `3 M log M - O(M) <= W_min^all <=
(66 A + o(1)) M log^2 M/loglog M` is unchanged. The draft is a correct but modest partial barrier. Its
proved content is: two exact identities (Theorems 1.1, 2.1); one GRH lemma on coherent labels with the
forced parity (Theorem 3.2); and explicit caps for one parity repair (Theorem 4.3). Each of these limits
one *design* or *certification* route. None of them bounds the class from below.

### Issues

| # | where | severity | issue | correction |
|---|---|---|---|---|
| 1 | §6 evidence table and its "Reading"; `designs.py` | **major** (evidence) | unit of a collision read as `u = <v,ell>/(m/4)` instead of `s` with `<v,ell> = s iota`; `iota = 3m/4` at about 40% of moduli, so units `i`, `-i` are mixed across `q`; the `GEN` column exceeds the genuine cap (`2.51` vs `<= 1.68` at `M = 50`), which is impossible | replace the table by R§3; genuine `~1.6 k log M` (at its cap), coherent `~0.8 k`, twisted 2-3.4 x coherent at `k = 1` and near random at `k = 2, 3`. "Genuine sit with random" and "twisted match coherent at `k = 1`" are withdrawn |
| 2 | Thm 1.1(2) heading; §0 table row 1.1; Author's summary; §5.1 "the parity forces *nothing*" | **major** (false as phrased) | the parity forces collisions: `v = (m_q/8) e_j` collides at every `q = +-1 mod 8` with `nu_j(q) = 0`, for every parity-respecting label (R§2.1); at `q = 7` this is `v = e_j` | "at `q = +-1 mod 8` the parity vetoes collisions of `v` with `<v,nu>` odd and doubles the collision probability of the others; for `v` with `gcd(v, m_q/2) > 1` it can force collisions (e.g. `(m_q/8) e_j`); these forced collisions are torsion-type and harmless for the height condition" |
| 3 | §0 (M2), §5.1, Author's summary "both failures are proved" | **major** (overstated) | "(M2) cannot work" is a meta-claim, not a theorem; Lemma 5.1 is single-modulus, `|V| < m_q/8`; §5.1 cites `p3lower.md` B.3 ("factor 2") and B.5 ("over all moduli at once"), both rejected as phrased by `referee_p3lower.md` (issues 2 and 8) | label (M2)'s failure **heuristic**: "the per-modulus and the two-term single-modulus counts give only average demands (Lemma 5.1, `p3lower` B.3(3)-(5), B.5); no counting argument forcing concentration is known" (this is `sharp.md` §6.7's sentence) |
| 4 | §0 (M1), §5.2 heading, Author's summary | **major** (headline) | Cor. 2.2 decouples `(ell_J, U)`, but the `2M` twin-column residues `ell_(f(x)) = eps_x Lambda_x`, `ell_(s(x)) = U_x - eps_x Lambda_x` are constrained through `Lambda = kappa - F`; (M1) is localised, not refuted, and the author's own heuristic (§5.2 last row) says the twin channel yields the random rate | "(M1) cannot constrain the `r - 2M` non-twin residues or the twin sums; it acts only through the `M` twin labels `Lambda` (M1')" |
| 5 | §0 item 2 | **major** (overstated) | genuine labels are "certified" in the §3 sense (`m_(v,s) | D_s(A_v)`), so the literal claim would imply growth for actual circles; only the crude certificate's *inability to prove membership* is shown; for twisted labels only the cap is shown | "the certificates available (additive, exact coherent, twisted) cannot *prove* membership at `O(M log M)`; whether coherent or genuine-size labels *are* members at `O(M log M)` is the open (I3)/`H` question" |
| 6 | Conventions (P3-all\* definition), §0 "Validity for actual circles" | **major** (class validity) | the union-bound "completion requirement" of `sharp.md` §6.1 is sufficient, not necessary; actual clusters violate it (`p3G.md`: sums 4052, 15190), so "actual clusters are members" is false for the class as written | define P3-all\* `= F^all_X` (all slab conditions (1.2)); use the union bound only for constructions. No proved result changes (all lower-bound-type uses are of the per-unit height condition, which is necessary) |
| 7 | Thm 2.1, Cor. 2.2, Rem. 2.3 | minor (novelty) | the decoupling is `sharp.md` Lemma 6.3(2) (fibre parametrised by non-twin labels, `u`, `t`; twin labels solved), with `c'` from Lemma 6.2, and `sharp.md` §6.6 ("with twins, `M` labels are always solved from the point design") | cite; the new content is the explicit basis `(ell_J, U, Lambda)` and identity (2.2) |
| 8 | Prop. 3.1 "Consequence"; §5.2 table row "additive ... yes if parities fit" | minor (hidden assumption) | at `q = +-1 mod 8`, labels (3.1) are parity-respecting only if `nu(q) in {0, lambda mod 2}`; for random primes this holds at **at most one** such `q` (exact, R§2.4), and under GRH it fails at all `q = +-1 mod 8 <= X` once `X >> log^2 p_max` (Thm 3.2's argument with `n = p_j p_k`); so the consequence with "all odd `q <= X`" is vacuous | take `Q = {q = +-3 mod 8}`: `W >= (1/8 - o(1)) M X/(h+2)`; growth conclusion unchanged |
| 9 | Remark 1.2 "a parity-based lower bound can only use `q = +-1 mod 8`" | minor (overstated) | at `q = +-3 mod 8` the parity fixes the **parity of the collision unit**, `s = <v, nu(q)> mod 2` (R§2.1, exact); the per-unit slab conditions (1.1), (1.2) and the per-unit height condition see this | restrict to all-unit quantities (`m_v`, union-bound sums); per-unit effects are at most a constant factor in the height condition |
| 10 | Cor. 3.3 | minor | (3.3) is a level-1 hypothesis, but the certificate bounds `m_(v,s)` at all levels; and for general Gaussian `zeta_j` `D_s(Z_v)` can vanish (e.g. `zeta_1 = zeta_2`, `v = e_1 - e_2`) | assume (3.3) at the top level `q^(A_q)` and `zeta_j/conj zeta_j` multiplicatively independent modulo units (as in Lemma 4.2) |
| 11 | §5.2 table, random fibre row | minor | `(2 + o(1)) log^2 M/loglog M` is proved only at `A = 1`: Prop. 6.5 gives `>= (2 - eps)` for every `X >= M`, Thm 6.4 gives `<= (2A + o(1))` | "between `(2 - o(1))` and `(2A + o(1))`" |
| 12 | §5.2 table, genuine row; §0 item 5.2 | minor (citation) | `p3G.md`'s "completion fails at poly sizes / `tau ~ log^2 M/loglog M`" is not accepted by `referee_p3G.md`: it holds for random prime sets and in a non-adaptive model, and selected primes defeat the computed closures | mark as "random prime sets / model; open for selected primes" |
| 13 | §5.3 (I2), "the only step proved here" | minor (roadmap) | twisted coherence is certified (global object, explicit caps) with `Norm zeta_j ~ X << p_j`, so "(I2) certified implies genuine size" is false for it; the chain (I1)-(I3) needs an extra branch "approximately certified implies cap growth plus (I3)" | state (I2) for exact coherence only; add the twisted branch as open |
| 14 | Thm 3.2, Rem. 3.4 | minor (attribution) | the argument is the standard GRH/Chebotarev bound behind lower bounds for pseudosquares (Ankeny, Lagarias-Montgomery-Odlyzko, Bach, Schinzel) | cite as a variant of the classical GRH pseudosquare bound |
| 15 | Def. 4.1 title "the best certified repair"; "admissible in P3-all\*" | minor (wording) | "best" is not proved; "admissible" means only that the residue data have genuine norms, not membership | "a natural repair"; "admissible residue data" |
| 16 | Cor. 4.4 "Nothing is certified for `k > k_1`" | minor | (4.3) is a sufficient condition obtained by lossy steps (`0.7/k <= 0.7`, `w(v) >= k tau`) | "(4.2) does not certify `v` of norm `k > c k_1` supported on columns of weight `~tau`" |
| 17 | Cor. 2.2 / Author's summary "every monomial is a combination of the same sparsity" | minor | `v = e_(s(x))` is 1-sparse, `vt` is 2-sparse; only `||vt||_1 <= 2||v||_1` holds | "of comparable norm (factor 2)" |

---------------------------------------------------------------------------------------------

## 1. Re-derivation, result by result

### 1.1 Theorem 1.1 (parity localisation)

*Item 1 (`q = +-3 mod 8`).* **Correct.** `m_Q = q^(A-1)(q - chi(q))` and `q - chi(q) = 4 mod 8`
exactly when `q = +-3 mod 8`. So `iota = log i` is odd there, and adding `iota u_j` flips the parity
of `ell_j`. Then `<v, ell> = <v, ell'> + iota <v,u>` at every level. The reduction of `iota` mod
`m_(q^a)` is `log i` at level `a`, because a generator reduces to a generator. The four events
`s in Z/4` are disjoint because `4 | m_(q^a)`. So `max_s a_(q,s)` and `m_v` are unchanged.

*What the statement does not say.* The units are permuted by a shift `<v, u(q)>` that depends on
`q`. Within the parity class, only shifts by `2 iota` (sign changes) remain free. Hence for a
parity-respecting `ell` at `q = +-3 mod 8`, a collision of `v` with `i^s` forces

```text
s  =  <v, nu(q)>  mod 2                                                    (R1)
```

because `<v, ell> = <v, nu> mod 2` and `iota` is odd. (Exact check in R§2.1.) For actual circles (R1)
says: `q | Im` or `q | Re` of `A_v` if `prod p_j^(v_j)` is a QR mod `q`, and `q | Re +- Im` otherwise.
Membership is defined by *per-unit* conditions (`sharp.md` (1.1), (1.2), and the per-unit height
condition), which see (R1). So Remark 1.2's "a parity-based lower bound can therefore only use moduli
`q = +-1 mod 8`" is true only for all-unit quantities. The per-unit effect on the height condition is
at most a constant factor (two units instead of four), so this does not matter for growth (issue 9).

*Item 2 (`q = +-1 mod 8`).* **The formal statement is correct.** Write `ell = nu + 2 ell''`. If
`gcd(v, m/2) = 1`, then `<v, ell''>` is uniform on `Z/(m/2)`. Since `8 | m`, all four elements of
`(m/4)Z/m` are even. This gives `8/m` or `0`, against `4/m` without the parity; exact enumeration for
`r = 2` and `q = 7, 17, 23, 31, 41` agrees (R§2.2).

**The heading "one bit, and only as a veto" is false** (issue 2). Let `nu_j(q) = 0`, so `ell_j` is even.
Then `(m/8) ell_j` lies in `(m/4) Z/m` for every even `ell_j`. So `v = (m_q/8) e_j` collides for every
parity-respecting label, while the unconstrained label `ell_j = 1` gives `(m/8) * 1`, an odd multiple
of `m/8`, which does not collide. Examples:

| `q` | `m_q` | forced monomial |
|---|---|---|
| 7 | 8 | `e_j` |
| 17 | 16 | `2 e_j` |
| 23 | 24 | `3 e_j` |
| 31 | 32 | `4 e_j` |

This is consistent with 1.1(2), whose probability statement assumes `gcd(v, m/2) = 1`. It is
inconsistent with the heading, with §5.1 ("at each modulus the parity forces *nothing*"), with the
§0 table entry for 5.1 (which drops the gcd hypothesis), and with the Author's summary ("at
`q = +-1 mod 8` it can only veto a collision").

For actual circles this is the familiar fact that a QR prime has `pi/conj pi` a square in `T_q`. Its
`(m/8)`-th power is then a fourth root of unity. Example: at `q = 7`, `29 = 2^2 + 5^2` gives
`D_3 = 2 + 5 = 7`.

These forced collisions are torsion-type: the forced moduli of `c e_j` divide `prod_(m_q | 8c) q`. That
product is `<= exp(d(8c) log(9c))`, which is negligible against `c w_j/2`. So the error does not affect
any growth statement, but the sentence must be corrected.

### 1.2 Theorem 2.1 and Corollary 2.2 (twin decomposition)

**Correct.**
* *Item 2.* `a_(y,f(x)) = sigma_(yx) + eps_x delta_(xy)` gives `<a_y, ell> = F_y + Lambda_y`.
* *Identity (2.2).* It follows from `v_f ell_f + v_s ell_s = v_s U + (v_f - v_s) eps Lambda`.
* *The map `v -> vt`.* It is unimodular (lower-triangular blocks with `+-1` on the diagonal), and
  `||vt||_1 <= 2||v||_1 <= 4||vt||_1`.
* *Item 4.* `sum_j a_yj nu_j = F_y + nu_(f(y)) mod 2`. Every kappa with this colouring is realised by
  every parity-respecting `(ell_J, U)`.

All of this was checked exactly on 40 random twin profiles, with random moduli and labels, and
exhaustively over small `v` (R§2.3).

Three remarks.
1. *Scope of "free".* The free coordinates are `ell_J` (`r - 2M` actual prime residues) and the
   *sums* `U_x = ell_(f(x)) + ell_(s(x))`. Both residues of every twin pair depend on `Lambda`. In terms
   of prime residues, the pair data therefore constrain `2M` of the `r` columns, not "the twin labels"
   as an abstract `M`-dimensional block. (M1) concerns prime residues. Corollary 2.2 refutes it for
   the non-twin columns only (issue 4).
2. *Sparsity.* `v = e_(s(x))` is 1-sparse, but `vt = (U_x, -eps_x Lambda_x)` is 2-sparse. So the
   summary's "same sparsity" should read "comparable norm (factor 2)" (issue 17).
3. *Novelty.* `sharp.md` Lemma 6.3(2) already parametrises the random fibre by the non-twin labels,
   the vector `u` (the `U`-coordinates) and the common shift `t`, with the twin labels solved. Lemma 6.2
   introduces `c'(v) = eps (v_f - v_s)`. The last paragraph of `sharp.md` §6.6 states the consequence
   in words. Theorem 2.1 is a clean basis-level restatement of these (issue 7).

### 1.3 Proposition 3.1 (additive labels)

**Correct as a statement.** The steps:
* *Pigeonhole.* `binom(r',k) > 2kH + 1` sums lie in `[-kH, kH]`, so two distinct `k`-subsets have
  the same sum. This gives `v = 1_S - 1_(S')` with `||v||_1 <= 2k` and `<v, lambda> = 0`.
* *Collisions.* `v` collides at `q^(A_q)` with the unit `i^(<v,u(q)>)` for every `q` in `𝒬`, so
  `m_v >= Pi_𝒬`.
* *Height.* `prod_s nu_s^(-1) = 2`, with `nu_s^(-1) = 1` for even `s` and `sqrt 2` for odd `s` (this
  is `sharp.md` (1.2); confirmed by `W - i^s conj W = 2iY, (1-i)(X-Y), 2X, (1+i)(X+Y)`). The per-unit
  height condition therefore gives `Pi_𝒬 <= 2 e^(2 w(v)) <= 2 e^(4k w_*)`.
* *Weights.* `W >= (floor(r/2) + 1) w_* >= (r/2) w_*`.

**The "Consequence" is vacuous as stated** (issue 8). For a member of P3-all\*, (3.1) at a modulus
`q = +-1 mod 8` forces `g_q lambda_j = nu_j(q) mod 2` for all `j`, since `iota` is even there. That
means `nu(q) = 0` if `g_q` is even and `nu(q) = lambda mod 2` if `g_q` is odd. If `lambda_j = lambda_k`
mod 2, then `(p_j p_k/q) = 1` at every such `q` in `𝒬`.
* With `r >= 3` such a pair exists. Under GRH, the argument of Theorem 3.2 (with `n = p_j p_k`, never
  a square or twice a square) excludes this as soon as `X >= c_1 (log 2n)^2 (loglog 3n)^2`.
* Exactly: for random primes (`r = 20, 60`; `p ~ 10^6, 10^7`; `X = 300, 1000`) all `nu(q)` at
  `q = +-1 mod 8` are distinct. **At most one** such modulus is compatible with any single `lambda`, and
  no pair `j < k` has equal symbols at all of them (R§2.4).

So the only admissible additive designs use (essentially) `𝒬 subset {q = +-3 mod 8}`. Then
`log Pi_𝒬 = (1/2 - o(1)) X` and `W >= (1/8 - o(1)) M X/(h+2)`. The growth conclusion (order `M^2/h`)
is unchanged.

### 1.4 Theorem 3.2 and Corollary 3.3 (exact coherence under GRH)

**Theorem 3.2 is correct under GRH.**

*Parity step.* At `q = +-1 mod 8`, `m_q/4` is even. Also, `log(z/conj z)` is even iff
`(Norm z/q) = 1`; this was re-derived for split and inert `q` and checked exactly by the referee's own
discrete logs. Together these give `(p Norm zeta/q) = 1` at every `q = +-1 mod 8`, `q <= X`,
`q` not dividing `n`.

*Chebotarev step.* Two points need adjusting.
* The equivalence "Frobenius `sigma` iff `q = +-1 mod 8` and `(n/q) = -1`" holds only for `q` not
  dividing `n`. When `q^2 | n`, the prime `q` may be unramified with Frobenius `sigma` while
  `(n/q) = 0`. The proof handles this correctly by counting more than `omega(n)` primes.
* `d_L <= 8 * 4n * 8n = 256 n^2`, so the stated `log d_L <= 3 log(32 n)` is valid.

*Error term.* With `x = c_1 (log 2n)^2 (loglog 3n)^2`, the main term `x/(4 log x)` is of order
`(c_1/8)(log n)^2 loglog n`. The error `O(x^(1/2)(log d_L + log x))` is of order
`sqrt(c_1)(log n)^2 loglog n`. So a large `c_1` works.

*Size step.* `p | t` gives `p | Norm zeta`.

Evidence (no proof uses it): for all 231,782 pairs `(p, zeta)` with `p < 3000` and `Norm zeta < p`,
the least admissible `q` is at most `2.7 (log n)^2` (R§2.4).

*Attribution* (issue 14). This is the standard GRH argument behind lower bounds for pseudosquares and
least non-residues (Ankeny; Lagarias-Montgomery-Odlyzko; Bach), here restricted to the class
`q = +-1 mod 8`. Remark 3.4 correctly calls it the rigorous form of the `sharp.md` §6.6 heuristic; it
should also cite it as classical in method.

*Remark 3.4, unconditional case.* The remark says the unconditional version "says nothing useful".
That holds for a threshold `X >= n^(c')` with `c' >= 1`, which the known Linnik/LMO constants give.
The referee did not examine whether two-character Burgess-type bounds give `c' < 1`. That would matter
when `p_j Norm zeta_j <= X^(1/c')`.

**Corollary 3.3: two minor gaps** (issue 10).
1. *Levels.* (3.3) is a level-1 hypothesis, but the certificate bounds `m_(v,s)`, which includes the
   levels `a >= 2` (`q^2 <= X`). The certificate needs coherence at the top level `q^(A_q)`.
2. *Non-degeneracy.* For general Gaussian integers `zeta_j`, `D_s(Z_v)` may vanish. For example, if
   `zeta_1 = zeta_2` then `Z_(e_1 - e_2) = |zeta_1|^2` and `D_0 = 0`. Then `v` collides at every level
   and modulus, and no certificate exists. Lemma 4.2 avoids this by taking the `zeta_j` to be Gaussian
   primes over distinct split primes. Corollary 3.3 needs the same, or multiplicative independence of
   the `zeta_j/conj zeta_j` modulo units.

Neither gap changes the conclusion "zero slack": a weaker hypothesis can only make certification worse.

*Scope of "never certify completion".* This concerns the union-bound sum, which is sufficient for
completion, not necessary (issue 6). It shows only that coherence plus the crude bound `m <= |D|` cannot
prove completion. It says nothing about whether completable coherent designs exist. The author says
this in the last sentence of Corollary 3.3, but §0 item 2 drops it (issue 5).

### 1.5 Definition 4.1 - Remark 4.5 (twisted coherence)

**Definition 4.1** is consistent.
* *`t(q) = 0` iff `q = +-3 mod 8`.* This follows from `(2/q)`, and `rho_0/conj rho_0 = i`.
* *Parity.* `log(rho_t/conj rho_t)` is odd iff `(ell_t/q) = -1`, so the twisted labels have parity
  `nu_j(q)`. They are therefore residue data with genuine norms. "Admissible in P3-all\*" should read
  "admissible residue data": membership also needs the slab conditions (issue 15).
* *"The best certified repair".* Nothing proves that this repair is best (issue 15).

**Lemma 4.2** is correct.
* *The identity.* `W - i^s conj W` equals `D_s(W)` times `2i`, `1 - i`, `2` or `1 + i`.
* *`D_s != 0`.* The valuation argument works because the `p'_j` are distinct split primes, different
  from the `ell_t`.
* *The check, extended.* The author checked the lemma at level 1 only. The referee checked it at
  **every** level `a <= A_q` with `X = 200`. This covers `3^4`, `5^3`, `7^2`, `11^2` and `13^2`, uses
  exact discrete logs in `T_(q^A)` built from scratch, takes all `v in [-2,2]^5 \ {0}`, and gives 20,398
  level-collisions. Every one satisfies `collision <=> q^a | D_s(W)`, and `D_s != 0` throughout (R§2.3).

**Theorem 4.3** is correct.
* *Groups `t >= 1`.* The prime powers in a group are pairwise coprime divisors of one non-zero
  integer, of size at most `sqrt2 |Z_v| ell_t^(|d|/2)`.
* *Group `t = 0`.* `q^a | D_(s-d)(Z_v)`, so the product divides `prod_(s') D_(s')(Z_v)`. In fact
  `|XY(X^2 - Y^2)| <= |Z|^4/4`, so the stated `2|Z|^4` is loose.
* *Sum over `d`.* `sum_(|d| <= k) |d|/2 = k(k+1)/2`.
* *Numerics.* On the same data the cap is never violated. The largest ratio `log m_(v,s)/cap` is
  `0.14`, so the cap is far from tight at these sizes.

**Corollary 4.4** is correct. The algebra was re-derived symbolically: solving the linear inequality
for `k` gives exactly (4.3) (R§2.4). The completion variant, `tau -> tau/4 - (log(4er)+1)/2`, was also
re-derived.

"Nothing is certified for `k > k_1`" overstates what the lossy reduction shows (issue 16). What is true
is this. Since `log|Z_v| >= (k/2) log X` (all `p'_j > X`), the cap exceeds the budget
`w(v)/2 ~ k tau/2` once `k >~ tau/(T log X)`, for `v` supported on columns of weight `~tau`. So the
certified range is `Theta(tau/(T log M))` for such `v`, as claimed in order of magnitude.

*Size of `T(X)`.* The GRH bound `ell_T <= c(log X)^2` and the resulting
`T = O(log^2 X/loglog X)` are right. Also `T(X) -> infinity` unconditionally: for any finite list of
primes there are `q` (CRT plus Dirichlet) at which all of them are residues.

**Remark 4.5** is correctly labelled heuristic, and its arithmetic is right. With `k_*` defined by
`T h k^2 = k log M/loglog M`, the demand per unit of norm at `k_*` is `log^2 M/loglog M`. Two further
remarks.
* If `T >= log M/(h loglog M)`, which is allowed under GRH, then `k_* < 1`. The cap never binds, and
  the heuristic gives the random rate from `k = 1` on.
* The supporting evidence cited ("twisted behave like coherent at `k = 1`") is an artefact of the
  bookkeeping bug. Corrected, twisted labels have 2-3.4 times the coherent demand at `k = 1` (R§3).

### 1.6 Lemma 5.1 and §5 (the mechanism)

**Lemma 5.1** is correct. Each `v` with `gcd(v, m_q/2) = 1` collides with probability at most `8/m_q`
(`4/m_q` at `q = +-3 mod 8`), and the union bound finishes. Its hypotheses make it nearly vacuous for
the question:
* `|V| < m_q/8 <= (X+1)/8`, whereas the monomials of norm `<= k` number about `(2r)^k`. With `X = 3M`
  and `r >= 2M`, this is not even all of norm 1.
* It is a single-modulus statement.
* The gcd hypothesis excludes the forced torsion collisions of R§1.1.

**§5.1, "(M2) as a counting argument cannot work"** (issue 3). This is not a theorem and should not
be in the "proved" column of §0 or in "both failures are proved".
* *"Parity changes their count by at most a factor 2 (`p3lower.md` Prop. B.3)."* `referee_p3lower.md`
  issue 2 shows this is false per modulus: the ratio is unbounded for `M <= m_Q < 2M`, and equals 4 and
  6.5 in finite cases. Only the additive bound B.3(3)-(5) holds.
* *"Counting over all moduli at once, the forced collisions are the ones of `p3lower.md` Prop. B.5."*
  `referee_p3lower.md` issue 8 shows that B.5 covers only two-term monomials at single moduli.
  Simultaneous (CRT) pigeonhole and longer monomials are "unproved generalisations".
* *Example.* Proposition 3.1 is itself a simultaneous-CRT pigeonhole, and it forces strong
  concentration, on additive designs. Whether counting forces concentration is therefore
  design-dependent. For unstructured designs, the full-CRT pigeonhole (`binom(r',k) > prod_q m_q/4`)
  gives only `W >~ r log(r/k)`, which is again the trivial count. This easy fact supports the author's
  intuition but does not prove the meta-claim.

The correct status is: **heuristic, and already stated in `sharp.md` §6.7** ("Pigeonhole forces only
average demands. The height condition is satisfiable for each monomial separately.").

**§5.2 (M1').** (M1') is a correct reformulation via (2.3), with one caveat. "At all but
`O(||w||_1 + 1)` moduli" uses `tau = O(log M)` and moduli of size `M^(Theta(1))`; this assumption
should be stated.

**§5.2 table.**
* *Additive row.* "Yes if parities fit" hides that they essentially never fit at `q = +-1 mod 8`
  (issue 8).
* *Random fibre row.* `(2 + o(1))` holds at `A = 1` only (issue 11).
* *Genuine row.* It cites `p3G.md` claims that `referee_p3G.md` did not accept: the claims hold for
  random prime sets and in a non-adaptive model, while selected primes defeat the computed closures
  (issue 12).
* *Last row (`D_A` + coherent free labels).* This is a heuristic *for* the mechanism (M1')+(M2). It
  contradicts the headline tone of "(M1) is false" (issue 4).

**§5.3 roadmap** (issue 13). Twisted coherence is a certified object (explicit global
`zeta_j, rho_t`, explicit caps) with `Norm zeta_j ~ X`, far below `p_j`. So "(I2) certified implies
genuine size" is false for it. The chain needs a third branch: approximately certified designs, with a
growing cap, plus an (I3)-type tail statement for them. The sentence "(I2) is the only step proved here"
should read "(I2) is proved for exact coherence under GRH".

---------------------------------------------------------------------------------------------

## 2. Exact checks (referee's own code; `referee_p3Dlower_checks/`)

All exact checks below print `ALL OK`. The arithmetic in `T_q` and `T_(q^A)` (`tq.py`, and
`check_twist.py`) was written from scratch. Collisions are cross-checked against direct divisibility
`q^a | D_s`.

### 2.1 `check_parity.py` (Theorem 1.1, and the counterexample to "only as a veto")

* **A.** Theorem 1.1(1): 5,277 random cases (`q = +-3 mod 8`, `q < 140`, `r <= 6`). Adding `iota u`
  shifts the collision unit by exactly `<v,u>`. It never creates or destroys a collision.
* **A'.** In every collision under a parity-respecting label at `q = +-3 mod 8`, the unit parity is
  `s = <v, nu> mod 2`. This is (R1) of R§1.1.
* **B.** Theorem 1.1(2), by exhaustive enumeration over all parity-respecting `ell in (Z/m)^2` for
  `q = 7, 17, 23, 31, 41` and all `v in [-5,5]^2` with `gcd(v, m/2) = 1`. The number of hits is
  exactly `8|P|/m` when `<v,nu>` is even, and `0` when it is odd.
* **C.** The parity forces collisions. For `q = 7, 17, 23, 31, 41, 47, 71, 73, 79`, the monomial
  `(m_q/8) e_j` collides for **every** even label, while the label `1` avoids it. Genuine
  confirmation: all 553 primes `p = 1 mod 4` with `p < 2*10^4` and `(p/7) = 1` have `7 | D_s(pi)` for
  some `s`. For example:
  * `29`: `(2,5) -> D = (5, -3, 2, 7)`;
  * `37`: `(1,6) -> D = (6, -5, 1, 7)`;
  * `53`: `(2,7) -> D = (7, -5, 2, 9)`.

### 2.2 `check_legendre.py`

For all odd primes `q < 200` and all `z = a + bi` with `0 <= a < 25`, `1 <= b < 25` and `q` not
dividing `Norm z` (26,264 cases):
* the parity of `log(z/conj z)` equals `[(Norm z/q) = -1]` (`sharp.md` Prop. 2.1(3), used in
  Theorem 3.2 and Definition 4.1);
* `iota` is even iff `q = +-1 mod 8`.

### 2.3 `check_twin.py` and `check_twist.py`

* **Theorem 2.1** (40 random twin profiles, `M <= 6`, `r <= 18`, random moduli
  `m in {8, ..., 72}`, 30 label vectors each):
  * item 1: the inverse map;
  * item 2: `kappa = F + Lambda`;
  * identity (2.2);
  * item 4: the colouring identity `sum_j a_yj nu_j = F_y + nu_(f(y))`, and realisation of an
    arbitrary kappa of that colouring from arbitrary parity-respecting `(ell_J, U)`;
  * unimodularity of `v -> vt`;
  * the norm bounds, exhaustively over `[-2,2]^r` for `r <= 7`.
* **Lemma 4.2 / Theorem 4.3 at all levels**, with `X = 200`, `T(X) = 4`, five twisted columns
  (`zeta_j` over the first five split primes above 200, genuine `p_j ~ 10^6`), and all
  `v in [-2,2]^5 \ {0}`:
  * 20,398 level-collisions, every one equivalent to `q^a | D_s(Z_v rho_(t(q))^(d_q(v)))`;
  * `D_s != 0` throughout;
  * the twisted parities equal `nu_j(q)`;
  * the cap (4.2) is never violated (largest `log m_(v,s)/cap = 0.14`).

### 2.4 `check_misc.py`

* **G** (Proposition 3.1, scope). Random primes:
  * `r = 20`, `p ~ 10^6`, `X = 300`: 28 moduli `q = +-1 mod 8`, all with distinct `nu(q)`. At most
    1 modulus is compatible with any single additive `lambda`. No pair `j < k` has equal symbols at
    all of them.
  * `r = 60`, `p ~ 10^7`, `X = 1000`: 80 moduli, all distinct. Again at most 1 modulus is compatible.
  * The log-product over `q = +-3 mod 8` is `146.2` of `276.3` and `491.8` of `955.6`, about one half.
    This gives the constant `1/8` in place of `1/4`.
* **H** (Theorem 3.2, evidence). For 231,782 pairs `(p, zeta)` with `p < 3000` and `Norm zeta < p`, the
  least prime `q = +-1 mod 8` with `q` not dividing `n` and `(n/q) = -1` is at most `2.70 (log n)^2`.
* **I** (Corollary 4.4). Solving `(2Tk + T + 4)(h log M + 0.7) + T(k+1) lam = tau` for `k` gives (4.3)
  exactly (sympy).

The author's `check_theorems.py` was re-run separately and prints `ALL OK`. Its Theorem 4.3 check is
level 1 only, which `check_twist.py` now extends.

---------------------------------------------------------------------------------------------

## 3. The evidence table (§6) recomputed

**The bug.** `designs.py` computes the "per-unit demand" as follows:

```text
hit = (V % (m/4) == 0);  cls = V // (m/4);  demand_u = sum_q [hit and cls == u] log q
```

The class `u` is the unit `i^s` only when `iota = m/4`. When `iota = 3m/4`, class `1` is `-i`. The
author's generators give `iota = 3m/4` at 13 of the 34 odd `q < 150`. So the reported maximum over `u`
is the maximum over *mixtures* of the units `i` and `-i`, not over units.

**Proof that the published numbers are not per-unit demands.** For genuine labels and `v = e_j`, the
per-unit demand is at most `log |D_s(pi_j)| <= log sqrt(2 p_j)`. On the author's own prime sets
(regenerated with the same seed by `gen_bound.py`):

| `M` | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|
| author's `GEN`, `k = 1` | 2.51 | 2.60 | 2.78 | 3.03 | 3.00 |
| hard cap, `log sqrt(2 * 2M^3)/log M` | 1.68 | 1.65 | 1.63 | 1.62 | 1.60 |
| true maximum (direct divisibility) | 1.66 | 1.65 | 1.63 | 1.61 | 1.60 |

**Reproduction.** `designs_ref.py` uses the referee's own discrete logs, but the same generator
convention, the same prime sets and the same coherent and twisted primes. In "raw-class" mode it
reproduces the author's `COH`, `TWIST` and `GEN` entries exactly at `k = 1, 2`. Examples:
* `M = 50`: `COH 1.63, 2.93`; `GEN 2.51, 4.56`; `TWIST 2.37, 5.04`.
* `M = 800`: `COH 1.45, 2.90`; `GEN 3.00, 6.17`; `TWIST 2.20, 6.19`.

`RAND` differs, because it is a fresh random draw.

**Corrected table** (evidence; unit `s` defined by `<v, ell> = s iota`). Setup as in the draft: `r = M`
columns; `p_j ~ M^3`; labels at all odd `q <= 3M`; demand counted over `q >= M^(1/2)` at level 1.
Entries are `max_v max_s demand/log M` over **all** monomials of norm `k`. (The author's `k = 3` used
only `e_j +- e_a +- e_b` with distinct indices.) The height budget is `w(v)/2 ~ (1.5 + o(1)) k log M`.

| `M` | RAND `k = 1,2,3` | COH (not in class) | TWIST (in class) | GEN (in class) |
|---|---|---|---|---|
| 50  | 2.82, 5.20, 6.43 | 0.92, 1.79, 2.61 | 1.85, 5.04, 6.77 | 1.66, 3.25, 4.80 |
| 100 | 4.56, 5.64, 8.43 | 0.86, 1.69, 2.49 | 2.24, 4.77, 6.87 | 1.65, 3.22, 4.78 |
| 200 | 3.56, 5.99, 8.49 | 0.83, 1.61, 2.38 | 2.17, 5.28, 7.91 | 1.63, 3.19, 4.75 |
| 400 | 4.84, 6.06, - | 0.79, 1.56, - | 2.35, 5.39, - | 1.61, 3.17, - |
| 800 | 4.94, 6.80, - | 0.78, 1.51, - | 2.67, 5.81, - | 1.60, 3.15, - |

**Corrected reading** (evidence only; the sizes are far from asymptotic).
* *Genuine labels* sit exactly at the genuine cap `~1.6 k log M`. This is zero slack, as Theorem 3.2
  and `p3G.md` predict. They are far below the random labels at every `k`, not "with the random ones".
* *Coherent labels* (norm `~X`) sit at `~0.8 k log M`, in line with their certificate
  `log(sqrt2 |Z_v|) ~ (k/2) log(cX)`. They are not at `1.5 k log M`.
* *Twisted labels* already have 2-3.4 times the coherent demand at `k = 1`, and the ratio grows with
  `M` (`TWIST` at `k = 1` goes from `1.85` to `2.67 log M`), as a growing `T(X)` would make it. They exceed the genuine labels at every `k`.
  At `k = 2, 3` they are near the random level. So "twisted labels match coherent ones at `k = 1`" is
  withdrawn. The `k >= 2` part of Remark 4.5's evidence survives.
* None of this enters a proof.

---------------------------------------------------------------------------------------------

## 4. Classes and validity for actual circles

| statement | class claimed | inputs valid for actual circles? | referee |
|---|---|---|---|
| Thm 1.1 | factored data, actual | yes: actual residues factor through `pi_j`, with parity `nu_j(q)` (`sharp.md` Prop. 2.1; R§2.2) | correct; heading of item 2 false (R§1.1) |
| Thm 2.1, Cor. 2.2 | twin profiles in P3-all\* | identities, valid for any twin profile (actual clusters need not have twins) | correct; a localisation of (M1), not a refutation |
| Prop. 3.1 | P3-all\*, `p_j > X` | uses only the per-unit height condition, which is necessary (`sharp.md` Prop. 6.1) | correct; at `q = +-1 mod 8` the hypothesis is unsatisfiable for generic primes |
| Thm 3.2, Cor. 3.3 | P3-all\* (GRH) | a statement about designs; genuine labels satisfy (3.3) with `zeta_j = pi_j` and the conclusion `Norm zeta = p` | correct; Cor. 3.3 needs top-level coherence and non-degeneracy |
| Def. 4.1 - Cor. 4.4 | P3-all\* | statements about one design family | correct; checked at all levels |
| Lemma 5.1 | P3-all\* | single modulus | correct; nearly vacuous for the question |
| "(M2) cannot work" | -- | -- | heuristic, not proved |

**The class definition** (issue 6). The Conventions define P3-all\* as `F^all_X` "with the per-unit
height condition ... and the completion requirement for `v` outside `L` (`sharp.md` §6.1)", and then
state "Actual clusters are members".
* `F^all_X` (`sharp.md` Def. 1.1) already contains the slab conditions (1.2) for every `v`. Those are
  necessary, and the per-unit height condition follows from them.
* The "completion requirement" of `sharp.md` §6.1 is the *union-bound* condition
  `sum_(v not in L) min(1, 4 m_v e^(-w(v)/2))` small. That condition is sufficient for the existence of
  angles, not necessary, and actual clusters violate it: the genuine sums on `p3G.md`'s verified prime
  sets are 4052 and 15190.
* With the completion requirement included, the class does not contain all actual clusters.

No proved result of the draft uses the completion requirement as a necessary condition:
* Proposition 3.1 uses only heights;
* Corollaries 3.3 and 4.4 use the union bound only to say what *certifies*;
* Remark 4.5 and the §5.2 heuristics use the height condition.

So the repair is purely definitional. Define P3-all\* as `F^all_X`, and call the union-bound sum a
construction tool. Every sentence of the form "completion fails" then means "the union bound fails",
which is not an obstruction to membership.

**What the draft says about actual circles.** Nothing beyond `two_thirds`. The draft says so, correctly.

---------------------------------------------------------------------------------------------

## 5. Novelty and attribution

* **Theorem 1.1(1).** The global form of `sharp.md` Remark 2.3(1), which the draft acknowledges. The
  "two elements of each parity" fact at `q = +-3 mod 8` is also used in `sharp.md` Prop. 6.5.
* **Theorem 2.1 / Corollary 2.2.** Essentially `sharp.md` Lemma 6.3(2): the fibre is parametrised by
  the non-twin labels, `u` and `t`, with the twin labels solved. Add Lemma 6.2 (`c'`) and the last
  paragraph of §6.6. The draft cites Cor. 2.2(5) and Lemma 6.2, but not Lemma 6.3(2) or §6.6
  (issue 7).
* **Proposition 3.1.** Acknowledged as `sharp.md` Prop. 7.1 transposed. It is Siegel's-lemma-type
  pigeonhole.
* **Theorem 3.2.** Standard GRH/Chebotarev, in the form used for lower bounds on pseudosquares and
  least non-residues. The specific restriction to `q = +-1 mod 8` (where the `mu_4` shift cannot
  repair parity) is the draft's own (issue 14).
* **Definition 4.1 - Corollary 4.4.** New as far as the referee can tell.
* **§5.1.** The statement "counting gives averages, growth needs concentration" is `sharp.md` §6.7
  ("Pigeonhole forces only average demands"). The draft presents it as a proved failure (issue 3).

---------------------------------------------------------------------------------------------

## 6. Corrected statement

> **`p3Dlower.md` (round 7), as corrected.** P3-all\* is `F^all_X(M, C)`: the slab conditions (1.2)
> for every `v`, and hence the per-unit height condition `m_(v,s) <= nu_s^(-1) e^(w(v)/2)`. The
> union-bound completion sum is a construction tool, not a membership condition. Actual clusters are
> members.
>
> No lower bound `W >= f(M) M log M` with `f -> infinity` is proved, for P3-all\* or for any natural
> subclass. The bracket `3 M log M - O(M) <= W_min^all <= (66 A + o(1)) M log^2 M/loglog M`
> (`X = 3M^A`) is unchanged.
>
> *Proved unconditionally.*
> 1. **(Thm 1.1, corrected heading)** Let `q = +-3 mod 8`.
>    * Every unconstrained label vector becomes parity-respecting after unit shifts `iota u`.
>    * Collision sets and all-unit moduli `m_v` are unchanged. Units are permuted by `<v, u(q)>`.
>    * The parity fixes the parity of the collision unit, `s = <v, nu(q)> mod 2`.
>
>    Let `q = +-1 mod 8`.
>    * A collision of `v` requires `<v, nu(q)>` even.
>    * For `gcd(v, m_q/2) = 1`, a uniform parity-respecting label collides with probability `8/m_q`
>      or `0`, against `4/m_q` without the parity.
>    * The parity also *forces* torsion-type collisions, such as `(m_q/8) e_j` whenever
>      `(p_j/q) = 1`. These are harmless for the height condition.
> 2. **(Thm 2.1, a restatement of `sharp.md` Lemma 6.3(2) in the basis `(ell_J, U, Lambda)`)** On twin
>    profiles:
>    * point and pair data depend only on `kappa = F(ell_J, U) + Lambda`;
>    * every parity-compatible `kappa` is realised together with every parity-respecting `(ell_J, U)`.
>
>    So pair data cannot constrain the `r - 2M` non-twin residues or the twin sums `U`. They act on
>    monomials only through the `M` twin labels `Lambda = kappa - F`. (M1) is *localised* to
>    (M1'), not refuted.
> 3. **(Prop. 3.1)** Additive labels `ell(q) = g_q lambda + iota u(q)` (`|lambda| <= H`) on a set
>    `𝒬` force `W >= (r/(8k))(log Pi_𝒬 - log 2)` whenever `binom(ceil(r/2), k) > 2kH + 1`.
>    * At `q = +-1 mod 8` such labels are parity-admissible only for exceptional prime sets. Under
>      GRH, never once `X >> log^2 p_max`.
>    * So the admissible version has `𝒬 subset {q = +-3 mod 8}`. It gives
>      `W >= (1/8 - o(1)) M X/(h+2)` for `H <= M^h` and `p_j > X`.
> 4. **(Def. 4.1 - Cor. 4.4)** Twisted coherent labels are parity-respecting residue data with the
>    explicit cap (4.2), quadratic in `k`. The cap was verified at all levels. The cap certifies the
>    height condition for norms `k <= k_1 = O(tau/(T(X) log M))`, and does not certify norms
>    `k >> k_1` on columns of weight `~tau`.
> 5. **(Lemma 5.1)** At one modulus, fewer than `m_q/8` monomials with `gcd(v, m_q/2) = 1` can all be
>    avoided by a parity-respecting label.
>
> *Proved under GRH.* **(Thm 3.2)** A label coherent with `log(zeta/conj zeta)` (up to `mu_4`) at
> every `q = +-1 mod 8 <= X`, with the forced parity, has `p Norm(zeta)` a square or twice a square,
> once `X >= c_1 (log 2n)^2 (loglog 3n)^2`. Hence `Norm zeta >= p`. Suppose moreover that the
> coherence holds at top level and the `zeta_j` are non-degenerate. Then the crude certificate
> `m_(v,s) <= |D_s(Z_v)|` has zero slack, and it cannot prove the union-bound completion. This does
> **not** say that such designs, or genuine labels, fail to be members at `O(M log M)`. That is the
> open (I3)/`H` question.
>
> *Heuristic, not proved.*
> * (M2): counting arguments give only average demands (`sharp.md` §6.7). No proof is given that
>   they cannot force concentration.
> * Every analysed design needs `tau >~ log^2 M/loglog M`. The status by design:
>   * random fibre: proved, with the constant between `2` and `2A`;
>   * additive: proved, and much worse;
>   * exact coherence: GRH, but only at the level of certificates;
>   * twisted coherence and `D_A` with coherent free labels: heuristic;
>   * genuine labels: random prime sets / non-adaptive model only (`referee_p3G.md`).
> * The roadmap (I1)-(I3) needs a branch for approximately certified (twisted) designs. (I2) is
>   false for them.
>
> *Evidence.* The §6 table is withdrawn: it was computed with mixed units. Corrected (R§3):
> * genuine labels: about `1.6 k log M`, at the genuine cap (zero slack);
> * coherent labels: about `0.8 k log M`;
> * twisted labels: 2-3.4 times the coherent value at `k = 1` (above genuine), near the random level at
>   `k = 2, 3`;
> * random labels: the largest values.
>
> No rate is read off.

---------------------------------------------------------------------------------------------

## 7. Files (`round7/referee_p3Dlower_checks/`)

| file | content | output |
|---|---|---|
| `tq.py` | referee's exact arithmetic and discrete logs in `T_q`, unit map `s` via `iota`, `D_s` | self-test in text (23,424 cases, 0 mismatches) |
| `check_parity.py` | Thm 1.1(1),(2) exact; forced unit parity (R1); counterexample to "only as a veto"; 553 genuine primes at `q = 7` | `check_parity_out.txt`: ALL OK |
| `check_legendre.py` | parity of `log(z/conj z)` vs Legendre symbol; `iota` even iff `q = +-1 mod 8` | `check_legendre_out.txt`: ALL OK |
| `check_twin.py` | Thm 2.1 items 1-4, (2.2), unimodularity, norm bounds, decoupling | `check_twin_out.txt`: ALL OK |
| `check_twist.py` | Lemma 4.2 at **all** levels with exact `T_(q^A)` logs; cap (4.2) | `check_twist_out.txt`: ALL OK |
| `check_misc.py` | Prop. 3.1 parity-admissibility; Thm 3.2 evidence; Cor. 4.4 algebra | `check_misc_out.txt` |
| `gen_bound.py` | genuine per-unit demand of `e_j` vs the hard cap `log sqrt(2 p_j)` on the author's prime sets | `gen_bound_out.txt` |
| `designs_ref.py`, `designs_ref2.py` | corrected evidence table (correct units) and the author's raw-class bookkeeping for comparison (`ref2`: chunked, `k <= 2`) | `designs_ref_out.txt`, `designs_ref2_out.txt` |
| `notes.md` | working notes (not part of the report) | -- |
