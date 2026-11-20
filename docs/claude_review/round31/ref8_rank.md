# Referee report on Part VIII.1 (rank recovery at inert primes) and its consequences (VIII.4), (VIII.5), field test

Referee block: rank. Source refereed: `astra_free_VIII/answer.md`, lines 1667–1708 (VIII.1, consequences, field test), with the
ten-mark data claim on lines ~1933. Own code: `ref8_rank/` (listed in section 8). Astra's scripts were not imported.

## 0. Verdict

| Item | Label | Short verdict |
|---|---|---|
| VIII.1 rank theorem (affine independence of enlarged rows; `M <= 1 + omega(N') + K_C`) | **Proved** (every line checked) | Correct as stated. It is also improvable: see Correction 1. |
| "For `C < 1` one may use the empty set" | **Proved, not sharp** | The empty set works for every `C < sqrt 2`, because `(1+i)` always divides `b_ij` (Correction 1). |
| (VIII.4) `Y = P_+(N') >= (2 - o_C(1)) M log M`; `M` bounded when `Y <= (2-eta) M log M` | **Proved** | Correct. The bounded subwindow is tiny compared with the open window `Y <= M^{2-delta}`. |
| (VIII.5) coefficient `2/k` under the exponent-`k` hypothesis | **Proved** | Correct. For an integer radius `R` it gives `1/2` exactly when `R` is squarefull (Correction 2). |
| Field test (transfer to `F_2[t]`, `F_3[t]` constant extensions) | **Proved** (transfer). Family remarks: consistent, partly checked | The theorem is function-field valid. By (F1) it therefore cannot by itself prove the uniform bound. |
| Real-data test | **Computed** | No failure in about 1.42 million windows (exhaustive `N <= 10^7` at `C <= 3`, exhaustive `N <= 10^6` at `C <= 9`, structured `N <= 10^13`), nor on any record cluster or family. |
| Novelty | — | New and unconditional. The record had `omega(N') > M/4` and `omega >= (M/2)(1-o(1))` (both for `C <= 1`), and `M <= omega(N') + O(1)` only conditionally (round 20, (UI-l)). |

**The most serious issue.** No step of VIII.1 is wrong. The issues are scope and sharpness. The constant `K_C` of (VIII.2)
grows like `C^4` up to a logarithmic factor: `K_C = 4168` at the record's 10-point constant `C = 8.24`, while the observed
defect `M - 1 - omega(N')` there is 4. The statement "for `C < 1` use the empty set" undersells the theorem: the threshold
is `sqrt 2`. Also, for `C < sqrt 2` the proof uses only the pair inequality, so angle-model fakes satisfy the same rank bound.
The theorem constrains integer circles and fakes alike in that range, and it is valid verbatim over function fields.

## 1. Setting and definitions

* `N = R^2` is a positive integer; `R = sqrt N` need not be an integer (record convention, `expert_note.md` §1).
* The marks `z_1, ..., z_M` are distinct elements of `Z[i]` with `|z_i|^2 = N`, lying on an arc of length `<= C sqrt R`.
* Let `g = gcd(z_1, ..., z_M)` in `Z[i]`. Put `z_i' = z_i/g` and `N' = N/|g|^2 = prod p^{e_p}`, `W = log N'`,
  `r = omega(N')`. For each `p | N'` fix `pi_p | p` and put `a_ip = v_{pi_p}(z_i')`.
* Write `g_ij = gcd(z_i', z_j')`, `b_ij = (z_i' - z_j')/g_ij` and `d_ij = sum_p |a_ip - a_jp| log p`.
* Let `S` be a finite set of primes `l = 3 mod 4`. The *enlarged row* of mark `i` is `a_i` followed, for each `l` in `S`, by the
  one-hot indicator vector of the class of `z_i' mod l`, taken in the norm fibre `{x in F_{l^2} : x^{l+1} = N' mod l}`.

## 2. Line-by-line check of VIII.1 (Proved)

1. **Normalisation and arc constant.** The normalised marks lie on `|z| = R' = R/|g|`. Every chord satisfies
   `|z_i' - z_j'| <= C sqrt R/|g| = (C/sqrt|g|) sqrt R' <= C N'^{1/4}`. The new constant is `C/sqrt|g| <= C`. Correct.
   Only chords are used later, so the theorem holds with `C` replaced by the chord constant `C' = max|z_i'-z_j'|/N'^{1/4}`.
   My code uses this sharper form.
2. **Only split primes remain.** Suppose an inert `l` divides `N'`. Then `l | z_i' conj(z_i')`, and since `l` is prime in
   `Z[i]` it divides every `z_i'`, contradicting `gcd = 1`. If `2 | N'`, then `(1+i)` divides every `z_i'`, again a
   contradiction. Hence `N'` is odd and every `p | N'` splits. Correct.
3. **Columns.** The normalisation gives `min_i a_ip = 0` and `max_i a_ip = e_p`, with `e_p >= 1`, so no column is
   constant. If `e_p = 0` after normalisation, `p` is simply not a coordinate. This is consistent.
4. **`|g_ij|^2 = exp(W - d_ij)`.** The valuations of `g_ij` are `v_pi = min(a_ip, a_jp)` and
   `v_{conj pi} = e_p - max(a_ip, a_jp)`; their sum is `e_p - |a_ip - a_jp|`. There are no other primes, and units do not
   change moduli. So `log|g_ij| = (W - d_ij)/2`. Correct, including units and the prime 2, which is absent.
5. **(VIII.3).** Item 1 gives `|g_ij| |b_ij| <= C e^{W/4}`, hence `d_ij >= W/2 + 2 log|b_ij| - 2 log C`. Correct.
   This needs `b_ij != 0`, that is, distinct marks; the hypothesis "distinct" should be stated. Marks differing by a unit
   are allowed: then `a_i = a_j`, `d_ij = 0`, and (VIII.3) is still true because `|b_ij| = |u - 1| >= sqrt 2` and the
   chord bound forces `R` to be tiny.
6. **Affine dependence to two measures.** Take `sum c_i x_i = 0` with `sum c_i = 0` and `c != 0`. Normalise the positive
   and negative parts to probability vectors `P` and `Q`. They have disjoint nonempty supports (because `sum c = 0`) and
   equal barycentres in every coordinate. In the indicator coordinates, equal barycentres means the residue
   distributions of `P` and `Q` mod `l` coincide. Correct.
7. **Pointwise inequality.** For `x, y` in `[0,1]`, `x + y - 2xy - |x - y| = 2 min(x,y)(1 - max(x,y)) >= 0`. Take `i ~ P` and
   `j ~ Q` independent; then `E[XY] = mu^2`, and `E|X - Y| <= 2 mu (1 - mu) <= 1/2`. Summing with weights `e_p log p` gives
   `E d_ij <= W/2`. Correct. Independence is a choice made in the proof, not an assumption.
8. **Norm fibre.** `Z[i]/l = F_{l^2}`, and conjugation acts as Frobenius, so `N(z mod l) = (z mod l)^{l+1} = N' mod l`.
   This value is nonzero because `l` does not divide `N'` (item 2). The norm map `F_{l^2}^* -> F_l^*` is a surjective map
   of cyclic groups, so each fibre has `(l^2 - 1)/(l - 1) = l + 1` elements. Correct.
9. **Collision probability.** `P(z_i' = z_j' mod l) = sum_t P_l(t) Q_l(t) = sum_t P_l(t)^2 >= 1/(l+1)` by Cauchy–Schwarz
   over the `l + 1` classes. Correct.
10. **Simultaneous divisibility.** A collision means `l | g_ij b_ij`, and `l` does not divide `g_ij` because `|g_ij|^2`
    divides `N'`. So `l | b_ij`. Distinct inert primes are pairwise coprime in `Z[i]`, so the product of the colliding
    `l` divides `b_ij`, and `log|b_ij| >= sum_{l colliding} log l` (the cofactor is a nonzero Gaussian integer). Taking
    expectations, `E log|b_ij| >= sum_S log l/(l+1)`. Correct.
11. **Strictness.** The two bounds give `W/2 >= E d_ij >= W/2 + 2(E log|b| - log C)`, and `E log|b| > log C` makes this
    impossible. The inequality in (VIII.2) is strict, so the contradiction is strict. If `C < 1`, then `log C < 0` and
    `S` may be empty. Correct.
12. **Dimension count.** The indicator block for `l` lies in the hyperplane "coordinates sum to 1", so its affine rank is at
    most `l`. It is also at most the number of residue classes hit, minus 1. Affine rank is subadditive over a
    concatenation of coordinate blocks. So `M - 1 = affrank(enlarged) <= affrank(a) + sum_S l <= r + K_C`. This gives
    `M <= 1 + omega(N') + K_C`. Correct.

**Correction 1 (sharper, Proved).** `N'` is odd, so every `z_i'` is congruent to a unit, hence to 1, mod `(1+i)`. Therefore
`(1+i) | z_i' - z_j'`, and since `(1+i)` does not divide `g_ij`, also `(1+i) | b_ij`. Combined with the inert part this gives
`log|b_ij| >= (1/2) log 2 + sum_{l colliding} log l`. The record already uses `|beta| >= sqrt 2`
(`round8/stability.md` Prop 3.4 and 4.1).

Consequently the theorem holds with (VIII.2) replaced by
`sum_S log l/(l+1) > log C - (1/2) log 2`. In particular:

* **For `C < sqrt 2` (with `C` the chord constant), the exponent rows of every normalised cluster are affinely independent,
  so `M <= omega(N') + 1`.**
* Adding the one-dimensional block "`z' mod 2`" (residues `1, i`) gives `E v_{1+i}(b) >= 3/2`. Then `M <= omega(N') + 2` for
  `C < 2^{3/4} = 1.682`.
* More generally, blocks mod `2^j` give `E v_{1+i}(b) >= 2 - 2^{1-2j}`.

Optimal `K_C` (sum of `l`), computed exactly (`kc.out`). The optimum by branch and bound equals the initial segment
`3, 7, 11, 19, ...` wherever the search finished.

| `C` | 1 | 1.3 | `sqrt2` | 1.7 | 2 | 3 | 4 | 5.9 | 8.24 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| `K_C`, Astra's (VIII.2) | 3 | 3 | 10 | 21 | 21 | 94 | 310 | 1301 | 4168 | 8684 |
| `K_C` with the free factor `(1+i)` | 0 | 0 | 3 | 3 | 10 | 40 | 94 | 381 | 1150 | 2411 |

Asymptotics (Computed and standard): `sum_{l <= x, l = 3 (4)} log l/(l+1) = (1/2) log x - 0.741 + o(1)`, computed to
`x = 4*10^5`. Hence the admissible `x` is about `4.4 C^2`, and `K_C ~ x^2/(4 log x) = Theta(C^4/log C)`.

## 3. Test on real data (Computed)

Method. For each `N`, all lattice points are generated exactly from Gaussian factorisations. I took every window of
consecutive points with arc constant `C = arc/N^{1/4} <= C_max`, with the start angle in `[0, pi/2)`. Each window was
normalised exactly, and for each I computed:

* `omega(N')`;
* the exact integer affine rank of the exponent rows;
* the chord constant `C'`;
* the enlarged rank with Astra's `S` (initial segment valid for `C'`), and with the improved `S` of Correction 1.

The script **asserts** the theorem in every window: enlarged rank `= M - 1`, `M <= 1 + omega + K_C`, and rank `= M - 1`
whenever `C' < sqrt 2`. As a negative control, it also tests `S` minus its last prime.

| Range | Windows | Failures |
|---|---|---|
| all `N <= 10^7` with every prime `= 1 mod 4`, `C <= 3` (814,182 circles) | 535,302 | 0 |
| all such `N <= 10^6`, `C <= 8` and `C <= 9` (87,881 circles) | 272,152 / 115,424 (`M >= 3`) | 0 |
| structured `N <= 10^13` built from the first 30 split primes, `r_2(N)/4 >= 384` (3,358 circles), `C <= 9` | 494,667 | 0 |

The negative control (one prime short of the condition) **never** produced a dependence. The residue blocks are far more
effective than the expectation bound requires.

**Largest observed defect `M - 1 - omega(N')`, by arc-constant threshold `T`** (all circles with `N >= 10^4`; tiny circles
such as `N = 25`, whose whole circle has small chord constant, give larger defects: 6 at `C = 7.66`, `N = 25`):

| `T` | `<= 2.5` | 3 | 4 | 5 | 5.9 | 8 | 8.3 | 9 |
|---|---|---|---|---|---|---|---|---|
| max `M - 1 - omega` | **0** | 1 | 1 | 2 | 3 | 4 | 4 (record 10-pt: 4) | 5 |
| witness `(N, M, omega)` | — | (707281, 3, 1) | same | (21025, 5, 2) | (21125, 6, 2) | (60125, 8, 3) | (1176852625, 10, 5) | (27625, 9, 3) |
| Astra's `K_C` | ≤63 | 94 | 310 | 753 | 1301 | 3885 | 4168 | 5823 |

* The first defect (rank or `omega`) occurs at arc constant `C = 2.8287` (`N >= 10^4`); in chord constant the minimum is
  `C' = 2.4903` (`N = 21125`). Over all `N` the minimum chord constant of a defect is exactly 2.0 (`N = 25, 125`).
* **No rank deficiency occurs below `C' = 2` in any scanned window**, whereas the proof gives independence only for `C' < sqrt 2`.
  Whether `C' < 2` forces affine independence without residue blocks is open (**Heuristic**, data only).

**Record clusters and families** (`families.out`; exact):

| cluster | `M` | `omega(N')` | aff. rank | `C'` | +mod 3 | +mod 7 | `z mod 2` block | Astra `K_C` |
|---|---|---|---|---|---|---|---|---|
| 10-pt, `N = 1176852625 = 5^3 13^2 17 29 113` | 10 | 5 | 5 | 8.2386 | 6 | **9 = M-1** | 5 | 4168 |
| 9-pt sub-cluster | 9 | 5 | 5 | 7.8798 | 6 | 8 | 5 | 3614 |
| 8-pt, `N = 32988780417125` | 8 | 6 | 4 | 5.6739 | 5 | 7 | 5 | 1150 |
| 7-pt, `N = 235370525` | 7 | 5 | 5 | 5.2695 | 6 | — | 5 | 880 |
| golden `n = 4, 7, 10, 13` (gcd `(1+i)^2` removed) | 8 | 5, 7, 8, 9 | **4** (all `n`) | 5.8977 | 5 | 7 | 5 | 1301 |
| axis Pell, 5 orbit members `a = 685549 ... 9.1*10^11` | 6 | 5, 6, 7, 8, 9 | **3** (all) | 4.0000 | 5 | — | 4 | 310 |
| 5-pt `C <= 3`, `N = 10414625` | 5 | 4 | 4 = M-1 | 2.4168 | | | | 63 |
| 4-pt `C <= 1`, `N = 412333805` | 4 | 4 | 3 | 0.881 | | | | 0 |

Astra's ten-mark data claim is confirmed: rank 5, then 6 with the mod-3 block, then 9 with mod 7.

In the two infinite families `omega(N')` grows while the affine rank stays fixed at 4 (golden) and 3 (Pell). Their rank
deficiency `M - 1 - rank` stays at 3 and 2 for every member, and in each member the blocks at `l = 3, 7` restore full rank.
The family members are not literal Cilleruelo–Granville constructions. I could not locate the CG six-point example in the
record (the cited `codex_uniform_bound_research.md` is not in the repository); the record says it has `C -> 4 sqrt 2`
after removing `(1+i)^2`.

## 4. Sharpness (part Proved, part Computed, part open)

* **`K_C` must grow with `C` (Computed).** The data give `max(M - 1 - omega) >= 1, 2, 3, 4, 5` at arc constants
  `2.83, 4.00, 5.18, 6.64, 8.75`. Any valid `K_C` is at least these values. The theorem's `K_C ~ C^4/log C` is far
  above them.
* **`M = omega(N') + 1` occurs (Computed),** for example the 5-point cluster on `N = 10414625` at `C' = 2.417`, and many
  pairs and triples. Equality for arbitrarily large `M` at bounded `C` would contradict the conjecture `M(C) < infinity`.
  It is therefore not expected, and none is known. In angle models (fakes) the full-rank profiles (Hadamard / Paley)
  are standard, but they have `r >= M`.
* **Without inert primes (Proved / open).** Below `C' = sqrt 2` no residues are needed (Correction 1). Between `sqrt 2` and
  2 the data show no deficiency, but no proof is known. From `C' = 2` (tiny circles) or `C = 2.83` (`N >= 10^4`) on,
  rank deficiency actually occurs. So **some** additional input is needed above about 2.5, and the residue blocks are
  one sufficient input.

## 5. Consequences

**(VIII.4) — Proved.** Active primes are split and at most `Y`, so `r <= pi(Y; 4, 1)`. Combined with `r >= M - 1 - K_C` and
`pi(y; 4, 1) ~ y/(2 log y)`, this gives `Y >= (2 - o_C(1)) M log M`.

For fixed `C` and `eta > 0`: if `Y <= (2 - eta) M log M`, then `pi(Y; 4, 1) <= (1 - eta/2 + o(1)) M`. This is
incompatible with `M - 1 - K_C` once `M >= M_0(C, eta)`, so the bound is uniform in `R` and depends only on `(C, eta)`.

Comparison with the record: Theorem 3.3(3) (`omega > W'/(4W'/M + 8 log^+ C)`, so `(1/4 - o(1)) M`) gave
`Y >= (1/2 - o(1)) M log M`, and Prop 3.4 (`C <= 1`) gave `(1 - o(1)) M log M`. VIII.4 improves the constant to 2. It does
**not** reach the open window `Y <= M^{2-delta}`: the subwindow it excludes has `log Y ~ log M`, a vanishing part of it.

**(VIII.5) — Proved.** Assume each prime active after normalisation has exponent `>= k` in the original `N`. The active
primes are `r >= M - 1 - K_C` distinct primes, so `log N >= k log((r+1)!) = (k - o_{C,k}(1)) M log M` (using
`p_n >= n + 1`). Using the first `r` primes `= 1 mod 4` instead gives the same leading term `r log r`, so `2/k` is the
correct output of the method. Inverting `x log x` gives `M <= (2/k + o(1)) log R/loglog R`. The hypothesis is on the
original `N`, which is the right one, since normalisation can lower exponents. It improves Theorem 1.1's `2/3` only for
`k >= 4`; it equals it at `k = 3` and is worse at `k = 2`.

**Correction 2 (meaning of "four-full").** In the record `N = R^2` is any integer. If one insists that `R` be an integer, then
`N` is a square and its exponents are even. Then:

* `N` squarefull is automatic, giving `k = 2` and coefficient 1, which is no information;
* "cubefull `N`" and "four-full `N`" both mean that `R` is squarefull, giving `k = 4` and coefficient `1/2`;
* `R` being `k'`-full gives `k = 2k'` and coefficient `1/k'`.

Astra's summary item 8 ("four-full integer norms: `2/3 -> 1/2`") is correct with this reading.

## 6. Field test — Proved (transfer); family remarks consistent

Every step of section 2 is valid in `F_p[t] ⊂ F_{p^2}[t]` with `|f| = p^{deg f}`:

* there is no finite ramification;
* odd-degree irreducibles `lambda` are inert, with norm fibre of size `p^{deg lambda} + 1`;
* collisions give `lambda | b`;
* the pair bound reads `d_ij >= W/2 + 2 log|b| - 2 log C`.

The inert contributions `sum (deg lambda) log p/(p^{deg lambda} + 1)` are about `log p/2` per odd degree, so they diverge.
At critical contact (`log C = 0`) one inert `lambda` suffices, and its block has affine dimension `<= p^{deg lambda}`.

Because the theorem transfers, every function-field family satisfies it automatically. For VI.5 the rows are the `m - 1`
nonconstant Walsh columns on disjoint prime blocks, so the affine rank is `m - 1 = M - 1` (Hadamard nonsingularity), which
is consistent. Astra's orthogonality remark for VII.1 I did not re-derive; it is not needed logically.

Consequence, not stated by Astra: since VIII.1 is function-field valid, and VI.5 has unbounded clusters, the rank bound
cannot by itself yield the uniform bound (F1). For `C < sqrt 2` the proof uses only the pair inequality, so the record's
angle-model fakes satisfy it too.

## 7. Novelty against the record and the literature

* **Record.**
  * Theorem 3.3 (`omega(N') > M/4`, `C <= 1`) and Prop 3.4 (`omega >= (M/2)(1 - M^{-1/3})^2`, `C <= 1`) are in
    `round8/stability.md`.
  * Prop 4.1(b) gives Gram vectors `f_x` in cut space with `<f_x, f_y> = -2 m_xy`. For `C' < sqrt 2` all `m_xy > 0`, so
    the `f_x` are affinely independent in dimension `Omega(N')`. This yields `M <= Omega(N') + 1`, not `omega + 1`: the map
    from exponent rows to cut vectors is not affine.
  * "Rank `>= M - 1`" appears only as the hypothesis (AI) on angle models (`round5/exact.md`).
  * `M <= omega(N') + O_C(1)` appears only as the *ceiling* of a conditional input (UI-l) (`round20/wild.md` Prop 6.8,
    Cor 6.9).
  * Lemma 3.1 "rank `>= n/4 - 2`" is quoted in `round31/README_VII.md` §1.6.

  **VIII.1 adds** two things: (i) affine independence in *exponent* space, which upgrades `Omega` to `omega`, via the
  `2 mu (1 - mu)` barycentre trick; and (ii) the unconditional extension to every `C` through inert residues. It is new
  and correct.
* **Literature.** I ran one web search and read only snippets, not papers. It surfaced:
  * Cilleruelo–Granville, *Close lattice points on circles* (Canad. J. Math. 2009): `R^{1/3}`-scale results;
  * the Cilleruelo–Córdoba bound (at most `l` points on arcs of length `sqrt 2 R^{1/2 - 1/(4[l/2]+2)}`);
  * a Bourgain–Rudnick lemma (`O(log m)` points on arcs of length `m^{1/4}`).

  No bound of the shape `M <= omega(N) + O(1)` appeared in what I saw. I did not read Chang's work. Absence of such a bound
  is therefore not established beyond these snippets.

## 8. Files (all under `scratchpad/round31/ref8_rank/`)

* `gauss.py`: exact Gaussian arithmetic, integer affine rank, residue blocks.
* `scan.py`, `scan2.py`: window scan with the theorem's assertions; `scan2` also handles non-primitive `N`.
* `scan3.py`, `tab.py`, `kcfun.py`: record dump and threshold tables.
* `families.py`, `run_families.py`, output `families.out`.
* `kc.py`, output `kc.out`: `K_C`.
* Outputs: `all_1e6_C8.out`, `all_1e7_C3.out`, `tab_all1e6.out`, `tab_struct13.out`.
