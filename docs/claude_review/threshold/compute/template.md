# Exact isotypic computation of shell vanishing orders (tau-type data), M = 3..9

Task: for `M = 3..9` and the shells `S_{s,t} = {k in Z^M : sum k = s, |k|_1 = t}`, compute the
maximal order `ord(c)` of a nonzero measure `c` on the shell (largest `m` with
`sum_k c_k P(k) = 0` for all polynomials `P` of degree `< m`), over ALL measures, decomposed into
`S_M`-isotypic components. Report the best ratio `ord/t`, the components attaining it, whether
anything beats the Vandermonde ratio or 2, and the trend in `t`. Problem statement:
`../tau_problem.md`.

* Scripts and raw data: `scratchpad/tau/compute/` (§8).
* Full per-shell tables: `../compute_tables_shell.md` (shells) and `../compute_tables_slice.md`
  (full slices).

Conventions.

* *Theorem / Lemma / Proposition / Corollary* means proved here, or cited with the place where
  the proof is written. *Data* means a finite exact-modular computation.
* Modular ranks give **rigorous upper bounds** on `ord` (rank mod p <= rank over Q). Lower
  bounds are rigorous only where an exact certificate or an explicit construction is cited (§5).
* `V(M) := 2 - 1/ceil(M/2)` is the Vandermonde ratio: `3/2, 3/2, 5/3, 5/3, 7/4, 7/4, 9/5` for
  `M = 3..9`.
* Partitions are written in exponent notation: `21^3 = (2,1,1,1)`, `1^M` = sign, `M` = trivial.

## 0. Verdict

**No computed shell beats the Vandermonde ratio, and none comes near 2.**

* **Coverage.** Every shell `S_{s,t}` with `0 <= s <= t` in these ranges was computed:
  `M = 3: t <= 30`; `M = 4: t <= 26`; `M = 5: t <= 24` (plus `t = 25..28` for `s <= 15`);
  `M = 6: t <= @T6@`; `M = 7: t <= @T7@`; `M = 8: t <= @T8@`; `M = 9: t <= @T9@`, plus the Vandermonde
  shell `(p,n) = (10,10)`, `t = 20`, `s = 0`.
* **Size.** That is {{CMD:stats.py nshell}} shells. Every shell was computed with all `S_M`
  components and two primes `P1 = 2^31-1`, `P2 = 2^31-19`, and the two primes give identical
  per-component rank profiles in every case.
* **The maximum ratio is V(M).** For every `M`, `max reg(S_{s,t})/t = V(M)` over the computed
  range, with @M9CAVEAT@. It is attained at the Vandermonde shell by the sign component alone,
  with `reg = C(M,2)`; every other component there is at most `C(M,2) - 1`. There are also a
  few ties (§3.1). No case exceeds `V(M)`, so no case exceeds 2.
* **Theorem 4.4 (computer-assisted).** `tau_sh(3) = tau_sh(4) = 3/2` and `tau_sh(5) = 5/3`,
  where `tau_sh(M)` is the supremum of `ord/t` over ALL shells and all `t`. The proof combines a
  new elementary bound, Proposition 4.3, which reduces each of these `M` to finitely many shells,
  with the data.
* **Trend (Proposition 4.3).** At fixed `M`, `ord/t` decreases beyond the Vandermonde shell and
  tends to 1, because `t <= reg(S^M(p,n)) <= t + 2^M - 3` for `M >= 4` and `p, n >= 1`. In the
  data the excess `max_s reg - t` is eventually periodic for `M = 3, 4, 5`, with maximum
  `1, 3, 5` respectively, and still grows slowly for `M >= 6`. So the shell quantity is
  maximised at small `t`, near the Vandermonde shell.
* **Full slices.** These carry the invariant that the reduction actually needs. Data on
  {{CMD:stats.py nslice}} slices (`M <= 9`, `t <= 24/20/18/15/13/13/12`, plus the `M = 8`
  Vandermonde slice): `reg Q^M(p,n) <= floor(phi_M(p,n))` holds in every case. This is an
  independent check of `../upper.md` Theorem 3.1 far beyond its `M <= 6, t <= 7` range.
  * `reg = C(M,2)` at the Vandermonde slices for `M = 3..8` (`upper.md` Theorem 4.4(a)).
  * Unlike shells, slices keep the ratio `V(M)` for arbitrarily large `t`. They reach it exactly
    at multiples of the Vandermonde point, as products allow.

## 1. What is computed

For a finite set `A` in a hyperplane `{sum k = s}` let
`reg(A) = min{d : every function A -> Q is the restriction of a polynomial of degree <= d}`.

* By duality (Lemma 2.1), `reg(A) = max{ord(c) : c != 0 on A}`.
* So the requested quantity for a shell is `reg(S_{s,t})`, and the ratio is `reg(S_{s,t})/t`.
* Since `k -> -k` maps `S_{s,t}` onto `S_{-s,t}`, only `s >= 0` is computed.
* Shells are written `S^M(p,n)`, with `p = (t+s)/2` (positive mass) and `n = (t-s)/2`
  (negative mass).

Two refereed notes (`../verify.md` §3.2 and `../upper.md` Remark 1.3) showed that the invariant
relevant for the reduction is the full slice `Q^M(p,n) = union_i S^M(p-i, n-i)`, not the shell,
and that `tau_sh <= tau`. The task asks for shells. §6 adds slices, computed with the same code.

METHOD_SECTION

### 2.4 Validation

* **Direct cross-check.** `sum_lam f^lam h_lam(d)` was compared, degree by degree, with the
  Hilbert function of the whole shell, computed directly from all monomials in `k_2..k_M`
  (`direct.py`, mod `P1`). The result: {{CMD:stats.py crosscheck}}. This covers `M = 2..5`
  with `t <= 8`, `M = 6` with `t <= 8`, `M = 7` with `t <= 7` and `M = 8` with `t <= 6`, all `s`.
  * Each `h_lam^{gen}` is at most `h_lam`, so equality of the sums forces
    `h_lam^{gen} = h_lam` (mod `P1`) for every `lam`.
  * So in these cases the generators are complete, independently of the ATY theorem.
* **Previous independent computations.** These are the direct, non-equivariant shell
  regularities of `../verify.md` (E1: `M = 3, 4, 5`) and `../upper.md` (`upper_checks/log_M*.txt`:
  `M = 3..6`). All @PREVN@ overlapping shells agree (`compare_prev.py`).
* **Two primes.** Every shell and slice was run with `P1` and `P2`. The full per-component rank
  profiles `h_lam(d)` agree in every case.
* **Parity split.** It was checked against the unsplit computation on 20 shells and slices
  (`test_parity.py`).
* **Sizes.** `sum_lam f^lam N_lam` equals the shell cardinality (closed formula) in every run.
* **Columns.** The final rank equals `N_lam` in every run, so the `c_T` are independent.

## 3. Results for shells

### 3.1 Maxima per M

{{TABLE:shell:maxima}}

MAXIMA_NOTES

### 3.2 The Vandermonde shells, component by component

{{TABLE:shell:vand}}

VAND_NOTES

### 3.3 Trend in t: max over s of reg/t

{{TABLE:shell:trendline}}

(For `M = 5` and `t = 25..28` the maximum is over the computed `s <= 15` only. By Corollary 4.4,
the omitted shells cannot exceed `5/3`.)

The excess `max_s reg(S_{s,t}) - t`, for `t = 1, 2, ...`:

{{TABLE:shell:excess}}

TREND_NOTES

### 3.4 Which components attain the maximum

COMPONENT_NOTES

TREND_SECTION

## 5. Exact certificates (lower bounds)

CERT_SECTION

## 6. Full slices Q^M(p,n) (the invariant the reduction needs)

SLICE_SECTION

## 7. Consequences for the problem

CONSEQ_SECTION

## 8. Files

FILES_SECTION
