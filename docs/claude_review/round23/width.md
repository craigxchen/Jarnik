# Round 23 — width.md: the weighted mean-width inequality (W) and tau_e(-K_{LM_5})

## 0. Summary and verdict

**Verdict: map item settled; target untouched.** `τ_e(-K_{LM_5}) = 10/3` is now **proved** (split and twisted forms).
The general weighted mean-width inequality (W) is **not** proved; it is no longer needed for `τ_e`.

1. **Theorem 2.5 (proved, elementary, self-contained).** For every `m ≥ 1`, `t_max(m) = max{ord_e f : 0 ≠ f ∈
   H^0(m(-K))} = ⌊10m/3⌋`. Hence `τ_e(-K) = 10/3`. The same proof gives `t_max(m) = 3m` for all `m` on `LM_4` and
   `LM_3 = dP_6`.
   * Mechanism: `S_m = mR_5 ∩ M` is the lattice-point set `T_5(m,m)` of `mΔ - mΔ`; slicing `T_n(a,b)` by a coordinate
     gives exactly the sets `T_{n-1}(a', b)`, `T_{n-1}(a, b')`. The Theorem B peeling lemma then yields a closed recursion,
     solved by explicit functions `F_n(a,b) = a + b + H_n(a,b)` with `H_5(a,b) = min(min(a,b) + ⌊max(a,b)/3⌋, 3 min(a,b))`.
     The one non-obvious step is the convex combination `(1/3)X + (2/3)Y` in Lemma 2.4.
   * It replaces the round-22 *evidence* (`t_max(m) = ⌊10m/3⌋`, `m ≤ 8`, by rank computations) with a proof for all `m`.
2. **(W) on the family `aΔ - bΔ` (proved, Cor 3.1)**, with equality exactly on `S_{3s}`. General (W) is Conjecture
   C_bal of `construct.md` at `M = 5` (Remark 1.1). It is open. New evidence: 150 certified instances (signed Minkowski
   sums of simplices) and 60 pair-box polytopes on which iterated peeling stays below the (W) value. I explain why the
   Theorem B template, run with mean-width slice bounds, cannot prove it (§3.2).
3. **Consequence (proved, with `toric.md` Thm 3.3).** For every finite family of auxiliary Laurent polynomials over
   `Q(i)`, used through Liouville with toric forced divisibility, the best exponent for arbitrary profiles is
   `δ_F ≥ 2/τ_e - 1/2 = 1/10`. That is exactly Cilleruelo–Córdoba, attained by the Vandermonde.
   * So **toric Liouville cannot improve Cilleruelo–Córdoba off a closed set**, for the twisted and split forms alike.
     The round-22 sentence that kill_toric F1 rejected as unsupported is now a theorem.
   * The cap `1/10` is uniform over `F`. The profile that witnesses it may depend on `F` (§5 item 2).
4. **Nothing here touches `A_5(δ < 1/30)`, `α^ess_e > 15/4`, or `(L_5^∘)`.** The result lowers a section cap; it does
   not open a lane. `M ≤ (2/3 + ε) log R/loglog R` remains the best bound. No growth or uniformity claim is made.
5. **No section beats (W) at any `m`.** The refutation branch of the task (a section of `mB` with order `> 10m/3`, for
   some `m ≥ 9`) is **impossible** by Theorem 2.5, so no search was needed.

## 1. Setup and statement of (W)

**Conventions.** *Theorem/Proposition/Lemma* = proved here with all quantifiers (or cited with the place of the
proof). *Exact* = finite exact integer computation; ranks over the prime field `F_p`, `p = 2^61 - 1`, are used **only**
in the direction "full rank mod `p` ⇒ full rank over `Q`" (so they certify upper bounds for `reg`); a rank deficiency mod
`p` is labelled *evidence*. Scripts are in `width_checks/`. Nothing here has been refereed.

**Objects** (as in `round22/toric.md` §1, §3).
* `M = {x ∈ Z^5 : Σ x_j = 0}`; for `A ⊆ [5]`, `x_A := Σ_{j∈A} x_j`. On `M`, `x_{A^c} = -x_A`.
* `B = -K` on `X = L̄_5^tw` (or on the split `LM_5`); `P_B = R_5 = conv(e_a - e_b) = {x : Σx = 0, Σ_j (x_j)_+ ≤ 1}`.
  `H^0(mB)` = Laurent polynomials with support in `S_m := mR_5 ∩ M`.
* For finite `S ⊂ Q^n`: `reg(S)` = least `d` with `Poly_{≤d}|_S = Q^S`. It depends only on `S` up to affine
  isomorphism, is unchanged by extending the field (a rank), and is monotone: `S ⊆ S' ⇒ reg(S) ≤ reg(S')` (restriction
  of a surjection to a subset stays surjective). A set of `N` collinear points has `reg = N - 1`.
* `toric.md` Lemma 3.1 (proved, killer-checked): `max{ord_e f : 0 ≠ f, supp f ⊆ S} = reg(S)`; in particular
  `t_max(m) := max{ord_e f : 0 ≠ f ∈ H^0(mB)} = reg(S_m)`, and `τ_e(B) = lim_m t_max(m)/m = sup_m t_max(m)/m`
  (Fekete; `t_max` is superadditive since orders and Newton polytopes add under products).
* `width_A(S) := max_S x_A - min_S x_A`.

**The inequality (W)** (round 22, judge check J4; `toric.md` §3 Remark (1)):

```text
(W)    reg(S)  <=  (1/6) * Σ_{A ⊂ [5], |A| = 2} width_A(S)        for every finite S ⊂ M.
```

For `S = S_m` every `width_A = 2m`, so (W) gives `t_max(m) ≤ 20m/6 = 10m/3`.

**Remark 1.1 ((W) is Conjecture C_bal at M = 5; derived).** `threshold/construct.md` §2.8 states Conjecture C_bal:
`ord(c)·π_bal ≤ E_bal w_P(supp c)` with `P` a uniform balanced cut and `π_bal = ⌈M/2⌉/(2⌈M/2⌉ - 1)`. At `M = 5`:
`π_bal = 3/5`; the balanced `P` are the 20 sets with `|P| ∈ {2,3}`, and `w_P = w_{P^c}` on a slice, so
`E_bal w_P = (1/10) Σ_{|A|=2} w_A` and C_bal reads `ord ≤ (5/3)(1/10) Σ_{|A|=2} w_A = (1/6) Σ_{|A|=2} w_A`, which is (W)
(with `ord ≤ reg` and Lemma 3.1 the two are equivalent). At `M = 4` C_bal reads `reg ≤ (1/4)Σ_{|A|=2} w_A`
(`= 3m` on `mR_4`). The refereed Theorem B (`construct.md` §2.4) is the uniform-cut version, giving `15m/4`.

**What is needed for the map.** For `τ_e(B) = 10/3` one does **not** need (W) for all `S`: one needs
`reg(S_m) ≤ 10m/3 + o(m)`. This round proves `reg(S_m) ≤ ⌊10m/3⌋` for every `m` (§2) by a slicing argument that is
specific to the polytopes `aΔ - bΔ`, and proves (W) for that whole two-parameter family (§3). The general (W) stays
open (§3.2).

## 2. Theorem: t_max(m) = ⌊10m/3⌋ for every m, hence τ_e(-K_{LM_5}) = 10/3

### 2.1 The sets T_n(a,b) and their slices

For integers `n ≥ 1`, `a, b ≥ 0` put

```text
T_n(a,b) := { y ∈ Z^n : Σ_j y_j = a - b,  Σ_j (y_j)_+ ≤ a }.
```

On `Σ y = a - b` the conditions `Σ y_+ ≤ a` and `Σ y_- ≤ b` are equivalent (`Σy_+ - Σy_- = a - b`), so `T_n(a,b)` is
the set of lattice points of the Minkowski sum `aΔ_{n-1} - bΔ_{n-1}` (`Δ_{n-1} = conv(e_1,…,e_n)`); it is the lattice
of Newton supports of `H^0(aH + bH')` in the notation of `toric.md` Prop 1.1(b). In particular

```text
S_m = mR_5 ∩ M = T_5(m, m).                                                        (2.1)
```

**Lemma 2.1 (slices; proved).** Let `n ≥ 2`. The slice `T_n(a,b) ∩ {y_1 = c}` is nonempty iff `-b ≤ c ≤ a`, and,
after deleting the coordinate `y_1`, it equals `T_{n-1}(a - c, b)` if `0 ≤ c ≤ a`, and `T_{n-1}(a, b - |c|)` if
`-b ≤ c < 0`. Moreover `T_1(a,b) = {a - b}` and `T_2(a,b)` consists of `a + b + 1` collinear points, so
`reg T_1 = 0`, `reg T_2(a,b) = a + b`.

*Proof.* Write `y = (c, z)`. If `c ≥ 0`: `Σ z = a - b - c = (a - c) - b` and `Σ y_+ = c + Σ z_+ ≤ a` iff
`Σ z_+ ≤ a - c`; this is `T_{n-1}(a-c, b)` (nonempty for `c ≤ a`: every `T_{n'}(a',b')` with `a', b' ≥ 0` contains
`a'e_1 - b'e_2`, or `(a'-b')e_1` if `n' = 1`). If `c < 0`: `Σ z = a - b + |c| = a - (b - |c|)` and `Σ y_+ = Σ z_+ ≤ a`; equivalently
`Σ z_- ≤ b - |c|`, which forces `|c| ≤ b`; this is `T_{n-1}(a, b - |c|)`. For `c > a`, `Σ y_+ ≥ c > a`; for
`c < -b`, `Σ y_- ≥ |c| > b`. `T_1`: `y_1 = a - b` (and `(a-b)_+ ≤ a`). `T_2`: `y_2 = a - b - y_1`, and by the slice
description `y_1` ranges over `[-b, a]`. `□`

### 2.2 The peeling lemma and its counting form

**Lemma 2.2 (peeling; proved).** Let `S ⊂ Q^n` be finite, `φ` a linear form, and `S = S_1 ⊔ … ⊔ S_r` the partition
into nonempty level sets `S_i = S ∩ {φ = c_i}` (the `c_i` distinct), listed in **any** order. Then
`reg(S) ≤ max_{1≤i≤r} (reg(S_i) + i - 1)`.

*Proof.* Induction on `r`; `r = 1` is trivial. Let `d` be the right side, `g : S → Q` arbitrary, and
`S' = S_2 ⊔ … ⊔ S_r`. Choose `p` with `deg p ≤ reg(S_1) ≤ d` and `p = g` on `S_1`. On `S'`, `φ - c_1 ≠ 0`, so
`h := (g - p)/(φ - c_1)` is a function on `S'`. By induction (order `S_2, …, S_r`),
`reg(S') ≤ max_{i≥2}(reg(S_i) + i - 2) = d - 1` (if `r ≥ 2`; note `d ≥ reg(S_r) + r - 1 ≥ 1`), so there is `q` with
`deg q ≤ d - 1`, `q = h` on `S'`. Then `P := p + (φ - c_1) q` has degree `≤ d`, equals `g` on `S_1` (second term
vanishes) and on `S'` (`p + (g - p)`). `□`

(Here `reg(S_i)` may be computed inside the hyperplane `{φ = c_i}`: restriction from `Q^n` to an affine subspace maps
`Poly_{≤d}` onto `Poly_{≤d}` of the subspace.)

For a finite multiset `V` of integers put `peel(V) := max_i (v_(i) + i - 1)` with `v_(1) ≥ v_(2) ≥ …`, and
`N_V(t) := #{v ∈ V : v ≥ t}`. Since the first `i` sorted entries are `≥ v_(i)`, `i ≤ N_V(v_(i))`, hence

```text
peel(V) ≤ F   ⇐   N_V(t) ≤ F + 1 - t  for every t ∈ V.                                (2.2)
```

`peel` is monotone: raising entries of `V` does not lower `peel(V)` (compare sorted sequences entrywise).
Lemma 2.2 applied with the order of decreasing upper bounds `ρ_i ≥ reg(S_i)` gives
`reg(S) ≤ peel({ρ_i})`.

### 2.3 The bound functions and the recursion step

Define, for integers `x, y ≥ 0`, the *excess* functions

```text
H_2(x,y) = 0,
H_3(x,y) = min(x, y),
H_4(x,y) = min( ⌊(x+y)/2⌋ , 2 min(x,y) ),
H_5(x,y) = min( min(x,y) + ⌊max(x,y)/3⌋ , 3 min(x,y) ),
```

and `F_n(x,y) := x + y + H_n(x,y)`. Each `H_n` is symmetric, and each `F_n` is nondecreasing in each argument
(raising `x` by 1 raises `x + y` by 1 and does not lower any of `min(x,y)`, `max(x,y)`, `⌊(x+y)/2⌋`, hence does not
lower `H_n`).

**Lemma 2.3 (recursion step; proved).** Let `n ≥ 3` and suppose `F_{n-1}` is nondecreasing in each argument and, for
all integers `0 ≤ p ≤ a`, `0 ≤ q ≤ b`,

```text
(R_n)    min( H_{n-1}(a-p, b) + q ,  H_{n-1}(a, b-q) + p )  ≤  H_n(a,b).
```

Let `V = V_n(a,b)` be the multiset `{F_{n-1}(a', b) : 0 ≤ a' ≤ a} ∪ {F_{n-1}(a, b') : 0 ≤ b' ≤ b - 1}`. Then
`peel(V) ≤ F_n(a,b)`.

*Proof.* Note `F_{n-1}(a-p,b) - (a-p) - b = H_{n-1}(a-p,b)` etc., so `(R_n)` with `a* = a - p`, `b* = b - q` reads

```text
(R'_n)   min( F_{n-1}(a*, b), F_{n-1}(a, b*) ) - a* - b*  ≤  F_n(a,b) - a - b     (0 ≤ a* ≤ a, 0 ≤ b* ≤ b).
```

Fix `t ∈ V` and use (2.2). By monotonicity `{a' ≤ a : F_{n-1}(a',b) ≥ t}` is an interval `[α, a]` (or empty) and
`{b' ≤ b-1 : F_{n-1}(a,b') ≥ t}` is an interval `[β, b-1]` (or empty).
* Both nonempty: `N_V(t) = (a - α + 1) + (b - β)` and `t ≤ min(F_{n-1}(α,b), F_{n-1}(a,β))
  ≤ F_n(a,b) - a - b + α + β` by `(R'_n)`; so `N_V(t) ≤ F_n(a,b) + 1 - t`.
* Only the first: `N_V(t) = a - α + 1`, `t ≤ F_{n-1}(α, b)`. `(R'_n)` with `(a*,b*) = (α, b)`: the min is
  `F_{n-1}(α,b)` (monotone), so `t ≤ F_n(a,b) - a + α`, i.e. `N_V(t) ≤ F_n(a,b) + 1 - t`.
* Only the second: `N_V(t) = b - β`, `t ≤ F_{n-1}(a, β)`. `(R'_n)` with `(a*,b*) = (a, β)`: `t ≤ F_n(a,b) - b + β`,
  so `N_V(t) ≤ F_n(a,b) - t`.
(`t ∈ V` excludes the case that both are empty.) `□`

**Lemma 2.4 (`(R_3)`, `(R_4)`, `(R_5)` hold; proved).** For all integers `0 ≤ p ≤ a`, `0 ≤ q ≤ b`, writing
`X := H_{n-1}(a-p, b) + q`, `Y := H_{n-1}(a, b-q) + p`, we have `min(X, Y) ≤ H_n(a,b)` for `n = 3, 4, 5`.

*Proof.* The pair `(X, Y)` is exchanged by the symmetry `(a,p) ↔ (b,q)` (each `H_{n-1}` is symmetric), and `H_n` is
symmetric, so we may assume `a ≤ b`. `X, Y` are integers.
* `n = 3`: `X = q ≤ b`, `Y = p ≤ a`; `min(X,Y) ≤ min(a,b) = H_3(a,b)`.
* `n = 4`: `X = min(a-p, b) + q`, `Y = min(a, b-q) + p`. (i) `min(X,Y) ≤ (X+Y)/2 ≤ ((a-p) + q + (b-q) + p)/2 =
  (a+b)/2`, so `min(X,Y) ≤ ⌊(a+b)/2⌋`. (ii) `X ≤ b + q ≤ 2b`, `Y ≤ a + p ≤ 2a`, so `min(X,Y) ≤ 2 min(a,b)`.
* `n = 5` (`a ≤ b`, so `H_5(a,b) = min(a + ⌊b/3⌋, 3a)`):
  (i) `Y = H_4(a, b-q) + p ≤ 2 min(a, b-q) + p ≤ 2a + p ≤ 3a`.
  (ii) `X ≤ 2 min(a-p, b) + q ≤ 2(a-p) + q` and `Y ≤ ⌊(a + b - q)/2⌋ + p ≤ (a+b-q)/2 + p`. Hence
  ```text
  min(X,Y) ≤ (1/3) X + (2/3) Y ≤ (2a - 2p + q)/3 + (a + b - q)/3 + 2p/3 = a + b/3,
  ```
  and since `min(X,Y)` is an integer, `min(X,Y) ≤ a + ⌊b/3⌋`. `□`

### 2.4 The theorem

**Theorem 2.5 (proved).** For all integers `a, b ≥ 0` and `n ∈ {2,3,4,5}`: `reg T_n(a,b) ≤ F_n(a,b)`. Consequently:
1. `t_max(m) = reg(S_m) = ⌊10m/3⌋` for **every** `m ≥ 1`: no section of `m(-K_{LM_5})` (split or twisted) vanishes at
   `e` to order `> 10m/3`, and the Vandermonde products attain `⌊10m/3⌋`.
2. `τ_e(-K_{LM_5}) = τ_e(-K_{L̄_5^tw}) = 10/3`.
3. (`k = 4, 3`) `t_max(m) = 3m` for every `m` on `LM_4` and on `LM_3 = dP_6`; `τ_e = 3` in both cases.

*Proof.* Induction on `n`. `n = 2`: Lemma 2.1, `reg T_2(a,b) = a + b = F_2`. Step `n - 1 → n` (`n = 3,4,5`): by
Lemma 2.1 the slices of `T_n(a,b)` by `y_1` are `T_{n-1}(a', b)`, `0 ≤ a' ≤ a`, and `T_{n-1}(a, b')`,
`0 ≤ b' ≤ b - 1`, each once; by induction their `reg` is `≤` the corresponding entry of `V_n(a,b)`. Lemma 2.2 and
monotonicity of `peel` give `reg T_n(a,b) ≤ peel(V_n(a,b))`, and Lemmas 2.3-2.4 give `peel(V_n(a,b)) ≤ F_n(a,b)`.

1. By (2.1), `t_max(m) = reg T_5(m,m) ≤ F_5(m,m) = 2m + min(m + ⌊m/3⌋, 3m) = 3m + ⌊m/3⌋ = ⌊10m/3⌋`. The lower bound
   is `toric.md` Thm 3.4(a) (killer-recomputed): with `m = 3s + r`, `r ∈ {0,1,2}`, `V^s·V_{123}^r ∈ H^0(mB)`
   (`V` = 5-point Vandermonde `∈ H^0(3B)`, order 10; `V_{123}` = 3-point Vandermonde `∈ H^0(B)`, order 3) has order
   `10s + 3r = ⌊10m/3⌋`. Split and twisted forms have the same `H^0(mB)` over `Q(i)` (and `reg` is field-independent).
2. `τ_e = sup_m t_max(m)/m = sup_m ⌊10m/3⌋/m = 10/3`.
3. `LM_4`: `P_B = R_4`, `S_m = T_4(m,m)`, `F_4(m,m) = 2m + min(m, 2m) = 3m`; lower bound `V_4^{s}V_{123}^r` with
   `m = 2s + r` (`V_4 ∈ H^0(2B)`, order 6), order `6s + 3r = 3m`. `LM_3`: `F_3(m,m) = 3m`, attained by `V_{123}^m`. `□`

**Remarks.**
* The proof is the refereed Theorem B mechanism (slicing + peeling, `construct.md` §2.4), but run on the exact slice
  family `aΔ - bΔ` instead of through a mean-width hypothesis: slicing `mR_5` with mean-width bounds on the slices
  loses (it returns exactly Theorem B's `15m/4`), because the slices `T_4(a,b)` with `a ≠ b` are far from tight for
  mean width (e.g. `reg T_4(0,b) = b` against the C_bal(4) value `3b/2`). The non-width bound `H_4 ≤ 2 min` is what the
  weight `1/3` in Lemma 2.4 exploits.
* The functions `F_n` are not merely upper bounds for the recursion: the exact recursion
  `ρ_2 = a+b`, `ρ_n(a,b) = peel(V_n(a,b))` (with `ρ_{n-1}` in place of `F_{n-1}`) satisfies `ρ_n = F_n` for
  `n = 3,4,5` on `[0,90]^2` (exact; `rho_verify.py`). This is not used.
* `F_5` is **not** always equal to `reg T_5`: `reg T_5(2,3) ≤ 7 < F_5(2,3) = 8` (full rank mod `p` at `d = 7`,
  `treg5.out`). On the diagonal it is sharp by part 1.

## 3. The inequality (W): what is proved, what is open

### 3.1 (W) on the family aΔ - bΔ (proved)

**Corollary 3.1 (proved).** For all `a, b ≥ 0`, (W) holds for `S = T_5(a,b)` (and every translate), with room:
`reg T_5(a,b) ≤ F_5(a,b) ≤ (5/3)(a+b) = (1/6) Σ_{|A|=2} width_A(T_5(a,b))`, and `F_5(a,b) = (5/3)(a+b)` iff
`a = b` and `3 | a` (then `reg = F_5` by Theorem 2.5, so (W) is an equality on `T_5(3s,3s) = S_{3s}`).
*Proof.* On `T_5(a,b)`, `x_A` (`|A| = 2`) takes the values `a` (at `a e_i - b e_j`, `i ∈ A`, `j ∉ A`) and `-b`, and
`-b ≤ x_A ≤ a` (since `x_A ≤ Σ y_+ ≤ a` and `x_A ≥ -Σ y_- ≥ -b`); so each of the 10 widths is `a + b`. For `a ≤ b`,
`F_5(a,b) ≤ a + b + a + b/3 = (5/3)(a+b) - (b-a)/3 ≤ (5/3)(a+b)`, with equality iff `a = b` and `⌊b/3⌋ = b/3`;
symmetric for `a ≥ b`. `□`

So in particular every section `f ∈ H^0(aH + bH')` (Newton polytope in `aΔ - bΔ`, `toric.md` Prop 1.1(b)) satisfies
`ord_e f ≤ F_5(a,b)`; on the diagonal `a = b = m` this is Theorem 2.5. (Off the diagonal `F_5` need not be sharp:
`reg T_5(2,3) ≤ 7 < 8 = F_5(2,3)`, exact upper bound; `treg5.out`.)

### 3.2 (W) for arbitrary finite S: open

**Reduction (proved).** (W) holds for all finite `S ⊂ M` iff it holds for all *pair-box* sets
`P(ℓ,u) := {x ∈ M : ℓ_A ≤ x_A ≤ u_A, |A| = 2}` with `ℓ_A = min_S x_A`, `u_A = max_S x_A` for some finite `S`
(i.e. with tight bounds). *Proof.* `S ⊆ P(ℓ,u)` with the same 10 widths, and `reg` is monotone under inclusion. `□`

**Status.** Not proved. What the record has:
* *Evidence (round 22 J4):* 470 sets (random, lattice points of random polytopes, random zonotopes, half-space pieces of
  `mR_5`), 0 violations; equality on all root zonotopes tested.
* *Evidence (this round, `pairbox_peel.py` → `pairbox1.out`):* 60 random pair-box polytopes (`|S|` up to 1843). The
  iterated coordinate-peeling bound (slice by `x_i`, then by `x_j`, best ordered pair `(i,j)`, exact-mod-`p` `reg` of the
  2-dimensional slices) never exceeded the (W) right side (maximum ratio `1.000`, attained twice, both at small sets with
  `reg` equal to the (W) value); it equalled the mod-`p` `reg` in 31 of the 40 cases where `reg` was computed.
* *Certified instances (this round, `gp_test.py` → `gp_test.out`):* (W) holds (full rank mod `p` at the (W) value, hence
  over `Q`) on the lattice points of 150 random signed Minkowski sums of coordinate simplices `±Δ_I`; equality in 5.
  These are finitely many certified cases, not a proof. So the
  peeling mechanism of §2 is not visibly too weak for general (W); what is missing is a closed-form (or certified)
  bound for 4-dimensional pair-box slices, which have up to 14 independent facet parameters instead of 2.
* *Why the Theorem B template does not suffice (derived; a statement about one proof template, not about (W)).* In the
  slicing proof of Theorem B one needs, for two slices `c < c'` of `S` by a primitive `u`, a lower bound
  `Σ_A μ_A width_A(S) ≥ (c' - c) + (g_c + g_{c'})/2` with `g_c` a *mean-width* bound for the slice. Writing each
  functional as `g + α u`, only functionals in the same class modulo `u` can produce the term `c' - c` (the centres of
  the slice ranges of different classes are independent). For `μ = 1/6` on the 10 pair functionals and `u = x_1`, the
  classes with two members are `{x_{ab}, -x_{cd} = x_{ab} + u}` for the 3 splittings `{a,b,c,d} = {2,3,4,5}`, giving
  coefficient `3·(1/6) = 1/2 < 1`; for `u = x_1 - x_2` likewise `1/2`, for `u = x_1 + x_2` it is `1/6`. Using
  `Σ_b x_b = -c` to recover the missing `1/2` costs exactly the singleton widths of the slices, and the remaining slice
  bound is `(1/3)Σ_{cuts} L` against C_bal(4)'s `(1/2)Σ_{cuts} L`. Concretely, slicing `S_m` by `x_1` with the C_bal(4)
  bound `(3/2)(2m - |c|)` on the slice `c` and peeling returns exactly `15m/4` (Theorem B's value), while peeling with
  the *true* slice values returns `⌊10m/3⌋` (exact for `m ≤ 5`, `slice_peel.py`). Theorem 2.5 succeeds because it uses
  the true slice values, through the non-width bound `H_4 ≤ 2 min`.

**Open point W.** Prove (W) for all pair-box polytopes by the same peeling with a
certified piecewise-linear bound for the 4-dimensional slices; or find a counterexample. This no longer affects
`τ_e`: Theorem 2.5 already gives `τ_e = 10/3`. What general (W) would add is a *single* profile (uniform weight on the
ten 2|3 cuts) witnessing the cap for all auxiliary polynomials at once (§5, item 2).

## 4. Exact computations (none is used in the proof of Theorem 2.5)

All in `width_checks/`, run under `prlimit --as=3500000000`. `reglib.py`: `reg` by increasing `d` until the
evaluation matrix of all monomials of degree `≤ d` has full rank mod `p = 2^61 - 1` (full rank mod `p` ⇒ full rank over
`Q`, so the returned value is a certified **upper** bound for `reg`; minimality is checked mod `p` only, i.e.
evidence).

| script → output | what | result |
|---|---|---|
| `Rn_check.py` → `Rn_check.out` | `(R_n)` (Lemma 2.4) and monotonicity of `F_n`, exhaustively on `0 ≤ p ≤ a ≤ 40`, `0 ≤ q ≤ b ≤ 40` | all hold; moreover `H_n(a,b) = max_{p,q} min(X,Y)` exactly (`n = 3,4,5`) |
| `rho_dp.py`, `rho_verify.py` → `rho_verify.out` | exact recursive peeling numbers `ρ_n(a,b)` (`ρ_2 = a+b`, `ρ_n = peel` of the slice values) | `ρ_n = F_n` on `[0,90]^2` for `n = 3,4,5`; `ρ_5(m,m) = ⌊10m/3⌋`, `m ≤ 90` |
| `slice_peel.py` → `slice_peel.out` | `S_m` sliced by `x_1`, true slice `reg` (mod `p`), peeling | slices `T_4(m-|c|, m)` have `reg` `3m, 3m-2, 3m-3, …`; peeling gives `3, 6, 10, 13, 16` = `⌊10m/3⌋`, `m ≤ 5` |
| `treg_check.py` → `treg4.out`, `treg5.out` | `reg T_n(a,b)` against `F_n(a,b)` directly | `T_4(a,b)`, `0 ≤ a ≤ b ≤ 7`: full rank at `F_4` (certified `reg ≤ F_4`), deficient at `F_4 - 1` (evidence `reg = F_4`); `T_5(a,b)`, `0 ≤ a ≤ b ≤ 4`: full rank at `F_5` in all 14 cases, deficient at `F_5 - 1` except `(2,3)` (full at 7) |
| `gp_test.py` → `gp_test.out` | (W) on lattice points of 150 random signed Minkowski sums `Σ ±Δ_I` (generalized-permutohedron type, `|S| ≤ 268`) | `reg_p ≤` (W) in all 150 (certified, since `reg_p ≥ reg_Q`); equality in 5 |
| `pairbox_peel.py` → `pairbox1.out` | 60 random pair-box polytopes: `reg` (mod `p`), (W) right side, iterated coordinate peeling | iterated peeling `≤` (W) in all 60 (max ratio 1.000); `reg =` iterated peeling in 31/40 computed |

Consistency with the record: `t_max(m) = ⌊10m/3⌋` for `m ≤ 8` (round 22, author mod `2^61-1`, killer mod `2^31-1`) and
`t_max(m) = 3m` for `k = 4`, `m ≤ 9`, are now special cases of Theorem 2.5.

*Side data (not claimed).* The same recursion one level up gives `ρ_6(m,m) = ⌊10m/3⌋` for `m ≤ 30` (but `ρ_6 ≠ F_5`
off the diagonal, e.g. `ρ_6(1,10) - 11 = 4 > 3`), and `ρ_7(m,m) > 10m/3` (e.g. `ρ_7(5,5) = 17`). A proof for `LM_6`
would need a closed form for `ρ_6`; I did not attempt it. Since `V_5` and `V_{abc}` pull back along a coordinate
inclusion, `τ_e(-K_{LM_6}) ≥ 10/3`, so the side data suggest `τ_e(-K_{LM_6}) = 10/3` (**evidence**).

## 5. Consequences for A_5 and the barrier map

1. **(proved) `τ_e(-K) = 10/3` closes toric open point 1 of round 22.** Combined with `toric.md` Thm 3.3 (proved,
   killer-checked: `inf_F δ_F = 2/τ_e(B) - 1/2`, using `τ = sup_m t_max(m)/m`, `ord(ff^*) = 2 ord f`,
   `Q_f - Q_f ⊂ wR_5`):

   ```text
   inf over finite families F of δ_F  =  2/(10/3) - 1/2  =  1/10.
   ```

   *Every finite family of auxiliary Laurent polynomials over `Q(i)` (any degree, any width, any number), used point by
   point through Liouville with toric forced divisibility off their zero sets, yields for points of arbitrary profile
   at most `A_5(δ)` with `δ > 1/10`, i.e. Cilleruelo–Córdoba; it never improves it off a closed set.* This holds for the
   twisted `L̄_5^tw` and the split `LM_5` alike (same fan, same polytope, same `H^0(mB)` over `Q(i)`). The round-22
   statement that kill_toric F1 marked as an overclaim ("toric Liouville does not even improve Cilleruelo–Córdoba")
   is now **proved**, in the all-profile, all-finite-family form. The Vandermonde products attain it (`δ_F = 1/10` for
   `F = {V}`: `max_A width_V(1_A) = 6`, `ord = 10`, `6/10 - 1/2 = 1/10`).
2. **Quantifier of the profile (constraint check).** Thm 3.3's minimax gives, *for each finite `F`*, a cut-weight
   profile `ℓ_F` at which no element of the semigroup generated by `F` beats `δ = 1/10`; `ℓ_F` may depend on `F`. The
   threshold `1/10` is uniform; the witnessing profile is not claimed uniform. If the general (W) were proved, the single
   profile `ℓ*` = uniform weight `1/10` on the ten 2|3 cuts (zero on the 1|4 cuts) would witness the cap for all `F`
   at once: by (W), `Σ_{|A|=2} width_f(1_A)·(1/10) / ord_e f ≥ 6/10` for every `f`.
3. **(proved) McKinnon–Roth base-locus Liouville at `e`** (`lit.md` Prop 1.3 / row 20) gives exactly
   `α ≥ τ_e(-K) = 10/3`, which in the dictionary `γ = 4/(1+2δ)` is `δ = 1/10`: the literature's section/base-locus
   constant at `e` *is* Cilleruelo–Córdoba. Seshadri-proportional bounds stay at `ε_e = 2`.
4. **No effect on the target.** `A_5(δ)` for `δ < 1/30`, i.e. `α^ess_e(-K) > 15/4`, and `(L_5^∘)` are untouched:
   Theorem 2.5 *lowers* a section cap (from `[10/3, 15/4]` to exactly `10/3`); it makes the window
   `α^ess ∈ (15/4, 4)` even further from every section argument. On **balanced** profiles the section cap is still
   `15/4` (`δ = 1/30`, `construct.md` Cor C3), attained by `V`; nothing here changes that.
5. **Barrier map, row 20 (proposed edit).** Replace "best exponent `δ = 2/τ_e(-K) - 1/2`, with `τ_e ∈ [10/3, 15/4]`" by
   "best exponent `δ = 1/10` for arbitrary profiles (`τ_e(-K) = 10/3`, round 23 Thm 2.5; minimax), `δ = 1/30` for
   balanced profiles"; status proved. Expert note §3, "Approximation reading": replace
   "`10/3 ≤ τ_e(-K) ≤ 15/4` (proved). Exact data … suggest `τ_e = 10/3`" by "`τ_e(-K) = 10/3` (proved, round 23):
   `t_max(m) = ⌊10m/3⌋` for all `m`". In round 22 README §2 the row "toric sections" becomes `δ > 1/10` (arbitrary
   profiles). The "intermediate-exponent lane" (`δ ∈ (1/30, 1/10]` for arbitrary profiles) is closed for toric
   Liouville.
6. **What this does *not* say.** It is a barrier for one method class (sections with toric forced divisibility), not a
   statement about the existence of A-points; it says nothing about non-section mechanisms, about `α^ess`, or about
   `(CB_e^norm1)`.

## 6. Constraints audit, labels, files, open points

**Owner constraints, item by item.**
* *Uniform exceptional sets from fixed-`S` Subspace Theorems:* none used. The only exceptional sets are the zero sets
  `Z(f)` of the auxiliary polynomials, inherited from `toric.md` Thm 3.2/3.3.
* *Gaussian gcd or denominator costs:* Theorem 2.5 is pure interpolation. Its consequence uses `toric.md` Thm 3.2, which
  keeps `N'^m` (first form) and the exact toric heights `h_{E_A} = ½ log n_A` of primitive tuples. Nothing is discarded.
* *Drifting common phase:* not involved. `T_5` quotients the phase exactly, and `τ_e` is a property of `(X, -K, e)`.
* *Template-dependent threshold treated as uniform:* the threshold `1/10` is proved as an infimum over all finite
  families. The witnessing profile is explicitly **not** claimed uniform in `F` (§5 item 2). No statement depends on the
  radius.

**Labels.**
* Proved: Lemmas 2.1-2.4, Theorem 2.5, Corollary 3.1, the reduction in §3.2, §5 items 1 and 3 (with `toric.md`
  Thm 3.3, killer-checked).
* Exact/certified: the instances in §4 where "full rank mod `p`" is reported.
* Evidence: the deficiency-mod-`p` statements, the general-(W) tests, and `ρ_6`, `ρ_7`.
* Derived (about a proof template, not a theorem about (W)): the coefficient-`1/2` remark in §3.2.
* Nothing is refereed.

**Novelty check (limited).** I grepped the record for the slicing proof. Theorem B (`construct.md` §2.4) is the same
peeling mechanism, applied with mean-width hypotheses. The exact-slice recursion for `aΔ - bΔ` and the closed forms
`H_3, H_4, H_5` do not appear in round 22 or `construct.md`. I did not search the literature on interpolation
regularity of lattice points of dilated root polytopes. A statement equivalent to `reg((mΔ - mΔ) ∩ Z^5) = ⌊10m/3⌋` may
well be known there, so treat it as "new to the record".

**Files** (`scratchpad/round23/`): `width.md` (this file), `width_notes.md`, and `width_checks/`:
* `reglib.py`
* `rho_dp.py`, `rho_explore.py`, `rho_verify.py` → `rho_verify.out`
* `Rn_check.py` → `Rn_check.out`
* `slice_peel.py` → `slice_peel.out`
* `treg_check.py` → `treg4.out`, `treg5.out`
* `pairbox_peel.py` → `pairbox1.out`
* `gp_test.py` → `gp_test.out`

The repository was not edited.

**Open points.**
1. General (W), i.e. C_bal(5), for all finite `S` (§3.2). Map value only. It would supply one uniform witnessing
   profile for the cap of §5.
2. `τ_e(-K_{LM_6}) = 10/3`? Evidence from `ρ_6(m,m) = ⌊10m/3⌋`, `m ≤ 30`. It needs a closed form for `ρ_6`. For `k ≥ 7`
   the peeling bound exceeds `10m/3`.
3. Off-diagonal: the exact `reg T_5(a,b)` (equivalently the vanishing thresholds of the classes `aH + bH'`) is not
   determined. `F_5` is an upper bound, and it is not exact at `(2,3)`. This is not needed for `τ_e(-K)`.
4. Unchanged: the target `α^ess_e(-K) > 15/4` (Question B / `(CB_e^norm1)`). This round neither moves nor obstructs it
   beyond the section class.
