# Referee report on `round6/sharp.md`: arithmetic and class definitions

Lens: (1) Def. 1.1 classes `F^char_X`, `F^all_X`, Prop. 1.2 (actual clusters are members) and
Prop. 1.3 (`two_thirds` inside); (2) Prop. 2.1 and Cor. 2.2 (norm-one group, forced parity, units,
twins, explicit preimage (2.2), levels); (3) Lemma 3.1, Def. 3.2 and Lemma 3.3 (dyadic assignment,
pair demand); (4) the negative claims of §7. I also re-checked the steps of `round5/lower.md` that
`sharp.md` uses: Prop. 1.2/1.3, Lemma S, Prop. 3.3 (= sharp Prop. 4.3), Lemma 3.1 (= sharp
Lemma 4.1), Lemma 5.2 (= sharp Lemma 5.1b) and Lemma 5.3 (= sharp Lemma 5.1a).

Method. Every step was re-derived by hand. Every testable claim was then recomputed with my own
code in `round6/referee_arith_checks/a_*.py` (exact integer or Gaussian-integer arithmetic unless
stated). The provided scripts were rerun, and their outputs are in `referee_arith_checks/rerun2/`.

* `check_parity.py`: output identical to the stored one.
* `check_realise.py`: output identical, 64 s.
* `check_design.py`: identical conclusions. The stored output was produced with `KMAX = 22`; the
  default `KMAX = 24` prints the range `2^25` with the same numbers.

The folder also contains `ref_*.py` and `rerun/`. These come from an earlier, unfinished referee
attempt, and I did not rely on them. Where they overlap with my checks they agree, and they
independently flag the Prop. 6.5 per-unit point and the `e^(5.3M)` count.

---------------------------------------------------------------------------------------------

## 0. Verdict

**Within this lens the results survive.** I found no fatal or major error. Everything below is
minor: quantifier slips, a lossy bound presented as if monotone, an overstated sentence and numeric
slips in heuristic paragraphs. None of them touches Theorem 5.1, Corollary 5.2 or Theorem 6.4.

| item | status |
|---|---|
| Def. 1.1; Prop. 1.2 (actual clusters lie in `F^all_X` for every `X`) | **correct**. Re-derived. Exact check on 15 genuine clusters (exponents up to 3, nonzero windings), all characters `||c||_1 <= 6`, all monomials `||v||_1 <= 4`, all units, same-level clause |
| Prop. 1.3 (`two_thirds` inside `F^char_X`, `X >= 2M`) | **correct**. `two_thirds` uses the configuration only through (2.5) and (3.6). The collisions it needs sit at `q^a < 1.25 M` |
| `F^char_(3M)` is a subclass of III-small (lead notes §1, §6) | **correct**. III-small's moduli satisfy `q^a < 2.5M`, and it strengthens with unit 1 only |
| Prop. 2.1 (1)-(6) | **correct**. Exhaustive over `Z[i]/Q` for all 184 odd prime powers `Q <= 1000`, items (4) for **every** `t` and **every** `n` |
| Cor. 2.2 (1)-(5) and (2.2) | **correct**. Includes general exponents, random windings, all levels, every `c` in `{-2..2}^M`, and exhaustive realisable sets |
| Remark 2.3(2) (quartic part) | **correct**; see the scope remark in §3.4 |
| Lemma 3.1 | **correct**. Exact multiplicity `<= 4` for all odd `q < 2^27`. `R(k) < 8.109` for all `k >= 20`. Item 3 holds for all `d <= 2*10^6` |
| Def. 3.2, Lemma 3.3(1),(2),(4) | **correct**. Independent class-level implementation (Paley `M = 108`, `X = 3M` and `3M^2`). The pair pattern is exactly that of Lemma 3.3(1) |
| max pair demand `<= (10 + o(1)) log M` | **true but lossy**. Recomputed: `5.7`-`7.2` times `log M` for `20 <= M <= 262144`. The chain from Lemma 3.1 to `mu_M` needs a `max omega` step (issue 3) |
| Prop. 7.1 (Siegel trap) | **correct**, given `v(c) != 0` (issue 8) |
| §7.2 / §7.4 (coherent designs) | **overstated** (issue 5); Prop. 7.4 needs its pair hypothesis for all units (issue 4) |
| §7.3 heuristic | numeric slips (issue 6); the qualitative conclusion stands |

---------------------------------------------------------------------------------------------

## 1. The classes, Prop. 1.2 and Prop. 1.3

### 1.1 Def. 1.1 and Prop. 1.2: re-derivation

**Units and windings.** Write `z_x = i^(u_x) prod pi_j^(a_xj) conj(pi_j)^(e_j - a_xj)`, and take
`theta_0 = theta*` (the arc start). Then (1.1) of `outside.md` gives `k_x = -u_x mod 4`. Any other
admissible `theta_0` differs by a multiple of `pi/2` and shifts all `k_x` together. So
`z_x = eps_0 i^(-k_x) prod(...)` with one common unit `eps_0`, and `rho_x(q) = conj(eps_0) z_x`.
A common unit cancels in `prod rho^c` because `sum c = 0`. Shifting a lift `phi_j` by `2 pi`
changes `k_x` by a multiple of 4, so `i^(-k_x)` is well defined.

**(e-char).**
* `U/V = e^(i c.delta)`.
* `G = gcd(U, V)` with `U/G = u_1 A_v` and `V/G = u_2 conj(A_v)`.
* `|G|^2 e^(w(v)) = N^(n/2)`.
* `D_eta = U/G - eta V/G` is nonzero (conjugate-primitivity) and divisible by `(1+i)`.
* `q^(a_(q,eta)) | D_eta` because `q` is prime to `G`.

Together these give `2 N^(n/4) |sin((c.delta - arg eta)/2)| = |U - eta V| = |G||D_eta| >= N^(n/4) e^(-w/2) sqrt2 m`.

**(e-all).** I checked the four identities `A - i^s conj A = 2iY, (X-Y)(1-i), 2X, (X+Y)(1+i)`.
Also `|D_s| = nu_s^(-1)|A_v||sin(<v,phi> - s pi/4)|`. A Gaussian factor `1 +- i` is prime to odd
`q`, so `q^a` divides the rational integer `D_s`.

**Same-level clause.** The factored level residue at `p_l` is `conj(eps_0) r_x` with `two_thirds`'
`r_x`, and `two_thirds` Lemma 3.3 (S|) applies.

*Exact check* (`a_cluster.py`, output `a_cluster_out.txt`).
* Data: 15 genuine normalised clusters (`N` up to `3.1*10^10`, exponents 1-3, `M = 5..7`), with
  windings up to `|k| = 5` occurring naturally.
* Integrality of `k` and the common unit `eps_0`: checked.
* `rho_x(q) = conj(eps_0) z_x mod q^(A_q)` for every odd `q <= 60`: checked.
* About 22,000 characters: checked exactly.
  * every unit: `D != 0`, `(1+i) | D`, `q^(a_(q,eta)) | D`, `Norm(G) prod p^|v_j| = N^(n/2)`, and
    `Norm(U - eta V) = Norm(G) Norm(D)`;
  * the float value of the sine agrees with the exact chord.
* About 16,500 monomials with all four `s`: checked.
* About 450 same-level pairs: checked.
* Cor. 2.2(4) with the actual windings: checked.

ALL PASSED.

### 1.2 Prop. 1.3

`two_thirds` Theorem 7.5 uses the configuration through two inequalities.
* (2.5) is the pair case `eta = 1` of (1.1) summed against the cut identity. I re-derived
  `sum log(sqrt2 m_xy) + Q/2 <= (M/8) W + K log C`.
* (3.6) is the pigeonhole. Its class counts depend only on the norms: `N mod q^a` for `q` not
  dividing `N`, and `N/p^e mod p^a` for same-level residues at `p | N`.

The collisions used satisfy `(q - chi(q)) q^(a-1) < M`, i.e. `q^a < 1.25 M <= X`. sharp.md says
`< 1.5M`, which is also true. The level data needed are at `p <= M + 1 <= X`. The chain lemma, the
prime sums and the normalisation `min a = 0`, `max a = e` are in (a). So Prop. 1.3 holds for every
member of `F^char_X`, `X >= 2M`. The ramified `sqrt 2` is the `2^(-1/2)` in (1.1).

### 1.3 Relations

* `F^all_X ⊂ F^char_X ⊂ A_X`: correct. The point and level residues have norms `N` and
  `N/p_l^(e_l)`, and (e) of `A_X` is (1.1) with `eta = 1`.
* `F^char_(3M) ⊂ III-small` (lead notes §1, §6): correct. III-small strengthens only with unit 1,
  at moduli with `m_(q^a) < 2M`, i.e. `q^a < 2.5M <= 3M`. Its labels are free, and `F^char` adds
  constraints.

---------------------------------------------------------------------------------------------

## 2. Prop. 2.1 and Cor. 2.2

### 2.1 Re-derivation

Prop. 2.1:
* (1) For split `q`, `T ≅ (Z/q^a)^*`. For inert `q`, `T ≅ mu_(q+1) x (trace-zero part of qO/q^a O) ≅ C_(q+1) x C_(q^(a-1))`.
  The orders are coprime, so `T` is cyclic. `4 | q - chi(q)`.
* (2) The kernel of `r -> r/conj r` is `(Z/Q)^*` (`q` odd), and fibre counting via `two_thirds`
  Lemma 3.2 gives surjectivity.
* (3) The index-2 subgroup of a cyclic group of even order is the squares.
* (4) Hensel.
* (5) `i = (1+i)/conj(1+i)` and `(2/q)`.
* (6) Surjective homomorphisms of cyclic groups map generators to generators.

All correct.

Cor. 2.2:
* `nu_j(q)` is forced. The map from `r_j` to `ell_j` is 2-to-1 (`r_j` and `-r_j`), and that sign
  is a common factor.
* `rho_x = B g^(kappa_x)` with `kappa = A ell - iota k`.
* `iota` is odd iff `q = +-3 mod 8`. This gives the `lambda_q k_x` term.
* (2.2) works because `A(e_f(x) - e_s(x)) = eps_x e_x`. The realised classes are `kappa - kappa_0`,
  so classes are defined up to a common shift, which is harmless.
* *Side remark.* Distinctness of the `2M` twin columns is **not** needed for (2.2), by linearity.
  It is needed for Lemma 6.2.

### 2.2 Exact checks

**`a_groups.py`** (output `a_groups_out.txt`; numpy, does not import `gres.py`). All 184 odd prime
powers `Q <= 1000`, 17 of them with `a >= 2`, were checked exhaustively over `Z[i]/Q`:
* `|T| = m` and `4 | m`;
* cyclicity: an element of order `m`, whose powers exhaust `T`;
* `r -> r/conj r` is onto `T`, with fibres of size `phi(Q)` equal to `r_0 (Z/Q)^*` for **every** `t`;
* the Legendre symbol of `Norm r` is constant on fibres and equals the parity of `log_g t`;
* item (4) for **every** `(t, n)`: there are exactly 2 such `r` iff `nu(t) = (n/q)`, and 0
  otherwise;
* `i` is a square iff `q = +-1 mod 8`, and `log -1 = m/2`;
* item (6): reduction is onto, generators go to generators, `nu` is compatible with independently
  chosen generators, and `log` reduces mod `m_(q^(a-1))`.

ALL PASSED. sharp.md's "89 prime powers `<= 400`" is confirmed. Its own script tests only 40 values
of `n` per modulus; this check covers all of them.

**`a_cor22.py`** (output `a_cor22_out.txt`). Two random profiles were tested.
* P1: `M = 6`, `e = 1`, 20 columns. P2: `M = 5`, exponents up to 2, 16 columns.
* Actual primes `p_j > 10^4`.
* All 30 odd primes `q <= 130` with their top powers (`3^4`, `5^3`, `7^2`, `11^2`).

Results:
* C1: the attainable `ell_j` are exactly `{ell = nu_j mod 2}`, each attained twice; `pi_j` has the
  forced parity; `nu_j(q) = [(q/p_j) = -1]`.
* C2: about 57,700 character/level tests. For every `c` and level `a`, `prod rho^c = i^s` iff
  `sum c kappa = iota s mod m_(q^a)`, with random windings `|k| <= 9`.
* C3: parity of `kappa` differences equals `sigma_q` differences including `lambda_q k`;
  `iota mod 2 = [q = +-3 mod 8]`.
* C4: the explicit preimage (2.2) with random windings and a random odd/even offset. The `r_j` are
  built by Hilbert 90 and a square root, and realise `kappa* - kappa_0` exactly.
* C5: exhaustive at `M = 3`, 7 columns, `q = 3, 5, 7, 13, 17`. The set
  `{A ell - iota k : ell = nu mod 2}` equals the parity-respecting set: sizes 8, 8, 64, 216, 512.

ALL PASSED.

### 2.3 What "free otherwise" means

`ell_j mod 4` and the higher digits are unconstrained by the norm alone: (4) holds for every `t`,
and the realisable set is the full parity class. Prime powers add only reduction compatibility. The
unit `-1` is always a square, and `i` is a square iff `q = +-1 mod 8`. All of this is confirmed.

### 2.4 Scope remark on "exactly one bit" (not an error)

The statement is exact for the class, whose data are residues with prescribed **norms**. Genuine
residues `pi_j mod q` satisfy more, through reciprocity.

* For inert `q = 7 mod 8`, `ell_j mod 4` (when `(p_j/q) = 1`) equals the quartic symbol of `-q`
  modulo `pi_j`. This is quartic reciprocity, `[pi/q]_4 = [-q/pi]_4`. It is determined by data at
  the modulus `p_j`, not at `q`.
* Checked in `a_cor22.py` C6 for `q = 7, 23, 31, 47, 71, 79, 103` against all `p <= 3000`. The
  same check confirms invariance under units and conjugation when `8 | m`, and that both values
  occur.

For `X < min p_j`, as in Theorem 5.1, this is invisible to local data at `q <= X`. So Cor. 5.2's
"every parity (quadratic reciprocity) consequence" is accurate. In the §7 regime `X >= p_j`,
quartic reciprocity among the varying primes would couple `ell_j(p_l)` and `ell_l(p_j)`. The class
does not encode this. A sentence would prevent "one bit" being read as "every reciprocity
consequence".

---------------------------------------------------------------------------------------------

## 3. Lemma 3.1, Def. 3.2, Lemma 3.3: the assignment and the demand

### 3.1 Lemma 3.1 (`a_assign.py`, output `a_assign_out.txt`)

* (1) holds for every odd prime `q < 2^27`: `2p' <= m_q`, `q < 8p'`, `p' <= (q-1)/2` for `q >= 5`,
  `p' != q`, and every `q` is assigned exactly once.
* (2) The exact maximum multiplicity is 4 in every block `k <= 26`. Prime 2 has 3 preimages (3, 5,
  11), prime 3 has 2 (7, 13). The maximum `n/t` is 3.70.
* R-S part: I recomputed `R(k)`. Its maximum over `20 <= k <= 5000` is 8.1071, and the limit is
  `1.51012/0.186235 = 8.1087 < 9`, matching sharp.md's 8.1047 for `k <= 2000` and 8.1306 for the
  tail. The two R-S inequalities were checked on the whole sieve range.
* (3) `Lambda'(d) <= 9 log d + 18.72 omega(d)` holds for all `2 <= d <= 2*10^6`. The true maximum
  of `Lambda'(d)/log d` is 7.37, at `d = 2`. The proof `sum_(p|d) 9 log(8p)` is correct, with
  `9 log 8 = 18.715`.

### 3.2 Def. 3.2 and Lemma 3.3

Re-derived:
* the CRT identification `Z/m_(q^a) ≅ Z/m_q x Z/q^(a-1)` is compatible with reduction;
* level-1 classes are injective modulo the parity because `2p' - 1 < m_q`;
* level `a*` is injective because `gcd(p', q) = 1` and `p' q^(a*-1) >= M`;
* random lifts keep the parity because `m_(q^(a-1))` is even;
* the probabilities in (3) and (4) hold: `m/2 - M + 1 > m/4` for `m >= 4M`, and
  `3/(m/2 - 1) <= 12/m`.

**`a_design.py`** (output `a_design_out.txt`) is an independent class-level implementation:
* canonical Paley `q_0 = 107`, `M = 108`, `b = 33`, 3531 actual primes `> 10^6`;
* every odd prime power `<= X`, for `X = 3M` (A = 1) and `X = 3M^2` (A = 2).

Results:
* D1: the top level has parity `sigma_q`, and every level is the reduction of the next.
* D2: unit-1 pair collisions occur **exactly** when Lemma 3.3(1) predicts them. There are none at
  large `q` or above `a*`.
* D3: `log m_(xy,1) <= Lambda'(d) + log d` for every pair. The maximum is `5.31 log M`.
* D4: at large `q`, pair collisions with units other than 1 number 9,407, below the Lemma 3.3(4)
  bound of 37,379.

`check_realise.py` was rerun (identical output). It covers the residue level end to end at
`M = 44`, `X = 2M^2`.

**Lemma 3.3(2).** The maximum of `q^(a*)/M^2` is 16 (claimed `< 32`). `log m_struct = (7.35 ... 8.00) M`
(claimed `<= 11 M`).

### 3.3 Max pair demand, recomputed

Here `D(d) = sum_(q < 4M, q <= X, p'(q) | d) log q * min(a*(q), A_q, 1 + v_q(d))`, in the worst case
over parity colourings (`a_assign.py` Part D, `a_demand_scan.py`).

| `M` | 44 | 64 | 256 | 1024 | 4096 | 16384 | 65536 | 262144 |
|---|---|---|---|---|---|---|---|---|
| `max D / log M` (A = 1 and A = 2) | 5.93 | 5.73 | 5.94 | 6.61 | 6.77 | 6.54 | 6.37 | 6.50 |

* The values for `M = 64 ... 16384` reproduce `check_design.py` exactly.
* Scanning every `M` in `[20, 2000)` and a grid up to 20,000 gives a maximum of 7.00 (`M = 111`,
  A = 1) and 7.16 (`M = 23`, A = 2).
* So sharp.md's "measured `<= 6.8 log M`" holds at the five sampled `M`, not uniformly. This is
  evidence only.
* The proved `(10 + o(1)) log M` is valid.

---------------------------------------------------------------------------------------------

## 4. The lower.md inputs used by sharp.md

* **Lemma S (sharp Lemma 4.2).** Re-derived line by line:
  * defect identity, window `|F - kappa n| <= 2(F - M)`, inversion;
  * Lemma 2.2: the coefficient of `Y^k` is `binom(d,k) sum c_x (-x)^(d-k)`, with Vandermonde on
    `<= d` points;
  * the small range: `f(4) = M/4 - 4`, `f(6) = (4/3)(M/3 - 6)`, concavity on `[8, 19M/48]`, and
    all three are `>= M/16` for `M >= 44`;
  * `kappa = 1` gives `|T_a| > (1 - 4 theta) M`, and `kappa = 2` gives `phi >= (1/2 - 6 theta) M`.

  Correct. `a_lowerS.py` (output `a_lowerS_out.txt`, `q_0 = 43, 47`):
  * Hadamard property;
  * exhaustive ternary support 4: `min F/M = 1.364, 1.333` against `17/16`;
  * 20,000 random larger supports;
  * amplitude `>= 2` gives `F >= 2M`;
  * pairs, columns and half-sums all have `F = M`.
* **Prop. 3.3 (sharp Prop. 4.3) and identity (3.1).** Verified on every tested `c` at `b = 33`:
  * pairs `>= b - 4`, attained;
  * columns `>= b + 2M - 8`; half-sums `>= b + M - 16`;
  * support 4: `>= 557` (needs `kappa_b M = 2.75`);
  * amplitude bound `hM(b-4)/2 + b`;
  * slack inequality `s_c >= (tau/2)(2B - r) - 1/2` with random weights.
* **Lemma 3.1 (sharp Lemma 4.1).** Twins `S_f(x) - S_s(x) = -2H(x, alpha x) e_x` with `2M`
  distinct twin columns; exact `rank S = M` (mod `10^9 + 7`); every column nonconstant. Saturation
  follows from the twins.
* **Lemma 5.3 (sharp Lemma 5.1a).** Re-derived:
  * the density of `c.delta/2` is at most `2/Delta`;
  * its range has length `n Delta/2`;
  * at most `2n Delta/pi + 2` bad intervals, each of length `<= pi eps`.

  Correct.
* **Lemma 5.2 (sharp Lemma 5.1b).** Correct in all its uses, but see issue 2.

---------------------------------------------------------------------------------------------

## 5. Section 7 (negative claims)

* **Prop. 7.1.**
  * Siegel's lemma (`m = 2` equations, `n = M` unknowns, height `max(1, max|lambda|)`) gives
    `||c||_inf <= (M max(1, max|lambda|))^(2/(M-2))`.
  * `sum c kappa = 0 mod m_Q` at every `Q`.
  * `sum_(Q in 𝒬) log q <= log m_(c,1) <= w/2 + (1/2) log 2`.
  * `w(v(c)) <= nW/2`, since `|v_j| <= (n/2) e_j`.

  All correct. `a_sect7.py` S1 finds a solution within the Siegel height in 100/100 random
  labellings (`M = 6, 7`). The one omission is issue 8.
* **§7.2.** As a statement about the **bound** `m^gen_c <= Y^(2n)/4` of `lower.md` Lemma 4.1, it is
  correct: `4 log Y > 4AW` exceeds the average pair slack `W/(2(M-1))`. That coherent designs need
  `|zeta_x|^2 > X` is also correct. But §7.2 does not show that coherent designs fail (issue 5).
* **§7.3.** Heuristic, and labelled as such. Its numbers are off (issue 6), but the qualitative
  conclusion survives and gets stronger.
* **Prop. 7.4.** The union-bound bookkeeping is correct, except that the pair hypothesis must cover
  units other than 1 (issue 4).
* **The Products-of-`k`-points remark.** Correct: pigeonhole over `binom(M+k-1, k)` multisets.

---------------------------------------------------------------------------------------------

## 6. Issues (all minor)

1. **Prop. 6.5 uses `m_(e_j)` (maximum over units) instead of the per-unit height condition.**
   * Prop. 6.1 bounds `m_(v,s) <= nu_s^(-1) e^(w/2)` for each `s` separately.
   * The product over `s` can be as large as `e^(2w)/4` even for genuine `A_v`, since `m_v` divides
     `Y(X-Y)X(X+Y) = Im(A^4)/4`.
   * At `q` in `Q_0` (`q = +-3 mod 8`) a hit of `e_j` has unit `s` of parity `nu_j(q)`, and both
     units of that parity are equally likely. So the hits are spread over at least two units, and
     `m_(e_j) > sqrt2 e^(w_j/2)` does not contradict Prop. 6.1 as written.
   * **Fix.** Split `Q_0` by `nu_j(q)`. One part has `sum 2/m_q >= lambda/4`. Fix the unit `s = 0`
     (respectively `s = 1`) on it. Then each `q` of that part gives a hit with that unit with
     probability `2/m_q`, independently. The Poisson tail for `kk` hits is still
     `exp(-(1 - eta_0 + o(1)) log M)`, because `lambda/4` is a constant. So the conclusion of
     Prop. 6.5 holds verbatim for `m_(e_j, s)`.
   * Nit: the number of non-twin columns is `(b-2)M - b`, not "at least `(b-2)M`".
2. **Lemma 5.1b (and `lower.md` Lemma 5.2): "its collision modulus can only drop" is false in
   general.**
   * If a deleted column's prime satisfies `p_l <= X'`, then `q = p_l` no longer divides `N'`.
   * (e-char) and (e-all) then impose collisions at `p_l` on **all** characters and monomials.
     Before, only same-level pairs were constrained there.
   * Harmless in every use, since `p_j > X`. **Fix:** add the hypothesis "every deleted `p_l > X'`",
     or "all `p_j > X`".
3. **The demand chain in Lemma 3.3(1) is not monotone as written.**
   * The per-`d` form `Lambda'(d) <= (9 + 26/loglog d) log d` stated after Lemma 3.1 is true for
     `d >= 3` but large at small `d`: about 315 at `d = 3`.
   * Bounding `Lambda'(d) + log d` by its maximum over `d < M` through this form exceeds
     `mu_M log M` for all `M < 1.19*10^7`. At `M = 2^20` the per-`d` route gives `22.6 log M`,
     while `mu_M = 19.9`.
   * **Correct route:**
     `Lambda'(d) + log d <= 10 log d + 18.72 omega(d) <= 10 log M + 18.72 max_(d<M) omega(d) <= mu_M log M`.
     This holds for all `M >= 30`: the primorial maximiser together with Robin's bound and the
     monotonicity of `log x/loglog x` on `[e^e, inf)`. The `max omega` inequality was checked
     numerically for every `3 <= M <= 2*10^6`.
   * The true values are `max_(d<M)(Lambda'(d) + log d) <= 6.94 log M` at `M = 2^4, ..., 2^20`.
   * Nit in Def. 3.2: the levels `2 <= a <= a*(q)` should read `2 <= a <= min(a*(q), A_q)`.
4. **Prop. 7.4: the pair hypothesis must cover units other than 1.**
   * The conclusion is about `F^char`, where (1.1) holds for every unit.
   * With moduli up to `N^A`, pair collisions with `eta = -1` or `+-i` can exceed `e^(W/4)` while
     the unit-1 moduli stay small.
   * Example (`a_sect7.py` S4): in a coherent design `zeta_x = n_x + i`, unit-1 pair moduli are at
     most `e^(4.7)` but unit `-1` pair moduli reach `e^(23)`.
   * **Fix:** require `log m^gen_xy <= C_1 log M`, or `log m^off_xy <= W/4 - 1`.
5. **"Coherent ones fail by §7.2" (§7.4, and the corresponding bullet in §0 item 6) is
   overstated.**
   * §7.2 shows only that the a-priori bound `2n log Y` certifies nothing for `X >= N^(1/4)`.
   * It says nothing about the actual collision moduli of a coherent design. These can be far
     below the bound.
   * In S4, `zeta_x = n_x + i` with `n_x^2 + 1` prime `> X` has unit-1 pair moduli dividing
     `n_y - n_x`, independent of `X` (checked directly at all `q^a <= X`). The bound gives 92.
   * Coherent designs plausibly still fail through units other than 1 and through support-4
     characters, whose `X`-smooth parts are typically `X^(Theta(1))`. In S4, 59% of support-4
     characters have `log m > (log X)/2`. But this is heuristic.
   * **Corrected sentence:** "the Lemma 4.1 bound certifies nothing for `X >= N^(1/4)`; we know no
     coherent design that works, and heuristically generic characters collide at `X`-smooth moduli
     of size `X^(Theta(1))`."
6. **§7.3 numbers.**
   * `#{c : ||c||_1 <= 4M}` is `e^((3.085 + o(1))M)`, not "about `e^(5.3M)`". The value 5.3 is the
     crude bound `2^n binom(M+n-1, n)` at `n = 4M`. Exactly:
     `sum_k 2^k binom(M,k) binom(4M,k)`, with exponent 3.085 (`a_sect7.py` S2).
   * The per-collision factor `(log 2/y)^kk` should be `(log 2)^kk`. Mertens gives
     `sum_(e^y < q <= e^(2y)) 1/q -> log 2`; S3 measures 0.685-0.690.
   * Corrected: `kk ~ 3.1M/log M` and demand about `3M^2/log M` at `y ~ M`.
   * Since the slack of such characters is about `M^(3/2) log M`, the heuristic failure point is
     already about `y ~ M^(1/2) log^2 M`, i.e. moduli `e^(M^(1/2+o(1)))`, which is inside the
     claimed `e^(Theta(M))`. Qualitatively unchanged.
7. **§6.6.** With `X = 3M`, `pi(X) log 2 ~ 3M log 2/log M`, not `M log 2/log M`.
8. **Prop. 7.1 needs `v(c) != 0`** for (1.1) to apply. This is automatic on the §4 profile (twins,
   rank `M`). For general profiles, add it as a hypothesis.
9. **Lemma 6.2's parenthetical "or any profile with `2M` distinct twin columns" needs `e = 1`**, or
   bounded exponents with a modified constant. The last step uses
   `|(c'^T a)_j| <= ||c'||_1/2 < G/4`, which fails for large `e_j`.
   * Counterexample: `M = 2`, with a column of exponent 10.
   * `v = (a_1 - a_2) + 10 e_j` lies in `(L + 10Z^r) \ L` with `||v||_1 = 2` and `G = 10 > 2||v||_1`.
   * All uses in sharp.md have `e = 1`, so they are unaffected.
10. **References.** Rosser-Schoenfeld (1962, Cor. 1), Robin (1983, `omega(n) <= 1.38402 log n/loglog n`),
    Siegel's lemma (e.g. Bombieri-Gubler, Lemma 2.9.1) and Burgess are named without precise
    references. Prop. 2.1 is standard (unit groups of `O/q^a`, the norm-one torus, Hilbert 90).
    sharp.md does not claim it as new, which is correct.

---------------------------------------------------------------------------------------------

## 7. Corrected statements (lens items only)

* **Prop. 6.5 (corrected).** Same hypotheses. With probability `-> 1`, some non-twin column `j`
  and some `s in {0, 1}` have `m_(e_j, s) > nu_s^(-1) e^(w_j/2)`. So Prop. 6.1 fails for
  `(e_j, s)`. Proof as given, restricted to the `q in Q_0` with `nu_j(q) = s`. Choose the `s` whose
  part carries `sum 2/m_q >= lambda/4`.
* **Lemma 5.1b (corrected).** As stated, under the added hypothesis that every prime of a deleted
  column exceeds `X'`. This holds in all uses, since `p_j > X`.
* **Lemma 3.3(1), demand line (corrected).**
  `log m_(xy,1) <= 10 log d + 18.72 omega(d) <= mu_M log M` for `M >= 30`, using
  `max_(d<M) omega(d) <= 1.3841 log M/loglog M`.
* **Prop. 7.4 (corrected).** Replace "every pair has `log m_xy <= C_1 log M`" by "every pair has
  `log m^gen_xy <= C_1 log M`", i.e. the maximum over all four units.
* **§7.2/§7.4 (corrected).** "Coherent designs: the a-priori bound of `lower.md` Lemma 4.1 is useless
  for `X >= N^(1/4)`. Whether some coherent design works is not decided here. Heuristically, generic
  characters collide at `X`-smooth moduli of size `X^(Theta(1))`."
* **§7.3 (corrected numbers).** `e^(3.1M)` characters, Poisson mean `log 2` per dyadic-in-`log`
  range, `kk ~ 3.1M/log M`.
* **Lemma 6.2 (scope).** "For the profile of §4, or any `e = 1` profile with `2M` distinct twin
  columns".
* **Prop. 7.1.** Add "with `v(c) != 0` (automatic when `rank A = M`)".

None of these changes Theorem 5.1, Corollary 5.2, Theorem 6.4, the constants 66 / `33 max(A,2)`,
or the `[1/33, 2/3]` bracket.

---------------------------------------------------------------------------------------------

## 8. Files (`round6/referee_arith_checks/`)

| script | checks | output |
|---|---|---|
| `a_groups.py` | Prop. 2.1 (1)-(6), exhaustive, all 184 odd `Q <= 1000`, all `(t, n)` | `a_groups_out.txt`: ALL PASSED |
| `a_cor22.py` | Cor. 2.2 (1)-(5), (2.2), windings, general exponents, exhaustive realisable sets; quartic part and quartic reciprocity | `a_cor22_out.txt`: ALL PASSED |
| `a_cluster.py` | Prop. 1.2 on 15 genuine clusters (e-char, e-all, same-level, units, Cor. 2.2(4) on genuine data) | `a_cluster_out.txt`: ALL PASSED |
| `a_assign.py` | Lemma 3.1 exact to `2^27`, `R(k)`, R-S inputs, `Lambda'` bound to `2*10^6`, Robin and its `max omega` form, demand to `M = 262144`, `m_struct` | `a_assign_out.txt`: ALL PASSED |
| `a_demand_scan.py` | demand over all `M` in `[20, 2000)` and a grid to 20,000 | `a_demand_scan_out.txt` (evidence) |
| `a_design.py` | Def. 3.2 and Lemma 3.3(1),(4) class-level, Paley `M = 108`, `X = 3M, 3M^2` | `a_design_out.txt`: ALL PASSED |
| `a_lowerS.py` | Lemma S (exhaustive support 4), Prop. 4.3 / (3.1), Lemma 4.1, slack bound, `q_0 = 43, 47` | `a_lowerS_out.txt`: ALL PASSED |
| `a_sect7.py` | Siegel (Prop. 7.1), `l^1`-ball count, Mertens factor, coherent example | `a_sect7_out.txt` |
| `rerun2/` | reruns of `check_parity.py`, `check_design.py` and `check_realise.py`, with timings | identical conclusions |

The repository `/home/user/Jarnik` and `round5/` were not modified.
