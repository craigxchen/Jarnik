# Round 31, Part III referee report: block "quartet"

Scope: `astra_free_III/answer.md`. The block covers the opening claims of Part III about the Paley double core (lines 658-662), sections III.1, III.2 (including "Exact tests"), III.3, and the planted fourteen-prime example the brief asks about. That example is printed in III.7 (line 898) and computed in `iii_quartet_extension.py`.

Scripts are in `round31/ref3_quartet/`. All are my own code, written independently of Astra's scripts.
- `genuine.py`: the identity on genuine clusters from `growth/families.md`.
- `genuine_batch.py`: 3077 genuine same-unit short-arc quartets.
- `exhaust.py`: my own kernel-lattice solver, used to reproduce both fourteen-prime computations.
- `paley.py`: the Paley columns built from lw.md §5.1, class sizes, the five-mark structure, the (III.3) lcm, and the 32-pattern count.

I also ran Astra's `iii_quartet.py` table.

## Verdicts in one line each

| Claim | Verdict |
|---|---|
| III.1 three-product identity, coprimality, `h3-h2` rule, `ABC | N`, lcm formula for `(A,B,C)` from the `d_ij` | CORRECT (proof checked; exact check on genuine clusters). NOT NOVEL: it is `growth/ptolemy.md` Theorem A(2),(4). |
| III.1 lattice test and recovery lemma (III.2) | CORRECT (proof checked line by line). |
| III.2 assertion 1 (a random assignment passes even one quartet with prob. `O(m^4 2^{-cm})`) | CORRECT; `c = 1/32` effectively. |
| III.2 assertion 2 (admissible assignments are pairwise at distance `> floor(r/64)`, hence at most `r!/floor(r/64)!` of them) | CORRECT. |
| III.2 assertion 3 (at most one of the `2^m` within-pair orderings is admissible) | CORRECT. |
| "Actual clusters on this profile must be admissible" | CORRECT under the one-unit convention and `Delta <= N^{-1/4}/2`. |
| "Every assignment still admits an angle fake" | CORRECT. The proof of round 8 Thm 6.1 is assignment-blind. The record already states it: `round16/barrier.md` Thm 2.3, "for every assignment of primes in `[P, P e^{1/(4r)}]`". NOT NOVEL. |
| Hasse count: each matching class has `r/8 + O(sqrt m)` columns | CORRECT (checked to `q = 79`). |
| Table `m = 7, 11, 19`: 0 surviving quartets | CORRECT (Astra's script rerun). |
| First parent (100049..100333): 156/45045 pass for k=2, 0/560560 for k=3; 16/12 factorable | CORRECT (reproduced exactly with own solver). |
| Planted parent: exactly 1 of 560,560 passes, relation `-1145443 A + 1927968 B - 782525 C = 0`; no five marks pass on any assignment | CORRECT (reproduced; the five-mark claim reduces to a combinatorial fact I checked). |
| III.3 (III.3) `lcm_k Q_k | B(ij|xy)` and elimination of a common mark | CORRECT. It is a negative result. |
| III.3 32-pattern count and log sizes `W/16`, `W/8` | CORRECT. |
| Significance | PERIPHERAL. No change to the growth rate, the constant 2/3, the uniform bound, Question A, Question B or route 5. |

## 1. III.1: the identity and the test

**Identity.** Put `[P_i,P_j] = g_ij t_ij` with `g_ij = prod p^{min}`. This works because `conj(P_i) P_j / p^{min}` is integral at every `p`. The Plücker relation `[12][34] - [13][24] + [14][23] = 0` then gives `g12g34 x - g13g24 y + g14g23 z = 0`.

At a prime, sort the exponents `h1 <= h2 <= h3 <= h4`. The low/high pairing has exponent `h1+h3`; the other two pairings have `h1+h2`. Hence the gcd has exponent `h1+h2`, and only `h3-h2` survives, on one matching. This gives pairwise coprimality and `ABC | N`.

Lcm formula. `D_k = X/A_k^2` with the same `X`, and the `A_k` are coprime, so `lcm(D_1, D_2, D_3) = X` and `A_k = sqrt(L/D_k)`. Correct.

Cofactor bound. `|t_ij| = sqrt(d_ij) |sin(star angle)|`, and the star angle is half the point angle. So `|t_ij| <= sqrt(d_ij) Delta/2`. This needs all points in one Gaussian unit class (the shortened one-unit convention, a pigeonhole loss of 4). Across unit classes the star angle shifts by multiples of `pi/4` and the bound fails. Astra states the convention, so this is not an error.

**Exact checks** (`genuine.py`). I computed all lattice points, took the best window and divided by the Gaussian gcd. Then, from the actual Gaussian factorisation, I computed `P_i`, `g_ij` and `t_ij` and checked four things: the signed identity, the `h3-h2` matching rule, the lcm formula, and `ABC | N`. The circles tested:
- 8-point record `R^2 = 32988780417125 = 5^3 13^2 17^2 89·109·557` (all 70 quartets; prime powers present);
- 7-point record `235370525` (35);
- golden `n=4`, `537606725` (70);
- `683617125425`, `412333805`, `10414625` (5, 1, 5).

All pass. On the same-unit quartets, `max |t|/(sqrt(d)·angle/2) = 1.000` (equality up to `sin x <= x`).

`genuine_batch.py`: 105 random circles with only split primes, 3077 same-unit quartets of span `<= 60 n^{-1/4}`. Every one satisfies:
- (III.1);
- `|t_ij| <= sqrt(d_ij)·span/2`;
- `A x + C z ≡ 0 (mod B)`.

So the genuine quartets pass the shortest-vector test. Same-unit quartets on arcs with `C <= 8` did not occur in 400 small circles, which is consistent with the scan data in families.md §6.

**Lattice test.** `gcd(A,B) = 1`, so the solutions are `x ≡ -C A^{-1} z (mod B)`. The lattice basis `(B,0), (r,1)` is correct. A nonzero `(x,z)` with `|x|,|z| <= U` has squared length `<= 2U^2`. Correct.

**Recovery lemma (III.2).** Write `AC'/(A'C) = u/v` in lowest terms. The common factor `k` is coprime to `Q` because `A`, `C`, `A'`, `C'` are coprime to `B` and `B'` respectively. Mod `Q` the two congruences give `k(u x z' - v x' z) ≡ 0`, so `Q | u x z' - v x' z`. The absolute value is at most `2U^2 max(u,v)`, so under (III.2) the congruence is an equality. Then `gcd(u,v) = 1` puts every prime of `uv` in `x x' z z'`. If all coefficients are coprime to the parent, then `u = v = 1`, so `A = A'` and `C = C'`. Correct.

**Novelty of III.1.** The identity is `growth/ptolemy.md` Theorem A(2),(4): `X_k` = product of 2|2 layer norms, with `c X_2 = a X_1 + b X_3`. Theorem B there identifies it with the Plücker relations. Astra acknowledges this ("ordinary Plücker relation"). New in III.1:
- the measured-divisor lcm formula (minor);
- the 2D-lattice recovery lemma (III.2).

## 2. III.2: the fixed-parent theorem

**Setting checked.**
- Ordinary rows of the lw.md §5.1 double core: `s1_a(x) = [x=a] + chi(a-x)` and `s2_a = s1_a` negated at `x = a` (border row deleted). My construction from the lw.md definition matches Astra's `ii_paley.py`.
- Every pair of ordinary rows is separated by exactly `m+1` of the `2m` columns: identity `sum_a chi((a-i)(a-j)) = -1`, plus one separating copy at `a = i` and one at `a = j`. Verified for `q` up to 79.
- The `2m` columns are distinct up to complement. Verified.
- Then `d_ij <= (P e^{1/(10r)})^{m+1}` and `sqrt N >= P^m`, so `t^2 <= d/(16 sqrt N) < P/8` and `|t t'| <= floor(P/4) = U`.

So genuine clusters on this profile (squarefree, narrow bin, one unit class, `Delta <= N^{-1/4}/2`) are admissible. Correct.

**Hasse class sizes.** The indicator formula `(1/8)(1 + sum eps_i eps_j chi chi + (prod eps) chi(prod(a-i)))` is correct: the odd terms cancel under the `±` symmetry. The quartic is monic and squarefree, so the curve has genus 1 with two rational points at infinity. Hence `|sum chi(f) + 1| <= 2 sqrt q`, and the stated `2 sqrt m + 1` is right (Weil's `3 sqrt q` would also do). Numerically the minimum class size over all quartets is `>= r/16` already for every `q` in 7..79:
- `q = 79`: 16 vs `r/16 = 9.9`;
- `q = 7`: 2 vs 0.875.

**Assertion 1.** Pair `s = floor(r/64)` A-slots with C-slots. Swaps within the pairs fix `B`, so `Q = B >= P^{r/16}`, and `u/v = prod (a_k/c_k)^{±2}` over the swapped pairs gives `max(u,v) <= (2P)^{2s}`. Since `|coeff| <= U < P`, the coefficients are coprime to `N`. So all passing members of an orbit share `A`, and `A` fixes the orbit element. The fraction is `<= 2^{-s}`; a union bound over `C(m,4)` quartets gives `O(m^4 2^{-floor(m/32)})`. Correct. This bound is `< 1` only once `m` is in the low thousands; "sufficiently large" is honest.

**Assertion 2.** Let `h <= r/64` primes move. The unmoved B-class primes give `Q >= P^{r/16-h}`, and `max(u,v) <= (2P)^{2h}`; then `3h < r/16 - 2` for `r > 128` gives (III.2).

Distinguishing quartet: orient `S, S'` to agree on at least half the marks. On the agreement set, `sum_{Ag} S = (sum S + sum S')/2 = O(1)`, so the common ones and common zeros are each `>= m/4 - O(1)`. That is exactly the step stated. A disagreement mark plus one common one and two common zeros makes `S` 2|2 and `S'` 1|3. Then the moved prime divides the matching product in one assignment and not the other, a contradiction.

The packing count is right: assignments with the same placement of the first `r-s` primes differ on `<= s` primes, so there are at most `r!/s!` admissible assignments. Correct.

**Assertion 3.** A within-pair swap at base `a` changes the quartet type only when `a` is one of the four marks, so at most 8 primes change type per quartet and `Q >= P^{r/16-8}`. The distinguishing quartet is `{a, y1, y2, y3}` with `chi(a-y1) = +1` and `chi(a-y2) = chi(a-y3) = -1`. Correct.

**Angle fakes.** In the proof of lw.md Thm 6.1 the primes enter only through `w_min`, the bin width (narrower is fine) and `Z(s)`. The margins come from `Phi(c)`, which depends on the profile only. So it is assignment-independent.

Bins: `[m^D, 2m^D]` has `~ m^D/(2D log m)` primes `≡ 1 (mod 4)` and `~14m` bins of width `1/(10r)`, so `D >= 3` suffices for `2m` primes in one bin. For the gap-only fake `P >= (40M/C_0)^4` is needed (`D >= 5` or so); for `LW(eps)`, `D >= 2/eps + O(1)`. Deleting the border row only discards a point.

The statement is already in the record: `round16/barrier.md` Thm 2.3, "for every assignment of primes in `[P, P e^{1/(4r)}]` there is a gap-admissible angle model". Correct, not new.

## 3. The computed examples (reproduced independently)

`exhaust.py`:
1. Kernel basis `(B,-A,0)`, `(-Ca,-Cb,1)` with `aA + bB = 1`.
2. Exact Lagrange reduction.
3. Exit if `|v|^2 > 3U^2`.
4. Otherwise bound `|beta| <= sqrt(3U^2 |v|^2/D)` and enumerate `alpha` from the three sup-norm interval constraints.

The `U` values (`U = ceil(prod of 8 largest/(16 isqrt N))`):
- first parent: 6294;
- planted parent: 62500076862.

| Parent | k per class | triples | passing | factorable (indiv / even) |
|---|---|---|---|---|
| 14 split primes 100049..100333 | 2 | 45,045 | 156 | 16 / 12 |
| same | 3 | 560,560 | 0 | 0 / 0 |
| planted (primes ~1.000001e12, log width 1.7e-7) | 3 | 560,560 | 1 | 1 / 1 |

The planted passing partition is exactly the three planted triples, with primitive solution `(-1145443, 1927968, -782525)`, matching Astra. In every passing case the determinant exceeds `3U^2`, so every solution is a multiple of the primitive one. That justifies testing only the primitive triple.

**"No five marks pass on any assignment"** (planted parent, `m = 7`). I checked that each of the 21 five-subsets contains exactly two three-per-class quartets, and that their unordered column partitions are always different (0 of 21 equal). The quartet histogram is 21 two-per-class and 14 three-per-class, as stated.

There is a unique passing prime partition, and an assignment is a bijection. So both quartets of any five-set would have to be mapped onto that same partition, which is impossible. This is equivalent to Astra's overlap-fingerprint check, so the claim is CORRECT.

Remark: the statement is that strong only because the passing partition was planted to be unique. It says nothing about parents where several partitions pass.

Astra's `iii_quartet.py` table (`m = 7, 11, 19`, seed 2026): 0 of 35, 330 and 3876 quartets survive. Reproduced.

## 4. III.3

(III.3) squarefree: a prime in `Q_k` puts `i,j` against `k,x,y`, so it separates `ij|xy`. The prime-power version (nested thresholds) I checked on 20,000 random exponent configurations (`m` = 5..8, exponents up to 6):
- 0 violations of `max_k min(v_p B(ij|kx), v_p B(ij|ky)) <= v_p B(ij|xy)`;
- the left side is positive in 3161 cases;
- equality holds in 2874 of them.

Telescoping of the raw cross quotient: trivial algebra, correct. The unit claim is correct for genuine data: a crossed cofactor is prime to `p` when `p` separates the pair. If `w = pi^k u` with `\bar pi ∤ w`, then `p ∤ Im w`.

The 32-pattern count is reproduced exactly: `(0,0)` 14 times, every other listed pair 2 times. So the common `B` has mass `2/32 = W/16`. The reduced numerator of `AC'/(A'C)` collects `(A,0)` and `(0,C)`, giving `4/32 = W/8`. Correct.

This is an honest negative statement. It covers only this particular elimination; it is not a proof that multiplicity cannot help.

## 5. Why the saving is in assignments and not in m

Fix the parent `N` (a set of `r = 2m` primes in one narrow bin) and the profile (the Paley double core). Assertion 2 says the set of admissible bijections from primes to cuts is a code in `S_r` with minimum Hamming distance `> floor(r/64)`. The packing bound then gives `|admissible| <= r!/floor(r/64)!`, but nothing in the argument shows the code is empty.

A growth counterexample needs only one cluster: one parent, one profile and one assignment. A bound on the number of assignments therefore implies nothing about `m`. Only "admissible set = ∅ for every parent, profile and assignment" would bound `m`, and Astra correctly says the argument cannot reach it ("A growth counterexample needs only one of those assignments").

There are two further restrictions:
1. **Profile.** The result is specific to the Paley double core. A counterexample may use any profile.
2. **Narrow bin.** The result uses the narrow bin in an essential way, through:
   - `|coefficients| < least prime` (coprimality in the recovery lemma);
   - `Q >= P^{r/16}`.

   Parents relevant to the `log R/loglog R` regime have primes spread over `[1, m^{O(1)}]` and may have prime powers. There, small parent primes can divide the cofactors, and recovery holds only away from `x x' z z'`.

The probabilistic assertion 1 likewise says only that a *random* assignment fails. The adversary chooses the assignment.

## 6. Novelty and significance

- **III.1 identity:** in the record (ptolemy.md Thm A/B).
- **Angle fakes for every assignment:** in the record (round16 Thm 2.3; round 8 Thm 6.1 proof).
- **New:**
  - the recovery lemma (III.2);
  - the observation that, at a fixed parent and profile, the assignments that admit integer stars with small cofactors form a code of linear minimum distance;
  - the `2^m` within-pair uniqueness;
  - the fourteen-prime computations;
  - the (III.3) no-amplification lemma.

Grep finds no assignment-counting statement of this kind in the record.

- **Significance:** peripheral. It moves:
  - not the growth rate;
  - not the constant 2/3, nor the 3 in `W >= 3 m log m`;
  - not Question A (k = 8 lemma) or Question B;
  - not route 5.

  The "integer test that angle fakes fail" is the familiar point that fakes carry no integer star matrix. The record already knows that integer Plücker data are not fake-blind but are not sufficient at `k <= 6` either (ptolemy.md Thm C: bounded-residue integer families at `k = 5, 6`). Astra's own III.11 also builds integer five-mark models.

  The planted `m = 7` example is a well-built control. It shows that local passing does not extend to five marks on that parent, but `m = 7` and a uniquely planted partition say nothing asymptotic.

**Corrections to wording, no change to the mathematics:**
- (i) The opening sentence "a constraint that the usual angle fakes fail" should say: for all but a `1/floor(r/64)!` fraction of assignments, no integer matrix with the prescribed divisibilities and cofactor bounds exists. Fakes fail it vacuously because they have no integer matrix.
- (ii) The `c` in assertion 1 can be stated as `1/32` (`s = floor(m/32)`).
- (iii) The angle-fake comparison should cite round16/barrier.md Thm 2.3, which states the every-assignment version directly.
