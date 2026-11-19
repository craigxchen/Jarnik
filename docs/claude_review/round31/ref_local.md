# Referee report, round 31, block "local": astra_lemma/answer.md §5

Scope: the claim that no nondegenerate, ray-generic member of the round-24 reflection-symmetric family (round24/search.md §5.1; referee_search.md §§3-4) descends to the twisted torus over Q. That covers the descent to two square conditions, the parameter recovery, the discriminant identity, the 2-adic table, and the split comparison at 3.

My own scripts are in `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_local/`:
- `chk1.py`: sympy identities.
- `chk2.py`, `chk3.py`, `chk4.py`: exact residue analysis at 2 and 3 for both forms.
- `chk5.py`: small-height search for the square conditions.

I also ran Astra's `check_component.py`. It passes.

## Verdict summary

| claim | verdict |
|---|---|
| Definitions of c, D, u, v, G agree with round 24 | CORRECT |
| v/u = 4(m^2-1)/(m^2-9) | CORRECT (sympy) |
| u(m)-u(-m) = -256 m^3 / ((m^2-9)^2 (m^2-1)^2) | CORRECT (sympy) |
| disc G = u (2D) [2(m+1)c/((m-3)^5(m-1)^3)]^2 | CORRECT (generic monic quartic discriminant substituted, exact) |
| Descent gives: u is a nonzero rational square, and G splits into distinct rational linear factors | CORRECT (checked line by line; one small step should be made explicit, see 2(c)) |
| m is rational | CORRECT |
| 2-adic table mod 256 and the conclusion of emptiness (twisted) | CORRECT (all odd residues mod 2^12 checked exactly, plus the even and non-integral cases) |
| Split form: -A, -B squares; 3-adic exclusion | CORRECT |
| "The two forms have different local obstructions" | CORRECT WITH CORRECTIONS. Each form is in fact obstructed at both 2 and 3 (see §4). |

Overall: **CORRECT WITH CORRECTIONS** (the corrections are minor and do not touch the theorem). Novelty: new; it settles Q24.2 negatively for this component. Significance: peripheral.

## 1. Identities (independent sympy, `chk1.py`)

**Definitions.** c, D, u, v and G match round24/search.md §5.1 and referee §4.1. Here D = m^4+70m^2-48m+9 is the odd-order quartic factor of disc_Y G, not the `D` of `L = Y(Y+e)^2(Y+D)`. The notation clash is harmless.

**Parameter identities.**
- `v/u` factors as `4(m-1)(m+1)/((m-3)(m+3))`.
- `u(m)-u(-m)` factors as `-256 m^3/((m-3)^2(m-1)^2(m+1)^2(m+3)^2)`.
- Limits at infinity: u → -1/16 and v → -1/4.

**Discriminant identity.** I substituted the coefficients of G into the generic discriminant of `Y^4+a3Y^3+a2Y^2+a1Y+a0`. The result equals `u·2D·[2(m+1)c/((m-3)^5(m-1)^3)]^2` identically. In factored form it is `-(m+1)^2(m+3)c^3 D / (2(m-3)^12(m-1)^8)`. This is consistent with the parities in the record.

**Root relations.** With `Y = -b^2/v`, the expression `Y^2-Y+u-b` is the b-quartic divided by v^2, and `G(-b^2/v)` vanishes modulo the b-quartic. So the Y_r are the roots of G, as claimed.

**Irreducibility.** D and c are irreducible over Q, so neither has a rational zero.

**Recovering m^2.** From v/u = k one gets m^2 = (4-9k)/(4-k). The value k = 4 is impossible, so m^2 is rational whenever u and v are.

## 2. The descent argument (checked line by line)

**(a) Rational affine coordinate on the X-line.**
- Forgetting u_4 is a Q-morphism of twisted tori, so the projected image is a Q-curve.
- On the X-line, u_1, u_2, u_3 are quotients of monic quartics, P_j over its formal partner. Their 14 places carry the distinct primitive vectors ±1_T, T a nonempty subset of [3].
- Pullback multiplies valuation vectors by the ramification index. Primitivity therefore forces index 1 above each such place, and uniqueness of each vector forces degree 1. So the X-line is the normalization of the projected image, hence a Q-curve.
- Im P_j is constant (record, round 24 §5.1(e) and referee c1), so u_j - 1 = const/Pbar_j vanishes only at X = ∞.
- Hence ∞ is the unique point over the projected identity. It is Galois-fixed, so rational, and X = λ(x - x_0).

Correct.

**(b) Places and conjugation.** The same argument makes the parametrizing conic the normalization of the full curve. Places are Q(i)-rational because Galois elements fixing i preserve the distinct ray vectors. Complex conjugation sends 1_A to 1_{A^c}. In the normal form the partner of ρ_S = a + ib is the formal conjugate a - ib. This is the definition of Pbar_j for complex members, consistent with the referee's repair in round 24, §4.2. Correct.

**(c) x_0 is rational.** Astra averages the 16 distinct projected values. In X these sum to 0, since Σ(±a_r ± i b_r) = 0.
- Distinctness is not argued in the text, but it follows from the hypotheses. Over each ρ_S (S a subset of [3]) lie the two conic places S and S∪{4}. A coincidence of ρ_4's X-value with another value would put three places in a degree-2 fibre.
- Alternatively, split the 16 values into two Galois-stable subsets: the 14 boundary values of the projected curve, and the pair {ρ_4, partner}. Their sums are 14x_0 + 2a_0/λ and 2x_0 - 2a_0/λ, and both are rational.
- Caution: averaging the multiset of the 30 place images does *not* give x_0. In X that multiset sums to 2a_0, not 0.

Correct, with this step to be made explicit.

**(d) The ratios a_r/λ and b_r/λ are rational.** We have x_0 + (a ± ib)/λ ∈ Q(i), and the two values are complex conjugates. Their half-sum gives a/λ ∈ Q. Their half-difference divided by i gives b/λ ∈ Q. This holds for complex a_r, b_r. Correct.

**(e) u is a nonzero square and G splits.**
- Since e_4 = t^2 = Π b_r, u = e_4/e_2^2 = ((t/λ^2)/(e_2/λ^2))^2, with t/λ^2 = (a_r/λ)(b_r/λ) ∈ Q.
- e_2/λ^2 = -½ Σ(b_r/λ)^2 is nonzero, being minus half a sum of rational squares that are not all zero.
- v = (e_3/λ^3)^2/(e_2/λ^2)^3 is rational. It is nonzero, since e_3 = 0 would repeat a b_r^2.
- The Y_r = -(b_r/λ)^2/(e_2/λ^2) are rational and distinct, and they are the roots of G.

Correct.

**(f) m is rational.** By §1, m^2 is rational. If m is not rational and m ≠ 0, then m and -m are conjugate, so u(-m) = u(m). This contradicts u(m)-u(-m) ≠ 0. The points m = -3, -1 and the zeros of c give uv = 0, and m = 1, 3 are poles. At m = ∞, u = -1/16 is not a square. Correct.

**From the square conditions to A and B.** u = A/(16(m-3)^2(m-1)^2), so A = -(m+3)c must be a square. The discriminant of a split polynomial is a square (Vandermonde). Its bracket factor is nonzero, since m ≠ -1 and c ≠ 0, so B = 2D is a square too. Correct.

## 3. 2-adic step (`chk2.py`, exact)

**Odd m.** A(m) and B(m) mod 256 depend only on m mod 8. I checked all odd residues mod 256:

| m mod 8 | A mod 256 | B mod 256 |
|---:|---:|---:|
| 1 | 192 | 64 |
| 3 | 64 | 128 |
| 5 | 0 | 192 |
| 7 | 128 | 0 |

This reproduces Astra's table. In every row one entry is a certified non-square: 128 has v_2 = 7, and 192 = 64·3 with 3 ≡ 3 mod 4. I repeated the check at precision 2^12 over all odd residues.

**Other m.**
- m even in Z_2: D is odd, so v_2(B) = 1.
- v_2(m) < 0: v_2(B) = 1 + 4v_2(m), which is odd. I verified that the numerator is odd.

Emptiness for the twisted form: **CORRECT**.

As a side check, a search over m = p/q with |p| ≤ 600 and q < 200 (`chk5.py`) finds no m with A a nonzero square. This is consistent with the result.

## 4. Split comparison at 3, and a correction

**The split claims.**
- -u square ⟺ -A square.
- disc G = (-u)(-B)·square ⟹ -B square.
- m ≡ 1, 2 mod 3: -B ≡ 2 mod 3.
- m ≡ 0, 3 mod 9: -A ≡ 6, 3 mod 9 respectively.
- m ≡ 6 mod 9: -B ≡ 54 mod 81, so v_3 = 3.
- v_3(m) < 0 or m = ∞: u ≡ v ≡ -1 mod 3, and G reduces to Y^4+Y^3+2Y^2+Y+1. Its only root in F_3 is 1, and it is not (Y-1)^4 = Y^4+2Y^3+2Y+1 mod 3.

All verified exactly over every residue mod 3^6. **CORRECT.**

**Correction to "the two forms have different local obstructions".** True as a description of the obstructions Astra exhibits, but misleading as a contrast: **each form is obstructed at both 2 and 3** (`chk2.py`, `chk3.py`, `chk4.py`, exact).

*Twisted form at 3 (needs A, B squares):*
- m ≡ 1, 2 mod 3: A ≡ 2 mod 3.
- m ≡ 0, 3 mod 9: v_3(A) = 1.
- m ≡ 6 mod 9: B ≡ 27 mod 81, so v_3(B) = 3.
- v_3(m) < 0: u ≡ -1 mod 3, which is not a square.

So the twisted family is already empty at 3 alone.

*Split form at 2 (needs -A, -B squares):*
- m even or v_2(m) < 0: v_2(B) is odd.
- m ≡ 1 mod 8: -B ≡ 192 mod 256.
- m ≡ 3 mod 8: -B ≡ 128 mod 256.
- m ≡ 7 mod 8: -A ≡ 128 mod 256.
- m ≡ 5 mod 8: -B = 64·(odd ≡ 5 mod 8). For example, at m = 5, -B = -64·67 with -67 ≡ 5 mod 8. This case needs precision 2^9, beyond 256; I checked it to 2^14.

So the split family is already empty at 2 alone.

Corrected statement: *"Both the twisted and the split Q-forms of this family have no nondegenerate members. Each is excluded by its square conditions at 2 and, independently, at 3."* This reinforces round 24's referee finding (§5.2) that the elimination of this component is not an anisotropy effect. The split/twisted distinction plays no role here.

## 5. Standing cautions

- No fixed-S Subspace theorem is used.
- No Gaussian gcd or denominator is discarded: the argument is exact algebra on one curve.
- No drifting phase is set to zero.
- No uniform threshold is claimed.

No violation.

## 6. Novelty

The record (round 24, search.md Theorem B and referee §4) proves finiteness of Q-members of M_C by monodromy, genus and Faltings, ineffectively, in both forms. The round-24 README (line 345) leaves "some Q-member of M exists" open as Q24.2, estimated at 20-30%. A grep for 2-adic or square-class obstructions to M finds nothing.

Astra's result is new. It is elementary and effective: it upgrades finiteness to emptiness for nondegenerate, ray-generic members of the parametrized family and its S_5-images, and so settles Q24.2 negatively. It reuses the record's descent (referee §4.2) almost verbatim. The new ingredients are the square class of disc G, the requirement that u be a square, and two local computations.

## 7. Significance

Peripheral. It changes neither the growth rate (2/3+ε) log R/loglog R nor any uniform bound. It does not change the status of Question A, Question B or route 5. The round-24 README already rates M's rational members as irrelevant to Question B: they could force at most finitely many curves into Z. Other components of the class (Q24.1), chain-place degenerations, and sporadic tuples are untouched, as Astra itself says in §6.

The only map value is that one component of the 15/4 twisted Pell class now has no Q-points at all, rather than finitely many unknown ones. And since the split form dies in the same way, there is still no anisotropy signal.
