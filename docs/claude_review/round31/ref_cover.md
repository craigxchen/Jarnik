# Referee report, round 31, block "cover"

Scope: `astra_lemma/answer.md`, §1 (heights and conventions) and §2.1 (the bounded-index theorem, proved with a
degree-I isogeny cover and McKinnon–Roth). The §2 theorem as stated (uniform over all I >= d_0) also uses (2) from §2.2,
which belongs to another block. I note below where that dependence enters.

Scripts (mine, independent of Astra's): `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref_cover/`
- `vol.py`: an exact Ehrhart count for the A_4 root polytope, plus Astra's combinatorial formula.
- `gauss.py`: exact Z[i] arithmetic.
- `lift_check.py`: computes the toric height place by place, the index, an HNF basis, and the exact rational lift.
- `lattice.py`: an integer echelon index routine, cross-checked against the gcd of the minors on 2000 random cases.
- `search2.py`, `search3.py`, `out_1e10.txt`: the evidence search.

## Verdict summary

| Item | Verdict |
|---|---|
| §1: height formula, columns, intrinsic index, bound (1) | CORRECT |
| (b) B^4 = 70 and h_B(y) = W | CORRECT (reproduced exactly) |
| (c) rational lift y = eps f_P(x), 4^4 translates | CORRECT (exact check on 152 random tuples) |
| (a) McKinnon–Roth Cor 4.2 and Thm 5.1 as used | CORRECT as used (relies on the record's full-text check; one citation precision) |
| (d) X_f over Q, B_f^4 = 70d, heights, distance, reality of targets | CORRECT, with expository corrections |
| (e) pushforward of exceptional sets; finitely many lattices | CORRECT |
| (f) thresholds B_5 > 17/5, B_7 > 94/25 > 124/33; A_5(1/31) normalisation | CORRECT. The uniform form for all I >= 7 is conditional on §2.2 (2) |
| Break attempt | No counterexample (weak evidence only) |

## 1. Section 1

**Primitivity.** A primitive tuple excludes inert and ramified primes from H. If an inert q divides H, then q divides
every z_i. If 2 divides H, then (1+i) divides every z_i. At a split prime p, primitivity at pi and at pi-bar forces
min_i a_ip = 0 and max_i a_ip = e_p.

**Formula for y.** z_i/z_5 = (u_i/u_5) prod g_p^{a_ip - a_5p}. This gives y_j = eps_j prod g_p^{b_jp} with eps_j in
mu_4. Note that i is in U(Q).

**Index.** I = gcd of the 4x4 minors. Changing the conjugate choice or the anchor acts unimodularly, so I is
intrinsic. I checked numerically that |det(HNF basis)| equals the gcd of the minors on every full-rank random case.

**Height of s.** For s = eps A/A-bar, the absolute height is D(s)^{1/2}: the only contributions are p^{k/2} at the
places pi-bar_p.

**Bound (1).** |s - 1| = |eps A - A-bar|/|A| >= |A|^{-1} = D^{-1/2}, because eps A - A-bar is a nonzero Gaussian
integer. CORRECT.

**Normalisation.** Delta = diam/R and H = R^2. The arc condition of A_5(delta) is "all mutual distances
<= R^{1/2-delta}" (expert_note §3), which is exactly Delta <= H^{-(1+2delta)/4}. The common phase is removed by taking
ratios, not set to zero. CORRECT.

## 2. (b) B^4 = 70 and the height

**Volume.** `vol.py` counts the lattice points of nR_5 = {x in Z^5 : sum x = 0, sum |x| <= 2n} exactly for
n = 0..6: 1, 21, 131, 471, 1251, 2751, 5321. The interpolated Ehrhart polynomial is
35n^4/12 + 35n^3/6 + 85n^2/12 + 25n/6 + 1. This is identical to round22/toric.md Prop 1.1(c), and
4!·35/12 = 70.

**Astra's count.** The off-hyperplane sum is sum_k C(5,k) sum_s C(s-1,k-1) C(s-1,4-k). I summed it symbolically: it
equals 35n^4/12 - 65n^3/6 + ..., with the same leading coefficient. Also sum_k C(5,k) C(3,k-1) = 70.

**Height.** `lift_check.py` computes, for each place pi | p of Q(i), the spread max - min of v_pi(y_i). It does not
use Astra's formula. On all 300 random primitive tuples (2–6 primes, exponents up to 3, random units), the product of
p^{spread} over places is exactly H^2. With local weights 1/2 the absolute toric height is therefore H exactly, and the
archimedean factor is 1 because |y_i| = 1. So h_B(y) = W, in agreement with toric.md and expert_note §3
("h_B = W'"). CORRECT.

## 3. (c) The rational lift

`lift_check.py` treats the 152 full-rank cases among the 300 tuples (the other 148 have rank < 4). For each it
computes:
- the column lattice and its HNF basis P, with |det P| = I;
- h_p = P^{-1} b_p, checked to be integral;
- x_k = prod_p g_p^{h_kp}, as exact Gaussian rationals.

It then checks that y_j / prod_k x_k^{P_jk} lies in {±1, ±i} for every j. Every case passed.

So y lies in mu_4^4 · f_P(T'(Q)), and only the 256 translates eps in mu_4^4 occur. Every denominator prime is
accounted for, since x is built from the same g_p. No Gaussian gcd or denominator cost is dropped. CORRECT.

## 4. (a) The McKinnon–Roth statements

I did not reopen the PDF. I rely on the record's full-text check, `papers_check/check.md` items 4–5 and
`read_mr.md` §3.2, which quote:

- **Cor 4.2 (p. 27).** "For any ample L and x in X we have beta_x(L) >= n/(n+1)·(Vol(L)/mult_x X)^{1/n}". Here X is
  irreducible of dimension n and x lies in X(k-bar).
- **Thm 5.1 (pp. 28–29).** X is irreducible over k. A point x_v in X(k-bar) is chosen for each v in S, and L is an
  ample Q-bundle. If sum_v beta_{x_v}(L) R_v > 1, then (5.2) holds: there is a proper Z such that, for every
  delta_v > 0, d_v(x_v, y) <= H_L(y)^{-(R_v+delta_v)} has only finitely many solutions y in X(k) \ Z(k). Z depends
  on R.

**Answer to (a): yes.**
- The target may be algebraic and non-rational.
- The approximants are k-rational, here Q-rational points of X_f.
- S = {infinity} is a single fixed archimedean place, and no primes enter S. The record's caution about growing S
  therefore does not apply.
- The exceptional set is a proper closed subset, depending on (X_f, B_f, q, R).
- mult_q = 1 because q lies in the open torus, which is smooth. The bound (4/5)(70d)^{1/4} follows.

**Precision.** The "finitely many solutions outside Z" form is (5.2) of Thm 5.1, not (5.1). Astra's s with
1/B_d < s < 1/g corresponds to R + delta with beta·R > 1. Each fibre point q needs its own application of Thm 5.1,
since S has one place and therefore one target. There are I targets, so the exceptional sets form a finite union. This
is implicit in "finitely many targets", and it is harmless.

**Normalisation.** MR heights are multiplicative relative to k. For k = Q this is the absolute height, so the
convention matches.

## 5. (d) The cover

**Model.** X_0 is the toric variety of the normal fan of the root polytope. This is the record's X', with
B = -K_{X'} ample. The Galois action is -1 on cocharacters and preserves the fan and the centrally symmetric polytope,
so B is defined over Q.

**X_f.** The isogeny f_P : U^4 -> U^4 is a Q-homomorphism, because U is a Q-torus and P is an integer matrix. The
translate t = eps lies in T(Q). So psi = t f_P is defined over Q, and the normalisation X_f of X_0 in Q(T') is defined
over Q. It is geometrically integral, since T' is.

**Correction (expository).** X_f can be described explicitly. It is the toric variety of the same fan Σ, read in the
finer cocharacter lattice N' = P^{-1}N, with the twisted Galois action. Its polytope is P^T(R_5) in M', and M' has
index I over P^T M. Hence B_f^4 = 70 I, which agrees with the general fact deg(psi̅^* B)^4 = deg(psi̅)·B^4. B_f is
ample as the pullback of an ample bundle under a finite morphism. CORRECT.

**Heights.** h_{B_f} = h_B ∘ psi̅ + O(1) by functoriality. The constant depends on (P, eps) only. CORRECT.

**Real points (the question asked).** All kernel points are real.
- Over C, U(C) = C^* with real structure z -> 1/z-bar, so U(R) is the unit circle.
- The fibre psi^{-1}(e) consists of x with x^P = eps^{-1}. Taking log|·| gives P·log|x| = 0, and P is invertible over
  Q, so |x_k| = 1 for all k.
- So all I fibre points lie in U(R)^4. They are tuples of roots of unity, hence algebraic, real, and smooth points of
  multiplicity one.

**Correction (expository).** Astra's "by properness" is better put as follows. f_P restricted to (S^1)^4 is an I-sheeted
covering map of compact tori: θ -> Pθ is a local isometry up to ||P||. So a real x with eps f_P(x) = y close to e is
within O_P(Delta) of exactly one fibre point. The étale local inverse gives dist(x, q) <= C_P·max_j |y_j - 1| <= C_P·Delta.
Any fixed projective distance is locally comparable to the torus coordinates there. CORRECT.

**The chain of inequalities.** If Delta <= H^{-1/g}, then
dist(x, q) <= C·H_B(y)^{-1/g} <= C'·H_{B_f}(x)^{-1/g} < H_{B_f}(x)^{-s} once H is large, because s < 1/g. Distinct y
give distinct (x, eps). CORRECT.

## 6. (e) Pushforward and finiteness

- psi̅ is finite, so the image of a proper closed subset is proper closed, of dimension at most 3.
- Over the finitely many (L, eps, q) with I in [d_0, I_*], the union is proper closed. Each L has finitely many
  index-I sublattices (HNF), there are 256 values of eps, and I values of q.
- The finitely many solutions outside Z are absorbed into H_0. A fixed y determines the primitive tuple up to a unit,
  and H = H_B(y).

So Z and H_0 depend only on (d_0, g, I_*(g)). CORRECT.

**Scope remark.** The hypothesis rank L = 4 is essential. Rank < 4 tuples lie in infinitely many proper subtori,
whose union is Zariski dense, so they are not absorbed by Z. §2.3 handles them separately; that is outside this block.

## 7. (f) Numerics and normalisation

- **B_5.** (4/5)·350^{1/4} = 3.4602. It exceeds 17/5 because 4.25^4 = 326.25 < 350, and 17/5 > 10/3.
- **B_7.** (4/5)·490^{1/4} = 3.7643. It exceeds 94/25 = 3.76 because 4.7^4 = 487.97 < 490, and 94/25 > 124/33 because
  94·33 = 3102 > 3100 = 124·25.
- **B_6.** (4/5)·420^{1/4} = 3.6216 < 15/4. So index 6 is not covered by this mechanism, and §2.1 has nothing to say
  about indices 1–6.
- **A_5(1/31).** A_5(1/31) requires Delta > H^{-(1+2/31)/4} = H^{-33/124}, i.e. g = 124/33 = 3.7576. Note
  124/33 > 15/4, so delta = 1/31 < 1/30 is genuinely inside the window the record needs.

**Conditionality.** For any fixed finite range d_0 <= I <= I_1 the theorem follows from §2.1 alone. For the uniform
statement over all I >= 7, the tail I > I_* uses (2) of §2.2, a Baker–Wüstholz bound that I did not check. If (2)
holds as stated, the conclusion "A_5(1/31) for full-rank tuples with I >= 7, off one proper closed Z" follows.

## 8. Attempt to break it

**Structural.** A five-tuple with columns in an index-7 lattice is exactly y = eps·f_P(x) with x in T'(Q). A
Zariski-dense family that approximates e with exponent below B_7 would give a Zariski-dense family upstairs violating
MR Thm 5.1. So a counterexample would have to refute McKinnon–Roth. Finitely many examples, or examples inside a
proper closed set, refute nothing.

Simple constructions are far from the threshold.
- Raising a tuple to a power multiplies g by about the degree.
- Disjoint-prime product lifts (x_1 near zeta_7, the other coordinates independent) give H = prod D(y_j) and
  Delta >= max D^{-1/2}, i.e. g about 8.
- 1-PS conic lifts have rank 1, outside the hypothesis.

**Exhaustive window search (evidence only).** `search3.py` took N <= 1e10 composed of the split primes 5–73 with
r_2(N)/4 >= 96. It enumerated all 5-subsets of 8 consecutive points, reduced each by its Gaussian gcd (H' >= 1e5), and
computed I and g_eff = log H'/(-log Delta).

Minimum g_eff by index:

| Index | Min g_eff |
|---|---|
| I = 1 | 4.71 |
| I = 2 | 5.92 |
| I = 3 | 5.57 |
| I = 4 | 5.16 |
| I = 5 | 5.61 |
| I = 10 | 5.93 |
| I = 12 | 5.12 |
| rank 3 | 5.13 |

No tuple with I >= 7 and g_eff < 6 was found. At heights this small, exponents are dominated by constants and say
nothing about the asymptotic 3.76. The search is consistent with the theorem but is not a test of it.

## 9. Novelty

Grep finds no prior use in the record of valuation-column sublattices or of isogeny lifts.
- The record's only cover discussion is MR Cor 8.9 on étale covers, dismissed in papers_check because X is simply
  connected.
- round25/coset.md uses isogeny charts for a different purpose.

Astra's point is that the relevant covers ramify on the boundary and that only a thin set (Serre type II) of
rational points lifts. This sidesteps that dismissal. It is new to the record. To a specialist it is a routine
application of McKinnon–Roth on a finite cover: "points in the image of a degree-d cover inherit beta scaled by
d^{1/n}", the higher-dimensional analogue of y = x^d on P^1.

## 10. Significance

Peripheral for the goal.
- It covers only full-rank tuples whose valuation columns lie in a proper sublattice of index >= 7 (>= 5 to beat 10/3).
  Generic many-prime tuples have index 1; for example, any tuple with four primes carrying the singleton cuts has
  index 1.
- It does not touch indices 1–6, and index 1 is the core of Question B.
- It does not change the growth rate, the uniform bound, Question A or route 5.
- What it does show is an unconditional essential exponent above 15/4 at e, but on a Zariski-dense thin set. That
  refines, without contradicting, the record's remark that no unconditional alpha^ess > tau is known at a rational
  point (that remark concerns all rational points off a closed set).
- Its value for (L_5°) rests on the reduction to indices 1, 3, 5 (§§3–4, other blocks). Even granted, the hard case,
  index 1, is untouched.
