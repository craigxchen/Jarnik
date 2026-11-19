# Referee 7L, block "search": exact and numerical data on H_5(D), five marks at contact above D/2

Scripts and outputs: `/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/`.
I installed `msolve` 0.6.5 and `Singular` 4.3.2 from Ubuntu packages. Everything is my own code. Labels:
**exact** (Groebner basis over Q, or a finite combinatorial enumeration), **exact mod p** (Groebner basis over
F_32003), **verified** (an isolated root refined to 60–100 digits; the Jacobian is nonsingular and the residual
is 1e-60 or less; Kantorovich-type, not interval-certified), **evidence** (random-start Newton).

## 0. Verdict

1. **Sporadic five-point tuples with 2h > D exist in characteristic 0, in the genuine conjugation model.**
   At D = 11, h = 6 and at D = 13, h = 7 I found real configurations (s_k = conj(rho_k)). They are distinct, primitive
   and squarefree, and they are verified. 4 of the 18 cut types at D = 11 have one, and 16 of the 24 at D = 13. The
   explicit D = 11 example is in §3.3. So M(11,0) >= 5 and M(13,0) >= 5 over C(t) with conjugation. These are new to the record:
   its odd-D data stop at "M(9,0) = 4 (evidence)". In Astra's normalisation h/D = 6/11 = 1/2 + 1/22 and
   7/13 = 1/2 + 1/26, so both exceed 1/2 + 1/30. In bridge.md's exponent, alpha = 2D/h = 11/3 and 26/7, both
   below 15/4.
2. **They are rigid, and they do not threaten non-density.** Expected dimension modulo affine maps is 0, and the
   Jacobians are nonsingular, so each is an isolated point and its ratio curve is one curve. Composition t -> t^k
   keeps the image. For fixed delta, the configurations with h = (D+1)/2 satisfy h >= (1/2+delta)D only when
   D <= 1/(2 delta). If each D carries finitely many, this is finitely many curves, and they lie in a proper closed
   subset. I found nothing that moves in a family above D/2.
3. **Exact negative results.** H_5(D) <= floor(D/2) for D <= 8 and D = 12, by Plotkin plus an exhaustive enumeration
   of cut types. At D = 7 every type (2 squarefree, 10 with multiplicities) is empty over C: Groebner bases over Q,
   with all coincidences saturated. **So H_5(7) = 3 exactly, and M(7,0) = 4**: the record had this only as evidence
   from 3000 numerical starts. At D = 9 the 18 non-squarefree types are empty (exact mod 32003). The unique
   squarefree type, K5 minus an edge, has **non-real** complex solutions: 218 distinct ones found, all with
   |s + conj(rho)| >= 0.48. No real solution came from 3000 real starts (evidence). So in the conjugation model
   H_5(9) = 4 (evidence); in the split model H_5(9) = 5 (verified).
4. **Above (D+1)/2 nothing exists in the data.** Every type with 2h > D other than odd D at h = (D+1)/2 has
   negative expected dimension. They are tested at D = 10 (K5: 1500 complex starts), D = 12 (no type exists) and
   D = 14 (8 types, 300 complex starts each). No nondegenerate solution was found, even over C (evidence).
5. **Answer to (4).** No configuration with h/D bounded above 1/2 as the primitive degree grows was found. What
   exists is a countable supply of rigid curves with h = (D+1)/2, at least at D = 11 and 13 (D = 15 inconclusive, §3.5), whose excess h/D - 1/2 = 1/(2D) tends to 0. Nothing here decides the unrestricted question.
   Nothing contradicts either answer. The data are consistent with the dimension-count picture: above D/2 there are
   only rigid sporadic curves.

## 1. Model, types and the dimension count (task item 3)

Squarefree split form: roots rho_k, s_k (k = 1..n). Point i is z_i = prod_k (t-rho_k)^{a_k^i} (t-s_k)^{e_k-a_k^i},
with primitive columns (min 0, max e_k). Contact h is equivalent, by Newton in characteristic 0, to
sum_k (a_k^i - a_k^0)(rho_k^m - s_k^m) = 0 for m = 1..h-1 and i = 1..4. That is 4(h-1) equations at rank 4.
The conjugation model is s_k = conj(rho_k). For distinct columns I proved this is the only possible real structure:
conj(rho_k) = s_l forces column l = column k, and conj(rho_k) = rho_l forces complementary columns, which are
identified. Gauge: the affine group (2 dimensions).

**Expected dimension modulo affine maps:** e = 2n - 4h + 2 (rank 4). All enumerated types have rank 4. Swept
dimension in T is at most e + 1, so density needs e >= 3, i.e. 4h + 1 <= 2D: Theorem D / (VII.2) exactly. For 2h > D:

* e >= 0 only for n = D odd and h = (D+1)/2, and then e = 0. These are rigid and never dense.
* Every non-squarefree type (n < D) and every even D has e < 0. Even D at h = D/2 + 1 has e = -2, and so does the
  K5 type at D = 10.

**Combinatorics (exact, `cuttypes.py`).** Pair degrees satisfy s_ij >= h, and sum_{i<j} s_ij <= 6D (Plotkin; this is
Astra's 3/5). Exhaustive enumeration of cut multisets up to S_5:

| D, h | 3,2 | 4,3 | 5,3 | 6,4 | 7,4 | 8,5 | 9,5 | 10,6 | 11,6 | 12,7 | 13,7 | 14,8 | 15,8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| squarefree types | 0 | 0 | 0 | 0 | 2 | 0 | 1 (K5 - edge) | 1 (K5) | 18 | 0 | 24 | 8 | 174 |
| non-squarefree (chains of nested cuts) | 0 | 0 | 0 | 0 | 10 | 0 | 18 | 10 | – | 0 | – | – | – |

Every non-squarefree type splits into a squarefree cut multiset with the same pair degrees, so `chains.py`
enumerates it completely from the squarefree list. The two D = 7 types agree with round30 search.md Cor. 2.2.
The single D = 9 type agrees with its 6 oriented types. The D = 10 type is K5, the function-field form of the
record's "ten-pair profile" (README_VIL: unique LP optimum 3W/5).

## 2. Record compilation (task item 1)

| source | points | D | h | 2h vs D | moves? | dense in T? |
|---|---|---|---|---|---|---|
| ff Prop 3.1/4.6 | 6 | 4 | 2 | = | yes | no (2h = D) |
| ff §3 (exact search) | 5 | 6 | 3 | = | yes | no |
| cube, ff Thm C | 8 | 8 (and 8k by composition) | 4 (4k) | = | 2-parameter (cube locus) | no: 56 five-subsets have alpha = 4, sweep hypersurfaces H_S and subtori S_ab|cd (bridge §4.2) |
| search Thm 5.1 (certified) | 8 | 12 | 6 | = | not computed | no |
| search §4 (certified) | 5 | 10 | 5 | = | expected dim 2 | no (swept dimension <= 3) |
| record M(7,1) >= 8, M(5,1) = 10 | >= 5 | 7, 5 | 3, 2 | < | yes | (5 points, D = 7, h = 3: e = 4) |
| search: M(7,0) >= 4, M(9,0) >= 4 | 4 | 7, 9 | 4, 5 | > | yes (four points) | four-point question |
| Astra VII.5 five-mark Pell | 5 | 8n+1 | 4n | < (h = (D-1)/2) | dense | yes; checked n = 1, 2 independently (`pell_check.py`) |
| Astra four-mark Pell | 4 | 6n+1 | 4n | > | dense | four marks only |
| Astra three marks | 3 | 3d | 2d | > | dense | three marks only |
| **this report** | 5 | 11, 13 | 6, 7 | **>** | **rigid (isolated)** | **no** |

Before this report the record had no five-point configuration with 2h > D in characteristic 0, real or split.
The record's integer five-point C -> 0 source (the cyclotomic family) is not a rational branch (bridge §5).

## 3. New computations (task item 2)

### 3.1 D = 7, h = 4: empty (exact)
`gb.py` uses Singular over **Q**. Variables rho_k, s_k and u with u·prod(rho_k - s_k) = 1. Gauge:
sum e_k(rho_k + s_k) = 0 and rho_1 - s_1 = 1. Then iterated saturation by all rho_k - rho_l, s_k - s_l, rho_k - s_l.
Both squarefree types and all 10 non-squarefree types give the ideal (1). Positive control (`ctrl4.py`): the
four-point D = 7, h = 4 types with 3 random slices give zero-dimensional saturated ideals of degree 58 or 64. So the
saturation does not kill genuine loci. The coefficient formulation (`coef.py`, cut groups G_tau, H_tau with
resultant saturation) agrees. **H_5(7) = 3** (attained: M(7,1) >= 8). With Prop 3.3 of ff, **M(7,0) = 4**.

### 3.2 D = 9, h = 5
* Non-squarefree: all 18 types empty (exact mod 32003, `nsfrun.py`).
* Squarefree (K5 minus edge {0,1}): complex Newton (`newton.py`) finds solutions. 218 distinct ones in 12000
  starts, and the count was still growing. Refined to 60 digits (`refine.py`): residual <= 1e-56, separation up to
  0.09, contact exactly 5. **No found solution is real** (min |s + conj rho| = 0.48). Real Newton (`realnewton.py`),
  3000 starts: none. The exact count did not finish within the session in either formulation (the coefficient-formulation Singular run was stopped unfinished) (Singular mod p
  saturation, or msolve with per-pair Rabinowitsch variables, which ran out of memory). So "no real solution at
  D = 9" is **evidence**. Real-Newton calibration: D = 8, h = 4 types give real solutions readily.

### 3.3 D = 11, h = 6: real solutions (verified)
Real Newton with 1000 starts per type finds real nondegenerate solutions in types 1, 6, 9 and 12 of 18. All are
refined to 60–100 digits with residual <= 1e-59, contact exactly 6, realness error 0 (`realref.py`). Jacobian
smallest singular value: 2.4e-4 (type 1), 2e-3 (type 9), 9.5e-3 (type 12). Several sit near a degeneration (one root
close to the conjugate of another, relative separation 3e-4); type 12 does not (separation 0.079).

**Explicit example (type 12).** The columns give, for each conjugate pair rho_k / conj(rho_k), which marks take rho_k:

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| marks with rho_k | {1} | {1,2} | {1,3} | {1,4} | {2,3} | {2,3} | {2,4} | {3,4} | {1,2,4} | {1,3,4} | {2,3,4} |

The 40-digit roots are in `d11_type12_roots.txt`, starting rho_0 = 0.73209780609671940197...+1.0i. The standalone
check `verify_d11.py` builds z_i = prod_{k: i in S_k}(t - rho_k) prod_{k: i not in S_k}(t - conj rho_k) for
i = 0..4. Every pair difference has coefficients of t^10..t^6 at most 1.3e-39 (the input precision), and its
t^5 coefficient is of order 0.1–1.5. So **D = 11, h = 6 exactly**, the norms agree to 2e-79 and min |Im rho| = 0.13.
The tuple is primitive (no column is constant), N is squarefree, and the points are distinct (rank 4). The numbers
are algebraic, being isolated roots of a system over Q. A PSLQ search found no polynomial of degree <= 30 with
coefficients <= 1e12 for an affine invariant, so there is no sign that the tuple lives over Q(i).

### 3.4 D = 13, h = 7 (verified); D = 10, 12, 14 (negative evidence)
* D = 13: real solutions in 16 of the 24 types (500 starts each). 35 refined roots, all with contact exactly 7 and
  realness error <= 1e-89. Example: type 2, separation 0.031, Jacobian sigma_min = 0.046.
* D = 10, h = 6 (K5, e = -2): 1500 complex starts, no nondegenerate solution.
* D = 12, h = 7: no combinatorial type exists, so H_5(12) = 6 (attained by the half-cube).
* D = 14, h = 8 (8 types, e = -2): 300 complex starts each, none.

### 3.5 D = 15, h = 8
Real Newton on 10 random types out of 174, 300 starts each: no real nondegenerate solution. About 90% of starts
did not converge, against about 75% at D = 13, so this is an underpowered, **inconclusive** negative. It is not
evidence against D = 15. (This D gives delta = 1/30 exactly.)

### 3.6 Table of H_5(D) in characteristic 0

| D | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H_5 (conjugation model) | 1 | 2 | 2 | 3 | **3 (exact)** | 4 | 4 (evid.) | 5 (6 excluded: evid.) | **6 (verified; 7 impossible)** | 6 (exact upper) | **7 (verified)** | 7 (8 excluded: evid.) |
| H_5 (split model) | same | | | | 3 (exact) | | **5 (verified)** | | >= 6 | | >= 7 | |

Upper bounds for D <= 6, D = 8 and D = 12, and "7 impossible" at D = 11, come from the Plotkin/type enumeration.
Lower bounds at h = floor(D/2) come from the record (M(D,0), M(D,1) tables, cube, half-cube) and from Astra's Pell
at D = 9.

## 4. Standing cautions
* *Primitive height.* Every configuration is checked primitive: non-constant columns, all 2n roots distinct, no
  gcd. The h = 6 at D = 11 is the primitive contact.
* *Uniform exceptional set.* Rigid curves give no uniform statement either way. Data at finitely many D cannot test
  density of a countable union of curves, and I claim nothing about it.
* *Template-dependent thresholds.* The ratio (D+1)/(2D) is D-dependent. These examples give, for each delta, only
  configurations with D <= 1/(2 delta), plus their compositions.

## 5. Novelty and significance
**Novel:**
* the first characteristic-0 five-point configurations with 2h > D, real model, at D = 11 and D = 13;
* the exact closure of D = 7 (H_5(7) = 3, M(7,0) = 4);
* the non-real split solutions at D = 9;
* the exhaustive type tables.

**Significance for the goal: modest.** They show that Astra's strict inequality h > D/2 is not vacuous: genuine
tuples exist above it. Any proof of the FF-A_5 statement must therefore allow an exceptional set containing rigid
curves with alpha < 15/4, as the cube already required at alpha = 4. They give no density and no growth-rate
improvement. They do not transfer to integers, since they are not over Q(i) (Astra's VII.7 specialisation lemma
needs Q(i) and density).

**Next step.** Certify the D = 11 type-12 root with interval Krawczyk. Then decide whether infinitely many odd D
carry real rigid solutions that are not compositions, and whether their ratio curves could accumulate to a dense
set. The first test: are the D = 11 and D = 13 curves contained in a common proper subvariety (e.g. the hypersurfaces
H_S of bridge §4.2)?
