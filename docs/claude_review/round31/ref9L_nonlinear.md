# Referee report, block NONLINEAR: Part IX.5, the PTE regression and the "comparison of scales" in IX.6

Scripts, all my own and none importing Astra's `ix_checks.py`, are in
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref9L_nonlinear/`:
`walsh.py` (coefficient bounds), `divis.py` (twin divisibility and the Gaussian identity), `quad.py` (quadratic
profile), `pte.py` (PTE regression). All arithmetic is exact (Python integers; Gaussian integers as integer pairs).

## 0. Verdict

**CORRECT.** Every step of IX.5 checks out: the counting lemma, the twin proof for any d, the partial-coverage
version and the circle proof. The regression numbers of IX.6 are reproduced exactly. The "comparison of scales"
paragraph agrees with `ref8L_attack.md`. I found no false statement. Below are one remark on scope, which should
appear next to the theorem, and some minor sharpenings.

* **Scope (should be stated).** If a profile is in the class (d, n), every bit is an F_2-polynomial of degree at
  most d in r = log_2 n variables. So the signature matrix has F_2-rank at most R_d(r) = sum_{j<=d} C(r,j) =
  O((log n)^d). Consequently, if n = O(s), the theorem only concerns clusters whose binary signature rank is
  polylogarithmic in s. All the extremal fakes in the record (Paley, skew-Hadamard, Legendre P2) and the IX.4
  relaxation have binary rank of order m. For them the class requires n >= exp(c m^{1/d}), and the bound is vacuous.
  This is the reconciliation with the fakes. Astra's degree argument for P2 in IX.6 is a special case of it.
* **Quantitative form of "n = O(s)".** The bound is o(s) exactly when n = o(s^{8/7}) for d = 2 (in both problems).
  For general d the twin needs n = o(s^{1/(1-2^{1-2d})}) and the circle needs n = o(s^{1/(1-2^{1-2b_d})}).
* **The p = 2 worry in the brief does not arise.** The divisibility q/Q_u | F_u never divides by 2^t. It uses the
  unnormalised integer Walsh sum of the bit, and that sum is exactly 0 whenever the coefficient vanishes.

## 1. The counting step (Proved)

**Statement.** Let V = F_2^r, n = 2^r. Let f_t(n) be the maximum size of a subset of V containing no affine
t-dimensional subspace. Then f_1(n) = 1 and B(B-1) <= 2(n-1) f_{t-1}(n/2). Hence f_t(n) <= c_t n^{1-2^{1-t}} with
an effective constant c_t (numerically c_2, c_3 -> sqrt 2 and c_4 -> 1.30 when the recursion is iterated with
equality).

**Check.** An affine 1-space is a pair {x, x+h}, so f_1 = 1. Fix h != 0. The set A ∩ (A+h) is a disjoint union of
pairs {x, x+h}. Let E_h be its image in V/<h>, which is a copy of F_2^{r-1}. Then |A ∩ (A+h)| = 2|E_h|. The full
preimage of an affine (t-1)-space of V/<h> is an affine t-space containing the direction h. So if E_h contained
an affine (t-1)-space, A would contain an affine t-space. Hence |E_h| <= f_{t-1}(n/2). Summing over the n-1 shifts
h != 0 counts each ordered pair of distinct elements once:

  B(B-1) = sum_{h != 0} |A ∩ (A+h)| <= 2(n-1) f_{t-1}(n/2).

Induction gives B^2 << n · n^{1-2^{2-t}}. The exponents are 1/2, 3/4, 7/8 for t = 2, 3, 4. The case t = 2 is the
Sidon bound B <= 1 + sqrt(2n). The lemma is used on labels, so the nonlinearity of the signature map plays no part.
In the applications, the pieces of the interval (or arc) and the residue classes of the winding numbers are fixed
before the flat is chosen. The lemma is applied to the set of labels inside one piece, and the factor for the
number of pieces depends only on C and d.

**Literature (search snippets only; I read no full text).** The search found Bonin–Qin, *Size functions of
subgeometry-closed classes of representable combinatorial geometries*, Discrete Math. 2000. A citing paper (Fox–Pham)
restates a Bonin–Qin theorem for F_2: every subset of F_2^r of size (2^{1-c 2^{-t}})^r contains an affine t-space.
The extracted text was garbled, so I could not confirm the constant c. That snippet has the same shape
n^{1-c/2^t} as the iterated Sidon bound here (c = 2). For t = 2, Sidon sets of size about 2^{r/2} exist, so the
exponent 1/2 is sharp. For t >= 3 the deletion method gives only n^{1-t/(2^t-1)} from below. This bound is not
needed. The induction above is self-contained.

## 2. Twin proof, bits depending on at most d affine forms (Proved)

**Statement (precise quantifiers).** Fix C > 0 and d >= 1. Let q >= 1, V = F_2^r, n = 2^r. For every exact prime
power p^f || q, prescribe two classes c_0(p^f), c_1(p^f) mod p^f and a bit β_{p^f}: V -> F_2. Each bit must be a
Boolean function of at most d affine forms on V; the forms and the function may depend on p^f. Call a set X ⊂ V
with distinct integers y_x (x ∈ X) admissible if y_x ≡ c_{β_{p^f}(x)}(p^f) mod p^f for every p^f and every x, and
all y_x lie in an interval of length C sqrt q. Then |X| <= (C 2^{2d-1} + 1) c_{2d} n^{1-2^{1-2d}}. The bound is
uniform in q, in the classes, in the forms and in the Boolean functions. For d = 2 it is B <<_C n^{7/8}.

The bound is uniform over all prescriptions, so it also holds a posteriori. If a configuration admits some
injective labelling into F_2^r for which every bit has this form, the bound holds.

**Check of each step.**

1. *Pieces.* Divide the interval into ceil(C 2^{t-1}) + 1 pieces of diameter < 2^{1-t} sqrt q, where t = 2d.
2. *A t-flat.* Suppose that the labels of one piece contain an affine t-flat. Identify it with F_2^t and put
   F_u = sum_x χ_u(x) y_x. For u != 0 we have sum_x χ_u(x) = 0. So F_u = sum_+ (y_x - o) - sum_- (y_x - o), where o
   is the minimum mark of the flat, and |F_u| <= 2^{t-1} · diam < sqrt q.
3. *Spanning.* Suppose the frequencies u != 0 with F_u != 0 lie in a hyperplane h^⊥ with h != 0. Then
   x -> y_{x+h} has the same Walsh transform as x -> y_x, because χ_u(h) = 1 on the support and F_0 is unchanged.
   Inversion over Q gives y_{x+h} = y_x. This contradicts distinctness. So a basis U of t frequencies with
   F_u != 0 exists.
4. *Fourier support.* Restricted to the flat, each affine form is affine. So each bit factors through a linear map
   F_2^t -> F_2^{d'} with d' <= d, and it is constant on the cosets of a subspace K of codimension d'. Its
   coefficient at u vanishes unless u ∈ K^⊥, and dim K^⊥ = d' <= d. A linearly independent set meets a
   d'-dimensional subspace in at most d' elements. So each p^f lies in Q_u for at most d of the u ∈ U, and
   prod_{u∈U} Q_u divides q^d.
5. *Divisibility.* Modulo p^f, y_x ≡ c_0 + (c_1 - c_0) β(x). So F_u ≡ (c_1 - c_0) S_β(u) mod p^f, where
   S_β(u) = sum_x χ_u(x) β(x) is an **integer**, and it is 0 exactly when the coefficient at u vanishes. No division
   by 2^t occurs, so the case p = 2 is identical to p odd. The prime powers are pairwise coprime, so
   q/Q_u | F_u. Exact test (`divis.py`): q = 8·9·25·(13 primes up to 53), r = 6, d = 2, random Boolean functions of
   2 random affine forms, marks by CRT, 300 random 4-flats, 4500 frequencies. All divisibilities held, and every
   bit's Fourier support on the flat had rank <= 2.
6. *Contradiction.* From 0 < |F_u| < sqrt q we get Q_u > sqrt q for each of the 2d elements of U, so
   prod Q_u > q^d. This contradicts step 4. Hence no piece contains an affine 2d-flat of labels, and §1 applies.

The case d = 1 gives n^{1/2}, consistent with IX.3 for linear signatures. Remarks: q may be arbitrary (small q
gives nothing new); prime powers are allowed; two classes that agree modulo p^{f-1} are allowed.

## 3. Partial coverage (Proved)

**Statement.** Fix C, d and 0 < ε <= 1/2. Suppose that the bounded-test condition holds on a set G of exact prime
powers with product q_g >= q^{1/2+ε}. The classes at the other prime powers are arbitrary. Then
B <<_{C,d,ε} n^{1-2^{1-t}} with t = ceil(d(1/2+ε)/ε).

**Check.** Form Q_u from G only. Step 5 gives q_g/Q_u | F_u, and |F_u| < sqrt q gives Q_u > q_g/sqrt q. We need
t(log q_g - ½ log q) >= d log q_g. Write x = log q_g/log q ∈ [1/2+ε, 1]. The required t is d x/(x - 1/2), which
decreases in x, so its maximum d(1/2+ε)/ε is attained at x = 1/2+ε. The chain of inequalities is strict
(prod Q_u > (q_g/sqrt q)^t >= q_g^d >= prod Q_u). At ε = 1/2 this recovers t = 2d. The pieces have diameter
< 2^{1-t} sqrt q, so the constants depend on C, d and ε. The statement "only a strict majority of the logarithmic
weight is required" is correct.

## 4. Circle proof (Proved)

**Setting.** Divide the cluster by its Gaussian gcd. The normalised norm N' contains only split primes, because
the parts at 2 and at primes 3 mod 4 are common up to units. The angular diameter is C N^{-1/4} <= C N'^{-1/4}.
Fix π_p, retain one unit η (a loss of at most a factor 4), and write z_x = η prod_p π_p^{a_xp} conj(π_p)^{e_p - a_xp}.
Then arg z_x = arg η - sum e_p φ_p + 2 sum_p a_xp φ_p.

* **The factor 2 is correct.**
* **β = arg η - sum e_p φ_p is common to the retained marks and is kept.** It cancels only in the combinations
  F_u with sum_x χ_u(x) = 0, that is, for u != 0. No drifting phase is set to zero.

Normalisation replaces a_xp by a_xp - min_x a_xp. The threshold bits of the normalised exponents are threshold bits
of the original ones, so the class (d, n) is preserved under normalisation and under passing to sub-clusters. The
proof uses one fixed factorisation of N' for all marks, so no gcd cost arises.

**Steps.** Fix t = 2b_d and m_0 = 2^{t-d}. Make pieces of angular diameter <= 2^{-d} N'^{-1/4} (about C 2^d + 1 of
them) and retain one class of k_x mod m_0. All these choices precede the choice of flat. On a t-flat,
F_u = sum χ_u θ_x = 2 sum_p φ_p sum_x χ_u(x) a_xp + 2π sum χ_u(k_x - k*). The first sum is m_0 b_{p,u}. The
second is divisible by m_0.

1. **Integrality of b_{p,j,u}.** The bit 1_{a_xp >= j} is constant on the cosets of K, with |K| = 2^{t-d'} and
   d' <= d. Each coset sum of χ_u is 0 or ±|K|. So m_0 = 2^{t-d} divides sum_x χ_u(x) 1_{a_xp >= j}.
   (Computed: `walsh.py`, all 312 such functions on F_2^4 with d = 2 and all 72 on F_2^3: no non-integral b.)
2. **(IX.5).** b_{p,j,u} = 2^{d-d'} S(u'), where S is the unnormalised Walsh sum of a Boolean function on
   F_2^{d'}. Parseval gives sum_{u'≠0} S^2 = k(2^{d'} - k) <= 4^{d'-1}, where k = |support|. So the sum of b^2 over
   u ∈ U ∩ K^⊥ is at most 4^{d-1}, and at most d' <= d such u exist. Cauchy–Schwarz and integrality give
   sum_{u∈U} |b_{p,j,u}| <= floor(2^{d-1} sqrt d) = b_d. Then sum_u log H_u <= sum_p log p sum_j sum_u |b_{p,j,u}|
   <= b_d W.
   * Computed exact maxima: d = 2: max 2 = b_2 (both over all 840 bases of F_2^4 and by the quotient formula);
     d = 3: max 6 = b_3; d = 4: max 12 < b_4 = 16.
   * So b_2 = 2 is sharp and b_d is not sharp at d = 4. This is harmless.
   * Astra's explanation is right: AND-type tables have three coefficients ±1, and their frequencies lie in a
     plane, so only two of them can be in a basis. Her table 0, 2, 3 with frequencies 2, 6, 8 is reproduced.
3. **(IX.6).** exp(i F_u/m_0) = prod_p (π_p/conj π_p)^{b_{p,u}} = G_u/conj G_u, where G_u takes π_p^{b} for b > 0
   and conj(π_p)^{|b|} for b < 0, and |G_u|^2 = H_u. Computed (`divis.py`): with 10 split primes, caps e_p <= 3,
   d = 2, t = 4, m_0 = 4 and 1500 flat/frequency pairs, the exact identity prod_x z_x^{χ_u(x)} = (G_u/conj G_u)^{m_0}
   held in Z[i] in every case.
4. **The gap.** |G/conj G - 1| = |G - conj G|/|G| = 2|Im G|/sqrt(H_u). Why Im G_u != 0: if Im G = 0 then
   G/conj G = 1, so F_u/m_0 ∈ 2πZ. But 0 < |F_u/m_0| <= 2^{t-1} 2^{-d} N'^{-1/4}/2^{t-d} = N'^{-1/4}/2 < 2π, a
   contradiction. (Alternatively, by unique factorisation G is never real when b_u != 0, because each prime
   occurs on one side only.) Im G_u is a rational integer, so 2 H_u^{-1/2} <= |e^{iF_u/m_0} - 1| <= N'^{-1/4}/2.
   Hence H_u >= 16 sqrt N'.
5. **Contradiction.** t elements give sum log H_u > tW/2 = b_d W, against (IX.5).

**Result.** M <= 4 (C 2^d + 1) 2^{2b_d - d} c_{2b_d} n^{1-2^{1-2b_d}}, uniform in N' and in the prescription. For
d = 2 (t = 4) this is M <<_C n^{7/8}. For d = 1 it is n^{1/2}.

**Angle models (Proved).** The proof uses three things:
* the linear system θ_x = β + 2<a_x, φ> + 2πk_x, which is part of the record's definition of an angle model (after
  retaining one class of the (π/2)Z part, at cost 4);
* distinctness of the lifts;
* one monomial gap per u, namely dist(<b_u, φ>, πZ) >= arcsin e^{-w(b_u)/2}.

The record's gap is "every monomial gap", dist(<v, φ>, (π/2)Z) >= arcsin e^{-w(v)/2} for all v, and it is
stronger. So the circle rate holds in every angle model satisfying the record's gaps, on this profile class.

**Reconciliation with the fakes (Proved).** A profile in the class (d, n) has binary signature rank
<= R_d(log_2 n). The Legendre P2 parity rows are owner singletons and have F_2-rank m on the selected marks.
Paley and skew-Hadamard profiles also have F_2-rank of order m (IX.4 gives AA^t = I over F_2 for m ≡ 3 mod 8).
So these fakes need n >= exp(c m^{1/d}), where n^{7/8} >= m, and there is no conflict. (Whether the P2 profile is
in the class is a question for another referee. The rank inequality settles it for any fixed d and any
n = exp(o(m^{1/d})).)

## 5. What class is this; examples (Proved / Heuristic)

**Definition.** Each local bit (each prime power in the twin; each threshold 1_{a_p >= j} on the circle) is
g(ℓ_1(x), ..., ℓ_d(x)), where the ℓ_i are affine forms on a label space F_2^r with n = 2^r, and the labels of the
marks are distinct. These are "linear d-juntas". They are much narrower than degree-d polynomials: a quadratic
form of rank about r is not a function of 2 forms. In particular the record's binary "two Walsh cores on x, x^3"
(trace forms of degree 2) are not in the d = 2 class.

**Examples.**
* (d = 1) Hadamard/Walsh profiles: bits u·x. Here the theorem reduces to the Sidon/sumset bound n^{1/2}.
* (d = 2) Astra's quadratic-form profile: the bits are x_i, x_i x_j (i < j) and all 2^r - 1 nonconstant parities
  u·x, so s = r + C(r,2) + 2^r - 1, which is about n. Computed (`quad.py`) for r = 3, 4, 5, 6 (n = 8, ..., 64,
  s = 13, 25, 46, 84):
  * The signatures are Sidon, with D = 1 + C(n,2) exactly. Astra's argument (h = x+y, then x_i + h_i x_j for a j
    with h_j = 1 determines {x, x+h}) is correct.
  * The real rank with the constant column is n (checked for r <= 5; for every r this follows from Walsh inversion
    over R). So the rows are affinely independent, and the sumset bound sqrt D = n and the rank bound of
    Theorem 1.2 (M <= s + O(1)) both allow the whole list.
  * The binary rank is 1 + r + C(r,2) = 7, 11, 16, 22, attaining R_2(r).

  The example is therefore caught by IX.5 and by neither earlier test, as claimed.

**Is there a natural reason for worst-case clusters to have this structure? (Heuristic)** I see none.
* The hypothesis is not implied by any known property of genuine clusters, and it is not detected by D or by
  affine rank. Its only invariant is low F_2-rank, about (log n)^d.
* The record's structural information points the other way. Every extremal fake (Paley, skew-Hadamard,
  Legendre) has F_2-rank of order m.
* The one setting where genuine clusters provably have low-rank Walsh shape is the binary model over F_2[t]
  (expert note, parity-code condition: "low-rank (Walsh) profiles, the shape of every genuine binary cluster").
  There the analogue of IX.5 fails, because Walsh inversion needs 2 to be invertible (see §6). So IX.5 shows that
  the characteristic-0 mechanism excludes large clusters of the shape that is genuine in characteristic 2.
* It gives no information about high-rank profiles, which is where an unbounded integer or circle cluster would
  have to live if it exists.

IX.5 is a correct theorem about a restricted, prescribed class. Astra does not claim it bears on the unrestricted
rate.

## 6. PTE regression and comparison of scales (Computed / Proved)

The families are recoverable from `ref8L_attack.md` §2.2. I rebuilt them independently (`pte.py`):
* the homogeneous node forms α_j T - β_j with weight α_j^2 and moments sum ε α = 0 (and sum ε β = 0 for r = 2);
* q = lcm of the node values, with every exact prime power kept, and a check that each prime power has at most
  two classes;
* the maximal-modulus families with q_2 computed from Lemma 1.1 of `ref8L_attack.md` (lcm over bipartitions of the
  within-side gcds).

Planes are 4-sets of signatures with zero F_2-sum, and 3-flats are found by doubling planes. I checked
(2L)^2 >= q and (4L)^3 >= q^2 globally, and also with L replaced by the diameter of each individual flat, which is
a stronger test.

| family (T = 1009 unless stated) | B | log q/log L | planes | 3-flats | inequalities |
|---|---|---|---|---|---|
| (4,1), weights (1,1,-1,-1) | 6 | 1.772 | 3 | 0 | hold (global and per flat) |
| (5,1) | 10 | 1.455 | 15 | 0 | hold |
| (7,2), the vectors of ref8L | 8 | 1.462 | 6 | 0 | hold |
| (8,1) | 70 | 1.137 | 7315 | 2520 | hold |
| parallelogram, maximal q_2 (a, b = 1009, 431) | 4 | 1.881 | 1 | 0 | hold |
| 5 projective marks, maximal q_2 | 5 | 1.877 | 1 | 0 | hold |
| 6 projective marks, maximal q_2 | 6 | 1.819 | 3 | 0 | hold |
| cube on 1009, 1013, 1019, maximal q_2 | 8 | 1.436 | 14 | 1 | hold |
| (4,1), (5,1), (7,2) at T = 10007 | 6, 10, 8 | 1.825, 1.542, 1.591 | 3, 15, 6 | 0 | hold |

All of Astra's counts are reproduced exactly: 3, 15, 6; 1, 1, 3; 14 and 1; 7315 and 2520. I did not recheck the
items of the same paragraph that belong to other blocks (the m = 43 relaxation counts 903 and 12,341, and the
Gaussian 2,925 and 1,818), except the Boolean-table count 0, 2, 3 with frequencies 2, 6, 8, which I reproduced.

The (k,1) families are themselves a live test of the d = 1 case. Their labels are ε ∈ {0,1}^k and their bits are
the coordinates. The per-flat inequality q <= (2 L_P)^2 is the d = 1 twin statement, and it held on every plane.

**Comparison of scales (Proved, consistent with `ref8L_attack.md`).**
* The four-, five- and six-mark constructions have limiting exponent 2. Prop. 2.1(4) of ref8L gives L ≍ sqrt q
  with a constant depending on the fixed nodes. So they are at critical scale with bounded B, as every bound
  allows.
* The eight-mark (7,2) family and the ten-mark (5,1) family have exponents 7/4 and 5/3 (ref8L table). So
  L = q^{4/7+o(1)} and q^{3/5+o(1)}, which are not C sqrt q for any fixed C. Astra's wording "node modulus grows as
  L^{7/4+o(1)}, L^{5/3+o(1)}" is right.
* At finite T the ratios are below their limits (1.462 and 1.455 at T = 1009; 1.591 and 1.542 at T = 10007), as
  ref8L §5 notes.
* Characteristic 2 (ref8L §2.4): the hyperplane family over F_{2^m}[T] has signature set equal to the first-order
  Reed–Muller (Sylvester–Hadamard) code. This is a linear code of dimension m+1, so it contains affine 3-spaces,
  and the analogue of IX.1 fails. The step that fails is the one Astra names:
  * In characteristic 2, (-1)^{u·x} = 1, so every F_u equals sum_x y_x.
  * Walsh inversion, which would derive y_{x+h} = y_x from vanishing coefficients, needs division by 2^t.
  * The size bound (degree, non-archimedean) and the divisibility bound both survive.

  The same is true of IX.5.

## 7. Labels

* **Proved:**
  * the counting lemma f_t(n) = O_t(n^{1-2^{1-t}});
  * the twin bound B <<_{C,d} n^{1-2^{1-2d}} (n^{7/8} for d = 2), uniform in all prescriptions, valid at p = 2 and
    at prime powers;
  * the partial-coverage bound with t = ceil(d(1/2+ε)/ε);
  * the circle bound M <<_{C,d} n^{1-2^{1-2b_d}} with b_d = floor(2^{d-1} sqrt d) (b_2 = 2, giving n^{7/8}),
    including the common phase, the factor 2 and the integrality of Im G_u;
  * validity in angle models satisfying the record's monomial gaps;
  * the scope remark: binary rank <= sum_{j<=d} C(log_2 n, j), so all record fakes are outside the class for
    n = exp(o(m^{1/d}));
  * the quadratic-profile example (Sidon, full real rank, caught by IX.5);
  * the comparison of scales.
* **Computed:**
  * integrality and the maxima of sum |b| (2, 6, 12 for d = 2, 3, 4);
  * 4500 twin divisibility tests and 1500 exact Gaussian identities, with no failures;
  * the quadratic profile for r <= 6;
  * all PTE plane and three-space counts and both inequalities, global and per flat.
* **Heuristic:** no natural reason for worst-case clusters to lie in the class.
* **Gap:** none in this block. The Bonin–Qin comparison rests on search snippets only.
* **False:** none.

Sources seen in search snippets only: Bonin–Qin, Discrete Math. 2000 (DOI 10.1016/s0012-365x(00)00108-4), via
https://www.bohrium.com/paper-details/size-functions-of-subgeometry-closed-classes-of-representable-combinatorial-geometries/813102332182003713-992,
and the restatement in Fox–Pham as summarised in the search results.
