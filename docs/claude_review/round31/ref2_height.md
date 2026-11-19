# Referee report, round 31, block "height": astra_lemma_II/answer.md §§II.1, II.2, II.3

Scope: lines 246-316 of the frozen `astra_lemma_II/answer.md`:
* II.1: the height cutoff (II.1), the bounded-gcd example, (II.2), and the corrected grid conclusion;
* II.2: the index-one count at the growth scale, and the two qualifications, including squarefree index in {1,2,3};
* II.3: boundary-branched covers, non-thinness, the cubic cover and (II.4), and the uniqueness calculation.

Part I §2, the theorem II.1 relies on, was refereed CORRECT in `ref_cover.md` and `ref_baker.md`, and I take it as
given.

Scripts are in `round31/ref2_height/`. They are my own code, and all runs used `prlimit --as=3.5e9`.
* `j_identity.py`: exact Gaussian-rational check of (II.4) and of the beta/beta-bar exponent formula.
* `modk.py`: Smith normal form of the formal k-th-power system.
* `binary_index.py`: exhaustive enumeration of 2^16 binary matrices and of 2^15 squarefree column sets.
* `gcd_example.py`: the u_j = t + c_j + i example.
* `conc.py`: a sanity check of (II.2) on Sylvester profiles.

Each script has a matching `*_out.txt`. I also read Astra's `check_part_ii.py`, and its assertions are consistent
with mine.

## Verdict table

| # | claim | verdict | novelty / significance |
|---|---|---|---|
| 1 | Cutoff (II.1): H_S <= C^{-g}N^{g/4}; exponent 47/50 at index >= 7; all H_S <= N at index >= 9; squareclass forces abs(G_S)^2 <= 4C^4 | CORRECT | calibration of §2 against J4's scheme; peripheral |
| 2 | Bounded-gcd squareclass example u_j = t + c_j + i | CORRECT (gcd norm is 2 or 8, not 1; harmless) | illustrative |
| 3 | (II.2): all but O_{C,eps}(m^4) ordered 5-tuples have log H_S / W within eps of 15/16 | CORRECT | mean 15/16 is in the record (round23 split.md, J4); concentration via the Gram bound is new |
| 4 | Corrected grid conclusion (Schwartz-Zippel, D m^4) | CORRECT | the record's J4 argument specialised to §2; peripheral |
| 5 | Growth scale W <= A m log m: all but O_{A,C}(m^4) tuples have index exactly 1, all unit classes | CORRECT (one harmless constant slip) | superseded by the parallel thread's O(m^3), ref2_index #5/#8 |
| 6 | Without growth: O_C(m^4) index failures for genuine clusters; every full-rank squarefree index is in {1,2,3} | CORRECT, reproduced exhaustively | elementary; not in the record |
| 7 | Covers étale over T are translated isogenies, and lifting needs columns in the image lattice (excludes I = 1) | CORRECT | standard; it corrects the coordinator's ramification-index argument |
| 8 | Ten-prime squarefree ratios are not thin (Cohen bound) | CORRECT; the cited bound exists in a stronger form (secondary sources only) | peripheral: ignores the arc condition |
| 9 | Cubic cover v^3 = J/b, identity (II.4), multiplicity 3, MR ratio unchanged | CORRECT, (II.4) verified exactly | abs(II.4) follows from the record's chord dictionary; no gain |
| 10 | Uniqueness: the cubic is the only formal block-power pattern | CORRECT, verified by Smith normal form | formal only |

No FALSE item and no GAP in this block.

## 1. The cutoff (II.1)

**Sufficient condition.** Dividing by G_S leaves the ratios unchanged. At each p it replaces e_p by the range of the
a_ip over S. Indeed, the gcd's pi-valuation is min a and its pi-bar-valuation is e_p - max a, so the exponent of p in
H_S is max - min. The diameter ratio Delta is scale-invariant, and chord <= angle gives Delta <= C N^{-1/4}.
Part I §2 gives Delta > H_S^{-1/g} off Z. So a contradiction (membership in Z) follows when
C N^{-1/4} <= H_S^{-1/g}, that is, H_S <= C^{-g} N^{g/4}. At g = 94/25, g/4 = 47/50. **Checked.**

**Lower bound.** Using d_ij <= log H_S, (1) gives H_S >= C^{-2} N^{1/2}, so H_S >= H_0 eventually. **Checked.**

**Index 9.** 70*9 = 630 > 625, so B_9 > 4. Then g/4 > 1 and every H_S <= N qualifies once N is large. **Checked.**

**Squareclass case.** (4) is 2^{-1/2} H_S^{-1/4} <= Delta <= C N^{-1/4}. This gives abs(G_S)^2 = N/H_S <= 4C^4.
**Checked.**

Astra's point is right and matters for the coordinator's question 1: "within a bounded factor of N" is the wrong
regime, because §2 at g = 94/25 needs a power saving N^{3/50}. The record's J4 argument (round23 `split.md` (iii)-(iv),
expert_note l.392) has exactly this shape: it uses (15/16 + eta)(1 + 2 delta) < 1.

## 2. Bounded-gcd example (`gcd_example.py`)

I checked each of the six choices of five 2-subsets, for t <= 3000 and for t = 10^6+17, 10^9+7, 10^15+3:
* the norms are equal;
* the five points are distinct;
* (Delta N^{1/4})^2 <= 324, so C <= 18 is absolute;
* two complementary pairs are always present, giving z_A z_{A^c} = N = z_B z_{B^c}, an exact multiplicative relation.

The polynomial gcd over Q(i) is 1. The Gaussian gcd at integer t has norm exactly 2 or 8, because (t+c)^2 + 1 is even
for odd t+c. After division the norm is odd (checked, 0 failures), so the tuple is primitive. The quadratic
coefficients have Im = -(2 sum_A c - 11), and the pair sums 1,3,7,4,8,10 are distinct.

**CORRECT.** The only point to note is that the gcd is bounded but not trivial, as Astra's wording allows.

## 3. (II.2)

I checked the proof line by line.

**Layers.** log H_S/W = 1 - eta, where eta is the weight of layers constant on S. The identity
1_{all equal} = 2^{-4} sum_{abs(J) even} prod sigma holds.

**Second moment.** For abs(J) = 2 or 4, E f_J^2 = sum_{l,l'} (w w'/W^2) x_{ll'}^{abs(J)}, with
x_{ll'} = m^{-1} sum_i sigma_il sigma_il'. Since abs(x) <= 1, this is <= norm(V^tV)_F^2/m^2 <= norm(V)^2 tr(VV^t)/m^2,
and tr(VV^t) = m.

**Gram bound.** Inside a piece, <V_i,V_j> = 1 - 2d_ij/W < 0. To see this, write psi for the half phase difference.
Then abs(y_i/y_j - 1) = 2 abs(sin psi) >= e^{-d/2}. The piece width bounds abs(psi) by N^{-1/4}/4, so d_ij > W/2.
The Gram matrix A is PSD, so Id - A is entrywise >= 0 with eigenvalues <= 1. Perron-Frobenius gives spectral radius
<= 1, hence norm(A) <= 2. Over the pieces, V^tV = sum_k V_k^t V_k, so norm(V)^2 <= 2L.

**Chebyshev.** Applied to the 15 nonempty even J, it gives an exceptional fraction O(L/(m eps^2)), that is,
O_{C,eps}(m^4) ordered tuples. **CORRECT.**

**Sanity check (`conc.py`).** On Sylvester-Hadamard profiles the empirical exceptional probability halves as m
doubles (0.22, 0.11, 0.057 at m = 64, 128, 256), and the mean is 0.9374.

**Novelty.** The mean 15/16 is in the record (round23 `split.md` (ii), with a positive-proportion conclusion). The
concentration statement with O(m^4) exceptions, via norm(V)^2 <= 2L, is new and slightly stronger.

## 4. Corrected grid conclusion

The margin is 47/50 - 15/16 = 1/400. Z is a proper closed set over Q, and in split coordinates it is cut out by a
Laurent polynomial in the y_j. Homogenising with z_5^D gives a nonzero form of degree D. Schwartz-Zippel on X^5 then
gives <= D m^4 tuples, and the gcd-atypical tuples add O_C(m^4). The squareclass count follows because
abs(G_S)^2 <= 4C^4 forces log H_S >= W - O(1), which is atypical. **CORRECT.**

This is the record's J4/Thm 1 argument (split.md (v)) applied to the index >= 7 theorem. It is peripheral.

## 5. The growth-scale index-one count

**First bound.** I checked e_p log p < 2W/s <= 4LW/m. The threshold count gives W_B > W(s-1)/s, and the half-interval
choice gives W_B <= W - e_p log p/2 (round8 `stability.md` Thm 3.3(2) has the same bound). Then:
* the small primes contribute <= 4LW/sqrt(m);
* e_p <= 8LA at p > sqrt(m);
* the discards (log H_S > 31W/32, and abs(<V_i,V_j>) > 1/128 using tr((VV^t)^2) <= 2Lm) cost O(m^4);
* typical pairs have d_ij >= 127W/256 > 31W/64;
* the barycentre argument with E abs(X-Y) <= 1/2 gives real affine independence on the large primes;
* the bounded minor is 24 ceil(8LA)^4.
All **checked**.

**Mod-l count.** The step "at most l rows per residue class per piece" is correct. One constant needs correcting: for
genuine points the gap is arcsin(e^{-w(v)/2}/2), not arcsin(e^{-w(v)/2}), because abs(u-1) = 2 abs(sin). The
inequality N^{-1/4}/(4l) < e^{-W/(2l)}/2 still holds for l >= 2. The first affinely dependent row lies in an affine
span of <= l^3 residues, giving O(L l^4 m^4). **CORRECT.**

**Significance.** It is superseded: the parallel thread gives O(m^3) with no growth hypothesis (ref2_index items 3, 5
and 8, refereed CORRECT).

**Genuine clusters without growth.** Index in [2,6] has a prime factor 2, 3 or 5, and index >= 7 is handled by item 4,
so the count is O_C(m^4). **CORRECT.**

**Squarefree index (`binary_index.py`).**
* The 2^16 binary 4x4 matrices have determinants {0: 42976, ±1: 10020 each, ±2: 1200 each, ±3: 60 each}, so the
  maximum is 3.
* The bordering identity abs(det5) = 16 abs(det4) holds on all 2^16 matrices.
* There are 15 classes of nonconstant 0/1 allocation vectors modulo complement. For all 2^15 subsets of them, the
  lattice index (the gcd of the 4-minors) is distributed as: rank < 4: 1306; I=1: 31227; I=2: 230; I=3: 5.
* So {1,2,3} is exact, and index 3 occurs. Astra's columns 12, 13, 14, 23 have det 2 (confirmed).

**CORRECT.** The claim is classical (the maximal 0/1 determinant of order 4 is 3). It is new to the record and
elementary.

## 6. Covers branched only at the boundary

A connected finite étale cover of G_m^4 over C is an isogeny once a base point is chosen (pi_1 = Z^4, plus Riemann
existence). With a rational point upstairs, the pointed group law descends by uniqueness of lifts. The cocharacter
lattice upstairs again carries the -1 action, so T'(Q) is N' ⊗ U(Q). Rational liftability of y t^{-1} then needs every
prime column outside supp(t) to lie in f(N'). With I = 1 and a unit (torsion) translate, this is impossible in degree
> 1. **CORRECT** and standard.

Astra's caveat is also right. The coordinator's argument, "ramification index e divides contact multiplicities", is
valid only on the sheet met by the section at a good model. The isogeny argument is the correct proof of the
conclusion at I = 1 for boundary-branched covers.

## 7. Non-thinness and Cohen's bound

**The pullback.**
* B -> B/B-bar has kernel G_m.
* The 4x10 pair-column matrix has trivial Smith invariants, since its column lattice is Z^4 (Part I §4, re-derived).
* So Phi is a surjective torus homomorphism with geometrically connected kernel, and Q(A^20)/Q(T) is regular.
* For any type-2 thin datum X -> T of degree >= 2, Q(X) ⊗_{Q(T)} Q(A^20) is a field. The pullback is therefore
  irreducible of the same degree, and preimages of thin sets are thin. A point of A^20(Q) over pi(X(Q)) lifts to the
  fibre product.

**The count.** Gaussian generators of primes 1 mod 4 with norm in [X^2/2, X^2] have abs(Re), abs(Im) <= X. Ten ordered
distinct choices give >> X^20/(log X)^10 integral points, against O(X^{39/2} log X). The tuples z_i = P_i prod_{e not
containing i} B-bar_e are primitive and distinct. **CORRECT.**

**Source check.** The statement that exists is stronger than the one Astra uses: for a thin Omega ⊂ A^n(Q), the
number of integral points with max abs(x_i) <= B is O(B^{n-1/2} (log B)^gamma), with gamma < 1 and the constants
depending on Omega. This is Cohen, PLMS 43 (1981), as presented in Serre, *Lectures on the Mordell-Weil Theorem*
§13.1 Thm 1, and in *Topics in Galois Theory* §3.4. I confirmed this **only through secondary sources** (search
snippets of Serre's book and of a 2026 survey, arXiv 2603.23334). The fetches of arxiv.org, Wikipedia and a PDF host
failed, because those hosts are not allowed in this session. I did not read Cohen's paper. The record never cites
Cohen.

**Significance.** Peripheral. The short-arc points have density zero in this parameter family, so the theorem says
nothing about whether they lift. Astra says this herself ("No angle constraint").

## 8. The cubic cover and (II.4) (`j_identity.py`)

**Identity.** I verified (II.4), J(u) = prod_{i<j} t_ij / H^3 with signed t_ij = det(P_i,P_j)/n_ij, exactly in
Q(i). The test used 300 random sets of ten Gaussian blocks with entries in [-60,60]. No primality or coprimality is
needed, and t_ij is integral each time.

**Hand derivation.** u_i - u_j = -2i t_ij n_ij/(P-bar_i P-bar_j), and prod abs(P_i)^4 = H^4.

**Rationality on the twist.** Inversion invariance holds: (-1)^10, and each index appears four times.

**Lift and geometry.**
* The lift is v = r/H when prod t = b r^3.
* The multiplicity at the identity is min(3,10) = 3.
* Normality follows from the hypersurface criterion: R1 because the singular locus is v = 0 over the codimension-2
  singular locus of the diagonal arrangement, and S2 because it is a hypersurface.
* abs(v) = O(Delta^{10/3}).
* McKinnon-Roth Cor 4.2 (n/(n+1)(Vol/mult)^{1/n}, as checked in `papers_check/check.md` item 4) gives (3*70)/3 = 70.
  The upstairs distance is still ≍ Delta, through the u-coordinates.

**CORRECT.** It gives no gain, and Astra says so.

**Novelty.** abs(II.4) is just the product over the ten pairs of the record's chord dictionary
abs(u_i - u_j) = 2 abs(t_ij)/sqrt(D_ij), using prod_{ij} D_ij = H^6. Only the signed, rational form, which is what
the cube-class lift needs, is new.

## 9. Uniqueness calculation (`modk.py`)

**Exponent formula.** I verified the formula exactly for random integer exponents with the homogeneity constraint
(300 cases): beta exponent v_ij + a_i + a_j, and beta-bar exponent v_ij - s_i - s_j - a_i - a_j, with the prefactor
(-2i)^{sum v} prod t^v.

**Solution count.** The 21x15 system (20 block conditions plus homogeneity) has Smith invariants 1^14 and 3. So the
number of solutions mod k is gcd(3,k) for k = 2..24 (checked), and they are generated by J's exponent vector
(a = -2, v = 1). That vector maps to (-3, ..., -3, 0).

**CORRECT.** As Astra states, it is a formal test only.

## Bottom line

The block is clean: nothing false, and no gap.

It is also peripheral. Item 1 is a calibration of an index >= 7 theorem that acts on a set the parallel thread
already shows has O(m^3) tuples. Item 5 is superseded. Item 6 is elementary. Item 8 ignores the arc. Item 9 gives a
lift without an exponent gain. Nothing changes the growth rate (2/3 + eps) log R/loglog R or the uniform bound, and the
status of Questions A and B and route 5 is unchanged.

The one useful clarification for the coordinator is the negative answer to question 1 as posed: the fixed-Z inclusion
needs H_S <= N^{47/50}, not H_S ≍ N. The only covers that lift index-one squarefree points are branched in the
interior, and the one exhibited (the cubic) has a multiplicity-3 target with no volume gain.
