# Round 31, Part VI: judge report on fixed fields, reflection and growth

Source: Part VI of `astra_free_VI/answer.md` (heading "# Part VI", lines 1258-1459; scripts `vi_profile_check.py`,
`vi_ternary_check.py`, `vi_squarefree_check.py`; coordinator's question `reply6.txt`; Astra's summary `final6.txt`).
Six referees/analysts each took one block:

| key | report | scope |
|---|---|---|
| ternary | `ref6_ternary.md` | VI.4: the F_3 construction, (VI.5), (VI.6), the certificates, legitimacy, rate, other p |
| squarefree | `ref6_squarefree.md` | VI.5: squarefree parents over F_2, the CRT trace condition, sharpness pairing with VI.2 |
| code | `ref6_code.md` | VI.1: the derivative (parity-code) condition, packing, the Paley/double-core exclusion, integer analogue |
| bounds | `ref6_bounds.md` | VI.2 (complete Hadamard cores) and VI.3 (full hyperplane designs with translation B) |
| port | `ref6_port.md` | analyst: can VI.4 be ported to Z or to characteristic 0, and which step fails |
| meaning | `ref6_meaning.md` | analyst: transfer of the record's growth bound to F_3, filters, Questions A/B, literature, wording |

The judge's scripts are in `judge_VI/`. They use the judge's own finite-field arithmetic and import nothing from
Astra or the referees:
* `ternary.py`: F_729 by log/antilog tables over a primitive sextic found by search, F_9 inside it, VI.4 with k=2.
* `gk3.py`: formal g for k=3 over the judge's own F_27 = F_3[y]/(y^3-y-1).
* `sqf.py`: VI.5 at m=4 over F_2[t]/F_4[t].
* `b8.py`: the bounds referee's 8-point cluster with a nonconstant translation.
* `pell.py`: the code referee's integer counterexample.

Outputs are in `*.out`. **Every claim labelled [J] is the judge's own and has not been refereed.**

Labels follow Parts III-V. "Proved" means a referee checked a complete proof line by line, or an exact computation
was reproduced with independent code. **No referee found an error in any claim Astra made, in the form she made
it, and no verdict on an Astra claim is GAP or FALSE.** The two FALSE verdicts in the summaries concern readings
the referees posed themselves: code (e), "VI.1 has an integer analogue", and meaning (e), "VI.1 is Mason in
disguise and abc would kill the integer fakes". I rule on both in §1.3 and redo the decisive computation.

## 0. Verdict in five lines

1. **Theorem (VI.4).** For p ≡ 3 mod 4, and in particular for **A^2+B^2 over F_9(t)/F_3(t)**, sup_D M_p(D,0) is
   infinite, with M_p(D,0) >= ((1/2)log_2 p - o(1)) log D/loglog D for all large D. That is doubly logarithmic in
   R. N is far from squarefree, and its multiplicities grow with m. The judge reproduced the four-mark instance [J].
2. **Theorem (VI.5).** Over F_2, sup_D M_2(D,0) is infinite even for **squarefree** N: m = 2^r points on a parent of
   degree 8m(m-1)(r-1). The judge reproduced m=4 [J]. Record question (c) is closed.
3. **VI.1, VI.2 and VI.3 are correct** and new in the F_2 model: a mod-2 derivative condition, Theta(m^2 log m)
   for complete Hadamard cores, and Omega(s^2) for full designs with translation B. They exclude construction
   classes (all twin-based fakes, complete cores, full designs) but say nothing about W(M) for unrestricted binary
   parents. **Growth scale over F_2 is open.**
4. **The filter grows (§3(f)).** No argument valid verbatim over F_3[t] with A^2+B^2 (all N) proves the uniform
   bound. The remaining algebraic escapes are characteristic 0, squarefree or bounded-multiplicity N in odd
   characteristic, p ≡ 1 mod 4, and ramified ("geometric") quadratic models.
5. **For the integer problem nothing changes:** the growth rate, the uniform bound, and Questions A and B. Parts V-VI
   are **not** evidence against the conjecture over Z (§3(g)). Score for Part VI on the goal: **2/10**.

## 1. Referee-corrected statements per block

### 1.1 ternary: VI.4 (referee: CORRECT; certificates corrected and improved)

**Proved**
* *Construction.* Let q = 3^r > 2^{m-1} with r odd, k = ⌈log_3(2m+1)⌉, V = { i·Σ_{b<k} x_b t^b : x_b ∈ F_q },
  n = q^k and h = n/q. Choose m pairwise nonproportional F_3-linear forms λ_j and signs ε_j, and put
  F_j(X) = Π_{v∈V}(X-v)^{1+χ(ε_jλ_j(v))}.
  * *Descent.* The 9-power map on coefficients fixes each F_j, so F_j ∈ F_9[t,X].
  * *Norms.* On F_9 the conjugation is the q-power map (r odd). It sends v to -v, and χ(-1) = -1, so
    F_j·F̄_j = P_V(X)^2.
  * *Distinctness* follows from the distinct kernels.
  * *(VI.5).* F_j = A(A^{(q-1)/2} - c^{(q-1)/2})^2 = P_V + 2c^{q-1}A - 2c^{(q-1)/2}A^{(q+1)/2}, with
    A = P_{ker λ_j} F_q-linearised. Hence pairwise X-degree <= (n+h)/2.
  * *(VI.6).* The expectation argument gives signs with common exponent degree g >= 2^{1-m}n > h. The condition
    q > 2^{m-1} is sufficient, not necessary.
  * *Specialisation.* X = t^L with L = 2(k-1)n+1 keeps the degree intervals of different X-powers disjoint, and
    gives D = (n-g)L with every difference < D/2.

  Every step was checked line by line by the ternary referee and again by the port referee.
* *Legitimacy.*
  * F_9(t)/F_3(t) is the constant extension. Its degree-1 place at infinity is inert, and A^2+B^2 is anisotropic
    there.
  * The points are monic of equal degree, so no two are unit multiples, and no conjugate of a point is in the
    cluster.
  * Dividing by the full gcd only tightens the condition (ff.md N2), so **no gcd cost is discarded**: the theorem
    uses the formal G, and the full gcd is larger.
  * Inert primes lie in the gcd. In the four-mark example every one of the 353 primes of N has even degree.
* *Other p.* The proof works verbatim for every p ≡ 3 mod 4 (referee: 3 points at D=174 for p=7).
* *p ≡ 1 mod 4: the mechanism fails formally (referee's lemma, proved).* Every half-system S of F_q^* has
  δ_S >= (q-1)/4. This forces g >= n/2+h/2, while two kernels give g <= n/2+h/(2q). Whether
  sup_D M_p(D,0) is finite for p ≡ 1 mod 4 is **open**.
* *Padding corollary (referee; [J] checked).* Multiplying all points by a common monic e keeps c=0 as long as
  deg e <= D - 2·maxdiff. The slack grows linearly in L, so M_p(D,0) >= m for **every** D >= D_1(m) = exp(O(m log m)),
  not only along D ∈ (n-g)Z.

**Computed**
* *Four marks* (q=27, k=2, forms x, y, x+y, x-y, all ε=+, X=t^11): formal g=92, full gcd 1038 = 92·11 + 26 (the
  extra 26 is t^26), D=6981, parent degree 13962, differences 3471, 3457, 3457, 3471, 3471, 3457.
  * Reproduced by the ternary and port referees.
  * **[J] reproduced** with `ternary.py` (own F_729, own primitive modulus; `ternary_m4.out`).
  * All coefficients lie in F_9. Norms are equal and lie in F_3[t] (computed as A^2+B^2 for z=A+iB, and the
    q-power conjugation was checked coefficientwise).
  * The points are distinct, and no point equals the conjugate of another.
  * With the formal G alone (D=7007) the differences are 3497, 3483, 3483, 3497, 3497, 3483 <= 3503.5, so the
    **formal certificate already suffices at L=11**. The extra t^26 brings the residual differences to Astra's
    numbers.
  * At L=10 (full gcd, D=6344) the maximum difference is 3185 > 3172; at L=9 it also fails. So L=11 is least,
    as the referee found (`ternary_m3_and_L.out`).
  * [J] also ran the three-mark instance (x, y, x+y; g=183, residual gcd 13, D=5993, differences 2483, 2469, 2483 <= 2996.5).
    It is trivial against the record's M_3(2,0)=4, but it is the smallest run of the mechanism at q=27.
* *Five marks.* Astra's certificate (forms x, y, z, x+y, x+z; g=1301, L=67, D=1,231,594, maximum difference
  615,589 < 615,797) is **correct** by exact expansion (ternary referee).
  * **Correction:** it is far from minimal. For k >= 3 the full gcd carries about 17,500 of extra degree at the
    F_{q^2}-rational places, where fibres of v ↦ v(α) have q^{k-2} elements (an exact formula, matched against
    actual gcds).
  * With that extra gcd, Astra's forms already work at L=37 (D=662,632).
  * The forms x_2, x_1, x_1+x_2, x_0, x_0+x_1-x_2 with signs (+,+,-,+,+) give g=1379, and **5 points at
    D=568,386** (L=32, full gcd 61,470). This was checked independently by gcd, exact division and norms.
  * **[J]** `gk3.py` reproduces both formal values, g=1301 and g=1379 (own F_27). The expansions were not redone.
* *Six marks* (referee, new): 6 points over q=27 although q < 2^5 (g=833 > h=729).
  * Certified with the formal G alone at D=6,880,250 (L=365, exact expansion).
  * D=3,582,084 (L=191) via the full-gcd formula only, so it is not independently confirmed.
  * q=27 cannot give 7 points.

**Corrected statements (structure of N; ternary, meaning).**
1. *Multiplicities.* N is far from squarefree.
   * In the four-mark example 338 primes are squared and 15 are simple; the radical has degree 7124 of 13962.
   * For k >= 3 the F_{q^2}-rational places carry multiplicities 4, 27, 31 (about q^{k-2}), so **multiplicities
     are unbounded as m grows**.
   * meaning: for this V, multiplicity q^{k-2} is forced at F_{q^2}-places for every polynomial specialisation
     X=T(t).
   * No squarefree or bounded-multiplicity variant is visible. The degree drop in (VI.5) needs the doubled
     exponents, and the {0,1} half-system version loses it.
2. *Residue fields.* 99.3% of deg N lies on primes of degree divisible by 6, whose residue fields contain
   F_{q^2} = F_729. In general the bulk of N sits on residue fields of size >> qL, which grow with m. Only the
   constant field is fixed.
3. *Rate.*
   * Astra's L gives log D ≈ 2k log q, hence **M_p(D,0) >= ((1/2)log_2 p - o(1)) log D/loglog D**. This is proved.
   * The constant log_2 p (with L ≈ (k-1)(n-h)/(g-h), q >= 2^m) is **plausible only**: distinctness below
     L = (k-1)(n-g) must then be checked separately.
   * In the model's units (W = deg N ↔ log R) this is m ≍ log W/loglog W, i.e. **loglog R/logloglog R**.
4. *Port.* Astra's "transfer" paragraph counts coefficient-Frobenius descent among the characteristic-p inputs. Port
   corrects this: descent is Galois descent and has an integer analogue (from O_K[i]). The single non-portable
   input is one parent with many same-degree right composition factors (§3(e)).

**Gap.** None in VI.4.

### 1.2 squarefree: VI.5 (referee: CORRECT; sharpness scope corrected)

**Proved**
* *Theorem (VI.5).* For d >= 3, m = 2^{d-1} and every integer L >= 4m(d-2), there is a squarefree N ∈ F_2[t] of
  degree 2(m-1)L carrying m points at c=0. The least degree is 8m(m-1)(log_2 m - 1). The referee checked line by
  line:
  * V = {a+(f^2+f)^2 : a ∈ F_2, f ∈ U} with U = span(t, ..., t^{d-1}). Every element is a square, with trace zero
    at even-degree primes.
  * (VI.7): deg R <= 2m(d-2).
  * CRT trace prescription Tr(T) = 1 + Tr(θ) = 1 + (e/2 mod 2), with deg T = L.
  * Pairwise coprimality of the factors T+θ+v, including against conjugates, and each factor squarefree.
  * Echelon bound deg P_J(1) <= 4(d-2)2^{d-2} <= L/2.

  The bound is attained, so the threshold L >= 4m(d-2) is exact for this V.
* *Sharpness pairing.* The marks realise exactly the m-1 nonconstant Walsh-Hadamard columns with disjoint block
  supports. With VI.2, the least degree **within complete Hadamard cores** is Theta(m^2 log m), with constants about
  24 apart (upper 4m^2 log_2 m, lower m^2 log_2 m/6). **Scope correction:** this is a statement about one profile
  class, not about M_2(D,0).

**Computed**
* CRT moduli of degree 6, 26, 60, 176; parents of degree 96, 896, 5760, 31744 (m=4, 8, 16, 32); 650 pairs. Also
  m=64: deg R=396, deg N=161,280, 2016 pairs (referee; also reproduced by bounds).
* **[J] reproduced m=4** with `sqf.py` (`sqf.out`):
  * V = {0, 1, t^4+t^2, t^4+t^2+1, t^8+t^2, t^8+t^2+1, t^8+t^4, t^8+t^4+1}; R = (t^2+t+1)(t^4+t+1), deg 6.
  * Exactly 16 of the 64 residues S (deg S < 6) meet both trace conditions. For every one of them, at L=16 and
    at L=17, the four points have one norm N ∈ F_2[t] of degree 96 (resp. 102), N is squarefree
    (gcd(N,N')=1), the points are distinct, and the pair differences are 20, 24, 24, 24, 24, 20 <= 24 = D/2
    (resp. 21, 25, ... <= 25.5).
  * At L=15 the pairs reach 23 > 22.5, so the threshold 4m(d-2) = 16 is exact.
  * Control: of the 48 residues violating a trace condition, only 12 give squarefree N.

**Gap.** None. Bookkeeping note (referee): as an F_2[t] product, N is the product of the m-1 block norms, not of
the norms of all 2m-2 factors (which would square it).

### 1.3 code: VI.1 (referee: CORRECT; integer analogue FALSE)

**Proved**
* *VI.1.* On every genuine F_2[t] parent with the norm of F_4(t)/F_2(t), every nonzero word of the anchored binary
  span of the split-prime parity rows has weighted weight >= K+1 > W/4, where K = ⌈D/2⌉. **[J] I read the proof
  and agree:**
  * R ∈ U_K gives ord_∞(R'/R) >= K+1.
  * In characteristic 2, (π/π̄)'/(π/π̄) = p'/p.
  * Inert primes and even exponents drop out, so R'/R = H'/H with H squarefree.
  * H' ≠ 0, because a squarefree nonconstant polynomial over the perfect field F_2 is not a square.
  * ord_∞(H'/H) = deg H - deg H' <= deg H.

  Sharpening (code): ord(R'/R) is even, so the weight is >= 2⌈(K+1)/2⌉ + deg H'. This is attained (K+2) on Part V.
* *Packing.* r <= (1-H_2(1/8))n + O(log n) < 0.457n, valid for arbitrary weights (Markov with radius W/8).
* *Exclusions over F_2[t], at any degree:*
  * the squarefree skew (Paley) double core for large m;
  * every fully flipped Hadamard core (any number of copies);
  * their residue designs;
  * single Paley cores (q ≡ 3 mod 8 from q=11, q ≡ 7 mod 8 from q=311).

  Passing: Walsh/Sylvester cores (rank log_2 m, the genuine Part V shape), t >= 3 permuted Paley bases, and
  low-rank two-core profiles.
* *Low-rank test (Astra).* Two Walsh cores on x and x^3 pass the weighted test. The polynomial condition forces
  e >= 4d-2, so W >= (2m-2)(4d-2) = O(m log m) is still allowed (code verified e_min = 4d-2 for d = 3, 5, 7, 9).

**Computed.** Rank 166 on 334 columns at m=167, with certificates at m=167, 191 and 251, reproduced independently.
Astra's script imports `ii_paley` from another folder and does not run from its own directory (cosmetic). The Part V
clusters (m = 4, 8, 16) have minimum weighted weight 8, 92, 556 against K+1 = 7, 57, 361.

**FALSE (a reading posed by the referees, not Astra's claim): "VI.1 has an integer analogue".**
* *Code (e) counterexample.* Take z_1 = γ^2(2+i), z_0 = z̄_1 with γ = x+iy and (x+2y)^2 - 5y^2 = ±1. Then
  |z_1 - z_0| = 2|Im z_1| = 2, and the parity word of z_1/z_0 is supported at 5 alone, while W = log N grows
  without bound.
* **[J] ruling: concur; decisive computation redone** (`pell.py`). For 12 consecutive Pell solutions,
  |z_1 - z_0| = 2 exactly, N = 5|γ|^4, and W = log N grows from 1.6 to 65.0. The parity support of z_1/z_0 is
  {5} (γ^2 contributes even valuations only).
* The mod-2 support condition therefore has no integer counterpart. The honest integer analogue (code, under abc
  over Q(i)) is a **Z-support** inequality:
  w_0(v) >= W/4 - (δ/(2(1+δ)))w_1(v) - log(|c|_1 C) - O_δ(1).
  This is the record's Lang-Waldschmidt support test. The skew double core passes it for every character
  (`round8/lw.md` Thm 5.1). It excludes the fully flipped families with b >= 4 (w_0/W = 0.222 at b=4), as LW
  Theorem L already does. Szpiro excludes none of them.

**FALSE (meaning (e); a reading, not Astra's claim): "VI.1 is Mason in disguise, so abc would exclude the integer
fakes".** **[J] concur**, for the same reason. VI.1 is a *mod-2* Mason inequality: d/dt kills squares, so only odd
support counts. Characteristic-0 Mason and abc over Z give only the support condition above, which the double core
passes. The record already found abc vacuous on squarefree blocks (`growth/conditional.md` Prop 8.3,
`round8/weakest.md`). "The fakes exist because Z has no derivation" should not be written.

### 1.4 bounds: VI.2 and VI.3 (referee: CORRECT; scope corrected)

**Proved**
* *VI.2 exact products.* Π z_i^{χ_v} = r_v^{m/2} and Π z_i^{(χ_v-χ_w)/2} = (r_v/r_w)^{m/4}, including units and
  common factors. This is the function-field form of the record's integral character.
* *Depth calculus.* depth(r^{2^a b}) = 2^a depth(r) for odd b. It forces the odd part of m to be 1, with
  depth(r_v) >= D/m and depth(r_v/r_w) >= 2D/m.
  * **Corollary (bounds, new):** complete Hadamard cores of order not a power of 2 (Paley, orders 12, 20, ...)
    never occur over F_2[t].
* *(VI.2).* D >= m(m-1)(log_2(m-1)-3)/6. The bound is vacuous for m <= 8 and gives 36.3 at m=16. The hypotheses
  actually used are:
  * (H1) the sign matrix plus the all-ones column is Hadamard;
  * (H2) m >= 4;
  * (H3) the f_v are monic, nonconstant, pairwise distinct, and f_v ≠ f̄_v;
  * (H4) the block products are critical at D = Σd_v.

  Disjoint supports are only one sufficient condition for these, and Part V (supports not disjoint) satisfies
  (H1)-(H4).
* *(VI.3).* D' >= (1/2)Σ_{variable pairs} h_v >= (Q/8 - 1/2)deg B for all-form hyperplane designs with
  translation B. The second half gives D' = Omega(s^2) via Theorem E at inert primes with P_J(B).
  * Every step was checked, including a missing one-line proof (supplied by the referee) that the reduced numerator
    is coprime to its conjugate.
  * 0 failures on 618,000 instances for the first half, and on 62,964 critical instances for the identity, the
    inert bound and the margin identity.
* **Corrected scope.**
  * VI.2 needs a **complete** orthogonal core. With k <= m/2 Hadamard columns the method gives only D >~ k log k,
    which is growth scale, so the quadratic cost comes entirely from completeness.
  * VI.3 needs **all** Q/2 forms. It is purely asymptotic: no threshold is given and the implied constant is about
    1/384. It proves s^2, not s^2 log s.
  * Sparse subfamilies with B divisible by the small odd-degree primes escape both halves of VI.3.

**Computed (bounds, new; improves the record's F_2 constants)**
* Nonconstant B beats B=1 by a constant factor of about 0.7:
  * 8 points at point degree **68** (B=t^2+t, Y=t^11+θ(t^2+t)) against 104;
  * 16 points at 472 (B=t^3+t, T=t^33) against 672 (and 720 for Part V);
  * 32 points at 2607 (B=t^4+t^2, T=t^86) against 3646 (and 3968).

  Every winning B is divisible by both degree-1 inert primes t and t+1, which is exactly the escape that the VI.3
  dichotomy predicts.
* **[J] reproduced the 8-point case** (`b8.py`, own F_4[t]). With V = F_2[t]_{<4}, all 8 forms with ℓ(B)=1, and all
  128 side choices modulo global conjugation, the least c=0 cluster has reduced point degree 68:
  * sides (0,1,0,1,1,0,1,0), full gcd degree 20, maximum difference 34 = D/2;
  * reduced points distinct and monic, with equal norms in F_2[t] of degree 136.
* Minor record inaccuracy (bounds): `ref5_lower`'s d=3 table omits the unit-twisted case T=0, where 4 points occur
  at degree 4.

**Gap (scope; not an error).** Incomplete orthogonal cores, sparse hyperplane subfamilies, non-hyperplane
transversals, independent per-prime choices (only VI.1 constrains them, and low-rate codes pass), and unions of
designs are all untouched. Whether full designs with B over a general V cost Theta(s^2) is a **conjecture**.

### 1.5 port: can VI.4 be ported? (analyst: no integer counterpart; failing step located)

**Proved**
* *Lemma 1 (port; elementary, classical in substance: Ritt/Engstrom-type uniqueness of right composition
  factors).* If char k does not divide deg f and f(A) = f'(A') with deg A = deg A', then A' = αA + β. So in
  characteristic 0 one parent cannot have m different same-degree right composition factors. In VI.4 the parent
  P_V = A_j^q - c_j^{q-1}A_j has one for every kernel, with f q-linearised (f'_Y constant). **[J] I checked the argument.** Let d = deg f be invertible. By a Tschirnhaus shift,
  f(Y) = c(Y+e)^d + (terms of degree <= d-2). Then P/c = (A+e)^d(1 + O(t^{-2 deg A})), so the polynomial part of the
  Laurent-series root (P/c)^{1/d}, which exists because char ∤ d, is A+e. The same holds for A' with ζ^d = c/c', so
  A' is affine in A. In VI.4, d = q is the characteristic power and the root does not exist. That is the exact step
  that separates characteristic p from characteristic 0.
* The char-0 version of the identity yields one reflection pair {z, -z̄}. Every C(t) cluster has at most 24 points
  (Part IV, cited).
* The other ingredients each exist in characteristic 0 or over Z: the character exponents, descent, common
  division and specialisation.

**Computed (routine).** None of the integer analogues produced a cluster:
* one odd character gives the trivial pair;
* several characters and the literal digit-vector port give width exponents 0.94-1.00 after exact Gaussian gcd
  removal, where a cluster needs 1/2;
* sign-pattern templates give at most 5 points (D=6..12);
* number-field roots and multi-parameter polynomial families reduce to one-parameter Z[i][u] templates, so they are
  <= 24 by Part IV.

The only escape left is non-polynomial families, which is the general genuine-angle problem itself.

**Corrected.** Astra's transfer paragraph is right in its conclusion (no integer construction and no growth
theorem). The non-portable input is "one parent, many same-degree decompositions", not Frobenius descent.

### 1.6 meaning: against the record and the literature (analyst)

**Proved**
* *(a) Transfer.* The record's growth bound transfers to F_3[t] with A^2+B^2:
  W >= 8F_3*(M)/M - 4(M-1)c = 3M log_3 M - O(M). Odd-degree primes are inert, even-degree primes split, Hensel
  lifting works, and no place ramifies.
  * Exact minima W_min = 4.8, 151.5, 440.2, 2516.1, 31326.2 at M = 5, 32, 64, 256, 2048. The same code reproduces
    the F_2 values.
  * The ternary clusters sit **exponentially** above this bound (W = exp(Theta(m log m))), the binary ones
    polynomially (W ~ m^2 log m).
* *(c) Questions A/B over F_3(t).* By Weil, every 5-subtuple of the ternary clusters is balanced:
  δ_S = 1/30 - 8/(15q) + O(q^{-3/2}) and H/w ~ 80/q. Checked numerically at q=243. Hence the F_9(t)/F_3(t)
  analogues of A_5(δ<1/30), (L_5°) and V_E(5) are false. This uses J4 + Schwartz-Zippel, which is
  characteristic-free but was not rechecked in characteristic 3. (L_8) fails asymptotically only, with the
  dictionary unchecked in characteristic 3.
* *(b) Escape clauses.*
  * Closed for all N: odd characteristic (p ≡ 3 mod 4) and reflection conjugation.
  * Closed only in characteristic 2: squarefree N.
  * Lemma P is still true but no longer a barrier.

**Not checked.** (d) Literature: no prior art was found for worst-case short-arc failure over F_q[t], q ≡ 3 mod 4
(search snippets only). The ingredients are classical: Ore linearised polynomials, Ben-Sasson-Kopparty-Radhakrishnan
subspace polynomials, and Π_{squares}(Y-s) = Y^{(q-1)/2} - 1. The judge did not search further.

### 1.7 Withdrawn or corrected in Part VI (collected)

1. Five-mark certificate: correct, but not minimal (568,386 is achievable). It ignores the extra full gcd at
   F_{q^2}-places for k >= 3.
2. The theorem's rate is (1/2)log_2 p · log D/loglog D. The constant log_2 p needs a separate distinctness check.
3. "Parent degree exp(O(m log m))": correct. Add that N is square-full in bulk, with multiplicities growing like
   q^{k-2} and the bulk of N on residue fields of size >> 2^m.
4. VI.5 "sharpness for the profile in VI.2": correct within complete Hadamard cores only, and conditional on VI.2.
5. VI.2/VI.3: correct. Their scope is complete cores and full designs. VI.3 is asymptotic with no threshold and
   gives s^2, not s^2 log s.
6. Transfer paragraph: descent is not characteristic-p specific.
7. `vi_profile_check.py` imports from another folder (cosmetic).

**No claim of Astra's is withdrawn.**

## 2. What is new relative to the record, and correct

* **VI.4: the first unbounded family in odd characteristic with a reflection-type conjugation**, at a fixed
  constant field, with the Gaussian form A^2+B^2 and an inert infinite place. It answers README_V §5 item 1 (ranked
  first there) **negatively**. It escapes Lemma P through nonconstant cofactor differences and 0/2 character
  exponents. grep finds no quadratic-character/linearised construction in the record. The ingredients are
  classical; the combination appears new (meaning, search snippets only).
* **VI.5: squarefree N over F_2** (README_V §5 item 3, question (c)), via a space of squares with trace zero at
  even-degree primes plus a CRT trace condition. New.
* **VI.1: a mod-2 derivative/parity-code condition** on every binary parent. It excludes every twin-based fake
  profile of the record from the binary model and explains the low rank (log_2 m) of genuine Part V clusters. New;
  the packing step is standard.
* **VI.2** (Theta(m^2 log m) for complete Hadamard cores; no complete core of order not a power of 2 over F_2[t])
  and **VI.3** (Omega(s^2) for full designs with any translation B). New in the model, and both close construction
  classes.
* **Referee by-products, new:**
  * 5 points at D=568,386 and 6 points at D=3,582,084 (formal-only 6,880,250) over F_3;
  * the exact full-gcd formula for k >= 3;
  * the padding corollary (all large D);
  * the p ≡ 1 mod 4 half-system obstruction δ_S >= (q-1)/4;
  * Lemma 1 as the failing step in characteristic 0;
  * the transferred F_3 growth bound;
  * smaller F_2 clusters with nonconstant B (68/472/2607);
  * the integer Pell counterexample to a parity analogue;
  * abc on character relations = the LW support test.

Not new: subspace polynomials (BKR, already cited in Part V), Lemma P (Part V), the LW support test and Theorem L
(round 4-8), and the margin identity.

## 3. Rulings

### (a) Is "sup_D M_3(D,0) is infinite for A^2+B^2 over F_9(t)/F_3(t)" a THEOREM?

**Yes.** Statement for the record:

> **Theorem (round 31, Part VI.4; referees ternary and port; judge [J] reproduced the four-mark instance).** Let
> p ≡ 3 mod 4, F_{p^2} = F_p(i), conjugation i ↦ -i on coefficients (the norm form A^2+B^2 of the constant
> extension F_{p^2}(t)/F_p(t), anisotropic, with the infinite place inert). Let M_p(D,c) be the largest number of
> distinct z_j ∈ F_{p^2}[t] of degree D with one norm z_j z̄_j = N ∈ F_p[t] (any multiplicities) and
> deg(z_i - z_j) <= D/2 + c. Then for every m >= 2 there is D_1(m) = exp(O(m log m)) such that M_p(D,0) >= m for
> every D >= D_1(m). Hence sup_D M_p(D,0) = ∞ and
>   **liminf_{D→∞} M_p(D,0) / (log D/loglog D) >= (1/2)log_2 p**.

*Exact rate.* The proved constant is (1/2)log_2 p (≈ 0.79 for p=3). log_2 p is plausible: it needs a distinctness
check below L = (k-1)(n-g).

*Upper side* (meaning; transferred 2/3 argument, not redone by the judge): W = deg N >= 3M log_p M - O(M). The gap is
**exponential**: in model units the construction gives m ≍ log W/loglog W, and the bound allows W/log W.

*First instances (p=3):* 4 points at D=6981 ([J] reproduced), 5 at D=568,386, 6 at D=3,582,084 (formal-only
6,880,250). The exhaustive table M_3(D,0) <= 4 for D <= 11 is no guide.

*Caveats, all mandatory when citing:*
1. *Multiplicities of N.* N is not squarefree. The generic prime pattern is π^2/π̄^0, and for k >= 3 the multiplicities
   at F_{q^2}-places are about q^{k-2}, **unbounded in m**. Squarefree or bounded-multiplicity N in odd characteristic
   is **open**, and this mechanism cannot reach it.
2. *Residue fields.* The constant field is fixed (F_9/F_3), but the bulk of deg N (99.3% in the four-mark example)
   sits on primes whose residue fields contain F_{q^2}, q > 2^{m-1}. So the result does not show that *small*
   residue fields fail to constrain clusters. It shows that fixing the constant field does not suffice.
3. *Common divisor.* The mechanism needs a common divisor of degree g >= 2^{1-m}n, which is most of the raw degree's
   slack. Dividing by it is legitimate and only tightens the cluster. The full gcd (formal G plus t^26 at k=2, plus
   about 17.5k at F_{q^2}-places for k=3) is larger still, and helps.
4. *Only p ≡ 3 mod 4* and only the constant extension. p ≡ 1 mod 4 fails formally (referee's lemma), and ramified
   ("geometric") quadratic models are untested.

### (b) Is the squarefree binary statement a theorem?

**Yes.**

> **Theorem (round 31, Part VI.5; referees squarefree and bounds; judge [J] m=4).** For d >= 3, m = 2^{d-1}, and
> every integer L >= 4m(d-2), there is a squarefree N ∈ F_2[t] of degree 2(m-1)L carrying m distinct points of
> F_4[t] (norm of F_4(t)/F_2(t)) at c=0. Hence sup_D M_2^{sf}(D,0) = ∞, with least supplied parent degree
> 8m(m-1)(log_2 m - 1), four times Part V's.

The threshold is exact for this V ([J]: L=15 fails at m=4). Rate m ~ sqrt(W/(4 log_2 W)) roughly, i.e. the same
order as Part V and far above the transferred bound 3m log_2 m.

### (c) The code condition (VI.1)

* **Correct** (code; [J] read the proof). Every nonzero word of the anchored parity span has weighted weight
  > W/4, so the rank is < 0.457n.
* **What it excludes**, over F_2[t] at every degree and for every prime assignment: the record's large squarefree
  skew (Paley) double core, all fully flipped Hadamard cores and their residue designs, and single Paley cores
  beyond small q.
* **What it does not exclude:** Walsh/Sylvester-type (low-rank) profiles, which is the shape of all genuine F_2
  clusters, and low-rate codes in general. It leaves W(M) ∈ [3M log_2 M - O(M), O(M^2 log M)] untouched, and the
  two-Walsh-core test still allows W = O(m log m).
* **Integer analogue.** There is **no** mod-2 analogue: the Pell family (§1.3; [J] recomputed) has parity support
  {5} with W → ∞.
* The analogue under abc over Q(i) is the Z-support inequality of §1.3, i.e. the Lang-Waldschmidt support test. It
  excludes fully flipped families with b >= 4 (as LW Theorem L already does under LW) and **not** the skew double
  core. Szpiro excludes nothing.
* **No new conditional growth-rate statement is worth recording.** The only addition is a one-line refinement of
  barrier 6: *abc applied to character relations (products of points, not of chords) is equivalent to the LW support
  test, and gives no growth improvement.* Whether simultaneous angle models survive abc outside the character lattice
  is heuristic only.

### (d) Scope of VI.2 and VI.3; growth scale over F_2[t]

* VI.2 is a theorem about **complete** Hadamard cores: Theta(m^2 log m), with no complete core of order other than a
  power of 2. VI.3 is a theorem about **full** hyperplane designs (all Q/2 forms) with any translation B:
  Omega(s^2), asymptotic, unspecified threshold.
* Both are statements about one construction class each, and **neither bounds M_2(D,0)**. The quadratic cost comes
  from the rigidity of completeness: with k <= m/2 columns the VI.2 method gives only k log k.
* **Status of growth scale over F_2[t]: open.**
  * Known: 3M log_2 M - O(M) <= W(M) <= 2M^2 log_2 M (proved, Part V). Empirically about 0.7 of that with B ≠ 1 (W = 5214 ≈ 1.0·M^2 log_2 M at M=32).
  * Every known family is at M^2(log M)^{0 or 1}. Nothing excludes M^{1+o(1)} beyond the transferred bound.
  * The open classes are incomplete Walsh cores, sparse hyperplane subfamilies with B divisible by the small inert
    primes, non-hyperplane transversals, independent per-prime choices, and unions of designs.
  * The judge's estimate that growth-scale binary clusters exist is low (well below even), with no evidence either
    way [J].

### (e) Did the port analyst find an integer counterpart or an obstruction?

**No integer counterpart. A precise obstruction, at the level of mechanism:**
* VI.4 is a *polynomial template* over the characteristic-3 coefficient field F_9(t), specialised at X=t^L.
* Its integer counterpart class is one-parameter (or weighted multi-parameter) templates in Z[i][u], including roots
  in totally real fields with Galois descent. Part IV bounds that class by 24.
* Its single non-portable step is **one parent with many same-degree right composition factors**. Lemma 1
  (Ritt/Engstrom-type) makes this impossible in characteristic 0. Equivalently, the Newton linearisation used in
  Part IV cannot be inverted in characteristic p.
* The other ingredients (descent, common division, specialisation, character exponents) all have integer
  analogues and do not help.
* Every computed integer analogue has width exponent 0.94-1.00, against the 1/2 needed.
* This is an obstruction to *porting the mechanism*, not a theorem about integer clusters: non-polynomial
  (genuine-angle) families are the general problem and are untouched.

### (f) The filter as it now stands

> **F3++ (judge's wording, combining meaning, port and the Part V filter F3+).** Call an argument
> *function-field valid* if every step holds verbatim, with deg in place of log and the inert degree-one place at
> infinity in place of the archimedean place, in **either** of the following settings:
> * (i) F_2[t] with the anisotropic norm of F_4(t)/F_2(t), for all N, **or for squarefree N only**;
> * (ii) F_p[t], p ≡ 3 mod 4, with the Gaussian norm A^2+B^2 of F_{p^2}(t)/F_p(t), for all N (prime powers allowed).
>
> Then **no function-field-valid argument proves the uniform bound M <= M(C).** Quantitatively:
> * an argument valid in setting (i) cannot prove M <= f(W) with f(W) = o(sqrt(W/log W)), even restricted to
>   squarefree N;
> * an argument valid in setting (ii) cannot prove M <= f(W) with f(W) = o(log W/loglog W). It can still prove, for
>   example, M << W^ε.
>
> Admissible (hence insufficient) steps include unique factorisation, exact norms, anisotropy and inertness at
> infinity, a reflection-type conjugation (setting ii), a fixed odd characteristic (setting ii), squarefreeness
> (setting i), residue counts and pigeonhole at small primes, cut/chord/Plücker identities, the chain lemma,
> Liouville with toric forced divisibility, Schwartz-Zippel, and Lemma P.

**What a proof over Z can still use** (none is known to suffice):
* characteristic 0 in an essential way. Part IV uses it through Newton linearisation and BMY, and Lemma 1 is its
  combinatorial shadow;
* archimedean order, curvature and real-place balance;
* Baker/transcendence input;
* squarefree or bounded-multiplicity N **in odd characteristic** (open, and VI.4's mechanism cannot reach it);
* features of ramified quadratic extensions (Q(i)/Q ramifies at 2; constant extensions are unramified everywhere);
* the absence of a coefficient Frobenius.

Growth-rate arguments are **not** filtered at the integer growth scale: (ii) filters only below log W/loglog W, and
(i) only below sqrt(W/log W). The record's W >= 3M log M transfers in both settings, and no function-field family
reaches it. Number-field Vojta and the implications A_5 ⇒ uniform and (L_8) ⇒ uniform are untouched.

### (g) What changes for the integer problem

Bluntly: **nothing**.
* **Growth rate:** unchanged at (2/3+ε) log R/loglog R, with the same constant 3. No Part VI step is valid over Z
  and gives a better order. VI.1's mod-2 gain does not exist over Z (Pell), and the integer analogue under abc is
  the known LW test.
* **Uniform bound M(C):** unchanged. Nothing is proved or disproved for any integer parent.
* **Questions A and B:** unchanged over Z. Their F_3(t) analogues (A_5(δ<1/30), (L_5°), V_E(5)) now fail too, with
  all 5-subtuples balanced, which only says that a proof must use number-field or archimedean input.
* **Are Parts V-VI together evidence against the uniform conjecture over Z?**
  * *For "some evidence":* every algebraic feature of Z[i] that holds verbatim in F_9[t] is compatible with unbounded
    clusters: a fixed odd characteristic, reflection conjugation, the Gaussian form, an anisotropic inert infinite
    place, exact norms, and unique factorisation. So does squarefreeness, in characteristic 2. This removes the hope
    of a "purely algebraic" proof and moves the burden onto characteristic-0 and archimedean input, which has never
    yet been made to work over Z.
  * *Against reading it as evidence (decisive in the judge's view):*
    1. Both mechanisms rest on additive (q-linearised) polynomials, which have no characteristic-0 counterpart. The
       characteristic-0 function-field model, which shares the reflection, the anisotropy, exact norms and
       isotriviality with Z[i], is **bounded** (Part IV, <= 24).
    2. The exact failing step is identified and impossible in characteristic 0 (Lemma 1). Every integer analogue
       of the mechanism reduces to a class that is already bounded.
    3. Characteristic-p counterexamples of this genre (Mahler for Roth; Lang and Easton for BMY) carry no weight on
       their characteristic-0 statements.
    4. The ternary clusters live at m ≍ log W/loglog W with square-full N and growing residue fields. That is
       nowhere near the integer growth scale and nothing like an integer parent.
  * **Verdict: not evidence against the conjecture over Z.** Parts V-VI are a methodological filter only.

## 4. Sentences of the expert note that should be updated (recommendations only; the record is read-only)

1. §5 item 12, header and the closing table row (line 865): add "Part VI (six blocks, `round31/README_VI.md`):
   **unbounded over F_3(t) with A^2+B^2**; squarefree unbounded over F_2(t)".
2. §5 item 12, Part V bullet *Theorem*: replace "N is not squarefree; the squarefree variant is **open**." with
   "N is not squarefree here; squarefree parents of degree 8m(m-1)(log_2 m-1) carrying m = 2^r points exist
   (Part VI.5, **proved**)."
3. §5 item 12, Part V bullet *Filter F3+*: replace "Not affected: ... arguments needing odd characteristic, a
   reflection-type conjugation, or squarefree N in an essential way." with "Not affected: growth-rate arguments at
   the integer scale; number-field Vojta and the implications A_5 ⇒ uniform, (L_8) ⇒ uniform; arguments using
   characteristic 0, archimedean order, Baker input, ramification, the absence of a coefficient Frobenius, or
   squarefree/bounded-multiplicity N in odd characteristic. Odd characteristic and reflection conjugation do not
   suffice (Part VI.4), and neither does squarefree N alone (Part VI.5, characteristic 2)." Then add F3++ from §3(f),
   including its two quantitative thresholds.
4. §5 item 12, Part V bullet *Not evidence against the conjecture over Z*: replace "In characteristic ≠ 2 with a
   reflection-type conjugation ... (analyst's Lemma P ...). Odd characteristic is **open**: nothing was constructed
   and nothing excluded." with "Lemma P (proved) explains only why the *characteristic-2* mechanism does not port.
   Over F_p(t), p ≡ 3 mod 4, with A^2+B^2 the model is unbounded (Part VI.4): quadratic-character multiplicities on
   the roots of a q-linearised subspace polynomial, coefficient-Frobenius descent and a large common divisor;
   M_p(D,0) >= ((1/2)log_2 p - o(1)) log D/loglog D; N square-full in bulk, multiplicities unbounded in m. The
   non-portable step is one parent with many same-degree right composition factors, impossible in characteristic 0
   (Ritt/Engstrom-type Lemma 1); the integer counterpart class (polynomial templates) is bounded by Part IV."
5. §5 item 12, *Theorem (function-field uniform bound)*: "Characteristic 0 is essential: the statement is false
   over F_2(t) with an anisotropic norm (Part V)" → add "and over F_p(t), p ≡ 3 mod 4, with A^2+B^2 (Part VI)".
6. §5 item 12, Part V bullet *Rate*: "M_3(D,0) <= 4 for D <= 11" → add "; nevertheless sup_D M_3(D,0) = ∞ (Part VI;
   first 5-point instance at D = 568,386)".
7. §5 item 10, *Function-field model*, last sentence: "Over the fixed finite field F_2 the model is nevertheless
   unbounded ... by an additive mechanism with no integer counterpart" → "Over the fixed fields F_2 (Part V) and
   F_3 with A^2+B^2 (Part VI) the model is nevertheless unbounded, by additive (linearised-polynomial) mechanisms
   with no integer counterpart".
8. §3 Calibration (line ~332) and §2 Calibration (line ~193): add "and over F_9(t)/F_3(t) (Part VI; all 5-subtuples
   balanced by Weil)".
9. §5 barrier 6 (abc on single relations): add "abc on character relations (products of points) is equivalent to
   the LW support test (round 31 Part VI, code referee): it excludes the fully flipped families with b >= 4 and not
   the skew double core; the stronger characteristic-2 parity condition (VI.1) has no integer analogue (Pell
   family z = γ^2(2+i), |z - z̄| = 2)."
10. Do **not** write: "false in positive characteristic" (p ≡ 1 mod 4 and ramified models are untested), "false
    for sums of two squares over function fields" (the A^2+tB^2 model is untested), "evidence against Z", "the
    growth rate over F_q is log W/loglog W" (only a lower bound is known), or "the fakes exist because Z has no
    derivation".

## 5. Open sub-questions, ranked by honest probability of affecting the integer problem

All probabilities are the judge's estimates [J], and all are small.

1. **Squarefree (or bounded-multiplicity) N over F_3 with A^2+B^2.** Is sup_D M_3^{sf}(D,0) finite? A bounded answer
   would be the first mechanism-level statement that survives every function-field filter and is Z[i]-shaped. It
   is the last algebraic route left to a uniform proof. P ≈ 0.005.
   *First step:* extend the exhaustive M_3(D,0) search with a squarefree filter to D = 12-14. In parallel, search
   odd sign functions (half-systems) on V whose products drop about n/4 in degree between different kernels.
2. **A ramified ("geometric") quadratic model** F_3(t)(sqrt d(t)), deg d even, nonsquare leading coefficient and
   ramified finite primes, so that conjugation is not a coefficient Frobenius. It mirrors the ramification of
   Q(i)/Q. P ≈ 0.004.
   *First step:* exact small-D search for M(D,0) in that model, and test whether any linearised-polynomial identity
   survives a non-constant conjugation.
3. **Binary growth scale.** W(M) = O(M log M) over F_2, or a proof of W(M) >= ω(M)·M log M. A proof using only
   Z-valid steps would improve the integer rate. P ≈ 0.003.
   *First step:* construct actual polynomials for the two Walsh cores on x and x^3 at e = 4d-2, i.e. the case VI.1
   still permits. In parallel, run bounds' least-D search for incomplete Walsh cores at m = 8, 16.
4. **p ≡ 1 mod 4.** Is sup_D M_p(D,0) finite for A^2 - dB^2 over F_{p^2}(t)/F_p(t)? This is less Z[i]-like
   (-1 is a square in F_p), so it is a calibration only. P ≈ 0.001.
   *First step:* exhaustive M_5(D,0) for D <= 8, and a search for odd sign functions beating δ_S >= (q-1)/4 with
   non-character exponents.
5. **abc outside the character lattice.** Can simultaneous angle models survive abc on relations that are not
   characters? P ≈ 0.001.
   *First step:* formalise code's union-bound heuristic for the skew double core at P >> M/δ and test it on the
   record's angle models.
6. **Exact ternary rate** (log_2 p versus (1/2)log_2 p, and any upper bound below W/log W in characteristic 3).
   Calibration only. P ≈ 0.0005.
   *First step:* check distinctness of the specialised marks at L ≈ (k-1)(n-h)/(g-h) for q=243, k=3.

## 6. Scores (0-10)

| block | score | reason |
|---|---|---|
| ternary (VI.4) | 3/10 | Thorough line-by-line check, independent reproduction, a 2.2x smaller five-point instance, a new six-point instance, the exact full-gcd formula, the structure of N, the honest rate, and the p ≡ 1 mod 4 obstruction. Model-level only. |
| squarefree (VI.5) | 2/10 | Clean confirmation up to m=64, exact threshold, scope of the sharpness claim corrected. Closes a stated record question; no bearing on Z. |
| code (VI.1) | 2/10 | Correct theorem, a sharpening, the exclusion list extended, and the decisive Pell counterexample showing the gain is characteristic-2 only; abc ⇔ LW support test placed correctly. |
| bounds (VI.2, VI.3) | 2/10 | Both proofs checked with large computations; hypotheses isolated; new corollary (no non-power-of-2 complete cores); smaller F_2 clusters with B ≠ 1 (judge reproduced 68). Closes classes only. |
| port | 2/10 | Located the exact non-portable step (Lemma 1), corrected Astra's attribution to descent, and showed every integer analogue fails or reduces to Part IV. |
| meaning | 2/10 | F_3 transfer of the growth bound with exact minima; Questions A/B fail in the F_3 model for all 5-subtuples; precise wording rules; the Mason/abc misreading refuted. |
| **Part VI, on the goal** | **2/10** | Two correct new theorems in the function-field model: unbounded over F_3 with A^2+B^2, the most Z[i]-like failure yet, and squarefree unbounded over F_2. Plus three correct class exclusions over F_2. The filter now covers everything algebraic about Z[i] that holds in F_9[t]. Growth rate, uniform bound, and Questions A and B over Z are unchanged, and none of it is evidence against the conjecture over Z. |
