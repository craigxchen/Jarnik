# Referee report, round 31, Part VII — block "ramified"

Scope: the ramified model in VII.1 (answer.md lines 1467-1520) and the ramified tables of VII.3 (lines 1606-1613).
Scripts (all my own arithmetic; her scripts were read only after my results were in):
`/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7_ramified/`

* `fp.py` — F_3[t] arithmetic (numpy; schoolbook or float-FFT products), gcd, squarefree decomposition, the order F_3[t,η].
* `f27.py` — F_27 = F_3[a]/(a^3-a-1), quadratic character χ(y) = y^13.
* `construct_q3.py`, `run_q3b.py` (log `q3.log`) — VII.1 with q = p = 3, k = 2 (no constant extension at all).
* `construct_q27.py`, `run27.py`, `run27d.py` (logs `m3.log`, `m4.log`, `odd.log`, `odd2.log`) — VII.1 with q = 27, k = 2.
  Each Frobenius orbit of labels is multiplied symbolically through the elementary symmetric functions of
  w = x0 + x1 t^3 (asserted to lie in F_3[t]), so the marks are computed entirely in F_3[t,η]. The marks are
  **expanded exactly**. Every big product is cross-checked by an independent reduction modulo a random degree-40
  polynomial.
* `ram_census.cpp` (log `census1_8.log`) — exhaustive census M(D,0) for η^2 = 2t^2+1, D = 1..10.
* `brute_small.py` — an independent brute force for D <= 5. It uses all leading coefficients, groups points by exact
  norm and runs a genuine max-clique search on the pairwise width relation, with no ultrametric shortcut.

Verdicts: **VII.1 (ramified model): CORRECT**, with one clarification about non-principality and two scope additions.
**Ramified tables of VII.3: CORRECT** (reproduced exactly for D = 1..10). **Novel w.r.t. the record: yes. Significance
for Z: nil beyond the filter.**

---

## (a) Set-up

* K = F_p(t)(η), η^2 = Δ, deg Δ = 2δ > 0 and lc(Δ) a nonsquare. At infinity, η = sqrt(lc)·t^δ(1+O(1/t)) with
  sqrt(lc) ∉ F_p, so the infinite place is **inert** (residue field F_{p^2}). For a + bη the leading terms
  a_d t^d and b_e sqrt(lc) t^{e+δ} cannot cancel. So the pole order is max(deg a, deg b + δ) and deg N(a+bη) = 2·(pole order).
  A^2 - ΔB^2 is anisotropic at infinity. Checked; also checked numerically (the leading coefficient of N is
  1, 1, 2 for the three leading types, never 0).
* Constant field: the condition "a finite ramified prime" means Δ is not a constant times a square, so F_p is
  algebraically closed in K. Because r is odd, lc(Δ) stays a nonsquare in F_q, so F_q(t,η)/F_q(t) is still a quadratic field
  extension, and σ: η ↦ -η extends F_q-linearly. Correct.
* Conjugation σ fixes constants, so it is **not** a coefficient Frobenius. Correct: this is the point of the model.
* Non-maximal order (square factors in Δ): **irrelevant for lower bounds.** Every cluster in F_p[t,η] is a cluster
  in the maximal order with the same norms and degrees. For the census, Δ = 2t^2+1 = 2(t-1)(t+1) is squarefree, so the
  order is maximal.
* **Qualification the record should carry (not in Part VII).** With inert infinity the place ∞ has degree 2. The
  class group of the affine ring then fits into 0 → Pic^0(K) → Cl(O) → Z/2 → 0. So **|Cl(O)| = 2h is even, and O is
  never a PID.** Concretely, for η^2 = 2t^2+1 (a conic, genus 0) Cl(O) = Z/2. The ramified primes above t±1 have
  degree 1, and they are **non-principal**: every norm has even degree, so no element has norm c(t-1). In Z[i] the ramified prime (1+i) is
  principal. The constant model F_9[t] is a PID. So the even-degree "ramified with inert infinity" model mirrors
  ramification but loses unique factorisation. It is a legitimate test of "ramification + geometric conjugation",
  but it is not a closer copy of Z[i] in every respect. The construction does not care: N_0 is coprime to the
  discriminant, see (c).

## (b) Descent and conjugation, step by step

* V = {η^p C^p Σ x_j t^{pj} : x_j ∈ F_q} = η·Δ^{(p-1)/2} C^p·(F_q-span of t^{pj}). It has q^k elements, all p-th powers in
  F_q(t,η), because Σ x_j t^{pj} = (Σ x_j^{1/p} t^j)^p. Correct.
* Coefficient Frobenius φ (x ↦ x^p on F_q, fixing t, η) sends the label with vector x to the label with vector x^p.
  The forms are over F_p, so ℓ(x^p) = ℓ(x)^p and χ(y^p) = χ(y). Hence the exponents 1+χ(ε_jℓ_j) are φ-invariant and
  F_j, P_V, G are φ-invariant. They therefore lie in (F_q[t,η])^φ[X] = F_p[t,η][X]. Correct. In code: exponents are asserted
  constant on every φ-orbit, and the elementary symmetric functions of each orbit are asserted to lie in F_3[t].
* σ(v) = -v for every label, since η^p ↦ -η^p. So V is anti-invariant, and σ(F_j) = Π(X+v)^{1+χ(εℓ(v))} =
  Π(X-v)^{1+χ(εℓ(-v))} = Π(X-v)^{1-χ(εℓ(v))}, because χ(-1) = -1 (q ≡ 3 mod 4). Hence F_j·σF_j = P_V^2. Correct.
* P_V = A^q - c^{q-1}A and F_j = A(A^s - c^s)^2, s = (q-1)/2. I re-derived both; deg_X(F_i - F_j) <= h(q+1)/2. Correct.
* g = Σ min_j exponents: my count gives **g = 92** for the four forms x, y, x+y, x-y with all signs + (q = 27). For
  three forms the best g is **183**; for two forms over q = 3 the best g is 5 (> h = 3), and for 3 or 4 forms over q = 3
  g <= h. So q = 3 gives only two marks.

## (c) Collisions, simplicity, squarefree R

* Two labels equal at a point P over t = α with C(α) ≠ 0 would force (ηC)^p w^p(α) = 0 with 0 ≠ w ∈ F_q[t]_{<k}.
  Since Δ | C we have η(P) ≠ 0, so w(α) = 0, so α ∈ F_{q^j} for some j <= k-1, so C(α) = 0. Contradiction. At zeros of C all labels vanish
  (they are divisible by C^p in the order, η being integral), and T ≡ 1. Correct.
* Simplicity: off C = 0 the point is unramified (Δ | C), t - α is a local parameter, and dv = 0 because v is a p-th
  power in the function field. So d(T - v) = T' dt, with T' = s C^p t^{s-1} or C^p, nonzero off C = 0 (t | C since k >= 2).
  Correct.
* **Independent proof that I prefer (and verified):** P_V is q-linearised, P_V = Σ a_i X^{q^i}, and each a_i is a
  polynomial in p-th powers. So d/dt kills the a_i and **R' = a_0·T'**, with a_0 = Π_{v≠0}(-v). Every prime factor of
  a_0T' lies over a zero of C, and R ≡ 1 mod C. So gcd(R, R') = 1. Numerically: gcd(R, C) = 1 in every run, T' | R'
  checked, and the full Euclidean gcd(R, R') = 1 was checked at L = 83 (deg R = 63341). Her `vii_cap2_check.py` uses the
  same identity (R' = c·e·A^{3u}).
* R ∈ F_p[t]: R is σ- and φ-invariant. In every run the η-part of R computed in F_3[t,η] is exactly 0.
* Cap two: N_0·Norm(G(T)) = R^2 **checked as an exact polynomial identity** in every run. With R squarefree, every
  prime exponent of N_0 is <= 2. "N_0 coprime to the discriminant" follows from R ≡ 1 mod C and Δ | C. Correct.
  **Remark:** the ramified primes therefore play **no role** in the construction. Ramification is present in the
  ring but is invisible to the mechanism.

## (d) Width at infinity

The pole order of a label is b = p(δ + deg C + k - 1). For L > b the bound 2 deg(z_i - z_j) <= (n+h-2g)L + b(n-h) holds
by the coefficient estimate. **I find it is attained with equality** in every computed case, so the threshold
L >= b(n-h)/(g-h) is exactly sharp for this family of T:

| model | m, q | g | predicted least L | computed: non-critical / critical |
|---|---|---|---|---|
| η^2 = 2t^2+1 | 2, q=3 | 5 | 45 | L=44 no, L=45 yes (D = 180) |
| η^2 = 2t^2+1 | 3, q=27 | 183 | 392 | L=391 no, L=392 yes (D = 214032) |
| η^2 = 2t^2+1 | 4, q=27 | 92 | 940 | L=930 no, L=940 yes (D = 598780) |

Padding (reserved primes of degrees 1, 2, 4, ... dividing C, so that R ≡ 1 mod each of them and e^2 N_0 keeps cap two) and the rate
computation (log threshold <= (2 log 2/log p + o(1)) m log m, giving liminf >= (1/2) log_2 p) were checked line by line.
They are identical to the referee-accepted VI.4 rate argument. Thresholds depend on m, and they are handled correctly: m is chosen as a
function of D.

## (e) Reproduction (own arithmetic, marks expanded exactly)

* **Her instance, ramified, η^2 = 2t^2+1, C = t^27 - t, T = 1 + tC^36 (L = 973):** g = 92, deg R = 709317,
  gcd(R,C) = 1, four **distinct** marks of degree **619801**, equal norms of degree 1239602, N_0·Norm(G(T)) = R^2 exactly,
  and actual pair-difference degrees {308815, 308773} <= 309900.5. Her certified upper bound 308815 is exactly the true
  maximum.
* **Same in the constant model (η^2 = -1):** pair degrees {307762, 307720}. This matches her 307762.
* q = 3, m = 2, η^2 = 2t^2+1, C = t^3 - t: everything checks for all s. The marks are critical exactly from L = 45, with D = 180,
  parent 2·180 = 360, and multiplicity profile {1: 180, 2: 90} (degrees of exponent-1 and exponent-2 parts, from my squarefree
  decomposition).
* The modular cross-check of all FFT products passed in every run.

## (f) Ramified tables (VII.3)

Definitions used, the same as hers: a point is z = A + Bη ∈ F_3[t,η], η^2 = 2t^2+1. Its pole order is D = max(deg A, 1 + deg B). A cluster
is a set of distinct points of pole order D with one norm N = A^2 - (2t^2+1)B^2 (any nonzero leading coefficient), such that each
pair difference has pole order <= floor(D/2). Unit multiples (±z) and conjugates count as distinct points. Non-primitive
points are allowed. "Squarefree" means N squarefree in F_3[t].

Method: "pole order of the difference <= H" is an equivalence relation (ultrametric), so M(D) = max over (high-coefficient class, N) of the number of
points. The sign flips A ↦ -A and B ↦ -B reduce the leading pair to (1,0), (0,1), (1,1). Run time: D = 9 in 16 s, D = 10 in 2.6 min.

| D | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| all parents (mine) | 2 | 3 | 2 | 4 | 3 | 4 | 3 | 4 | 4 | 4 |
| squarefree (mine) | 1 | 2 | 2 | 3 | 3 | 4 | 3 | 4 | 3 | 4 |

These are **identical to hers for all D = 1..10.** The independent clique brute force (all leading coefficients) agrees for
D <= 5. Example witnesses (N, high→low coefficients): D = 10 squarefree, four points, N = 100000001020111010020. D = 9 all
parents, four points, N = 2000000000201000200 (not squarefree). Her script is the same algorithm in numpy, so the agreement is a
reproduction, not a second method. The D <= 5 clique search is the genuinely different check.

Observation: the ramified values never exceed the constant-model values recorded for D <= 11 (M_3 <= 4), so small D
shows nothing special about ramification.

## (g) Does this close README_VI §5 item 2? Models not covered

* **Yes, for the question as posed.** README_VI §5 item 2 asks about F_3(t)(sqrt d), deg d even, nonsquare leading
  coefficient, ramified finite primes, and "whether any linearised-polynomial identity survives a non-constant
  conjugation". It does: the labels in η·F_q[t^p] are anti-invariant under the geometric σ, and P_V = A^q - c^{q-1}A is
  untouched. The model is unbounded, even with exponent cap two. The "ramification" and "conjugation not a Frobenius" items in the
  record's list of what a proof over Z may still use are now filtered, **at the cap-two level**. Positive
  characteristic is still essential (linearised polynomials, Frobenius descent from F_q to F_p, dv = 0), as Part VII
  itself says.
* **Not covered by Part VII as stated, but covered by my computation (and the proof transfers verbatim):** a
  **ramified infinite place** (deg Δ odd). The valuation at ∞ is half-integral (pole order max(deg a, deg b + deg Δ/2)),
  and every step (descent, anti-invariance, R' = a_0T', the coefficient estimate with b = p(deg Δ/2 + deg C + k - 1))
  is unchanged. Computed, q = 27, m = 3, forms x, y, x+y:
  * **Δ = t**, so O = F_3[u] with u^2 = t and σ: u ↦ -u. This is a **PID** with one principal ramified finite prime (u), a single
    place at infinity and units F_3^*: arguably the closest Z[i]-shaped model of all. Non-critical at L = 384,
    **critical at L = 385** (predicted 384.75). D = 210210, all checks pass.
  * Δ = t^3 - t (genus 1): non-critical at L = 398, **critical at L = 399** (predicted 398.25), D = 217854.
  The record (ref6_ternary.md line 220) had excluded A^2 - tB^2 as "leaving the model". It is now also unbounded at
  cap two.
* **Not covered:** p ≡ 1 mod 4 (χ(-1) = +1 breaks F_jσF_j = P_V^2; the record's half-system obstruction applies);
  infinity split (square leading coefficient: infinite unit group, not Z[i]-like, irrelevant); **squarefree parents**
  in any ramified model (open; the ramified squarefree table is small-D only); and p = 2 ramified (Artin–Schreier)
  models, which nobody has touched.

## Novelty and significance

* Novelty: real relative to the record. README_VI §§3, 5 list ramified models as untested, and nothing in the record
  constructs clusters in a geometric quadratic extension. The tables are new. The mathematics is a modest variant of
  VI.4: replace the label space by η^pC^p·F_q[t^p] and choose T so that the T - v separate. The new idea is the
  separation device (labels are p-th powers, so dv = 0, and collisions are confined to C = 0 where T = 1). It is clean and correct.
* Significance for the integer goal: **none directly.** It removes "ramification" and "conjugation is not a
  coefficient Frobenius" from the list of algebraic features that might power a proof over Z, but only in the
  presence of positive-characteristic additive structure. The residual Z-relevant question (squarefree or
  bounded-multiplicity N without linearised polynomials) is untouched. The cap-two clusters have doubled primes
  carrying a fraction 1 - O(1/m) of deg N. The growth rate (2/3+ε) log R/loglog R and M(C) over Z are unchanged.

## Corrected / added statements

1. VII.1 (ramified), as stated: **CORRECT.** Addendum: it holds verbatim for deg Δ odd (ramified infinity), e.g.
   O = F_3[u], u ↦ -u. Computed for m = 3; the proof transfers line by line.
2. Add to the model description: with inert infinity, |Cl(O)| = 2h(K) is even, odd-degree ramified primes are
   non-principal, and O is never a PID. The unbounded construction is insensitive to this (N_0 is prime to Δ).
3. The width bound in VII.1 is exact (attained) for T = 1 + C^p t^s and T = 1 + tC^36. The stated L-threshold is
   the true least L for these T.
