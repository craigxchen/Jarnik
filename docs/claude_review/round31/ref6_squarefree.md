# Referee report, round 31, Part VI, block "squarefree" (VI.5 and the "sharpness for the profile in VI.2" paragraph)

Source: astra_free_VI/answer.md lines 1412-1441. Astra's script vi_squarefree_check.py was read but **not used**: all
computations below come from my own code, ref6_squarefree/sqf.py. It implements F2[t] as Python ints and F4[t] as pairs
(A,B) = A+B*theta, theta^2 = theta+1. Irreducibility uses a Rabin test, the absolute trace uses repeated squaring,
CRT uses Fermat inverses, and the residue of S at each prime is chosen at random among residues of the required
trace. The drivers are run2.py and `sqf.py 7`. Outputs are in run2.out and d7.out.

## Verdicts

| Claim | Verdict |
|---|---|
| VI.5 theorem: for m=2^(d-1), d>=3 and every L >= 4m(d-2), there is a squarefree N in F2[t] of degree 2(m-1)L carrying m points at c=0. Least degree 8m(m-1)(log2 m - 1); O(M^2 log M) for every M | **CORRECT** (checked line by line; reproduced for m = 4, 8, 16, 32, 64, even and odd L) |
| COMPUTED: CRT modulus degrees 6, 26, 60, 176; parent degrees 96, 896, 5760, 31744; 650 pairs | **CORRECT** (reproduced exactly) |
| Sharpness paragraph: the marks have precisely the complete Walsh-Hadamard profile (m-1 nonconstant columns) with pairwise disjoint block prime supports; with VI.2 this gives Theta(m^2 log m) within complete Hadamard cores | **CORRECT** for the identification and the disjointness. **Theta(m^2 log m) is conditional on VI.2**, which is refereed in another block. Its scope is a parent whose profile consists of exactly those m-1 blocks, after any common factor is removed. Constants: upper 4m^2 log2 m against lower m^2 log2 m / 6, a factor of about 24 |

## Line-by-line check of VI.5

1. **V and dim V = d.** f -> (f^2+f)^2 is F2-additive on U = span(t,...,t^(d-1)). Its kernel is U ∩ F2 = 0. Images of
   t^j have degree 4j, so the image has basis degrees 4,...,4(d-1), and every nonzero image has positive degree, so
   1 is not in the image. Hence dim V = d. Checked: |V| = 2^d distinct elements.
2. **Squares.** a+(f^2+f)^2 = (a+f^2+f)^2 because a is in F2. Checked: every v has derivative 0.
3. **Trace zero at even-degree primes.** Tr(g^2) = Tr(g), Tr(f^2+f) = 0, and Tr(1) = e mod 2 = 0. Checked at every
   prime dividing R, for all v.
4. **(VI.7), deg R <= 2m(d-2).**
   - R | H: a prime dividing 1+(f^2+f)^2 = (f^2+f+1)^2 divides a factor of H.
   - An even-degree prime dividing (f^2+f)^2 = f^2(f+1)^2, with f != 0 in U, divides f or f+1, so it has degree
     <= d-1. Then F2[t]_{<d-1} surjects onto the residue field. Multiplication by t is invertible because p != t (p
     has even degree), so U surjects too.
   - The residue field is F_{2^e} with e even, so it contains both roots of X^2+X+1. Two distinct f in U then give
     factors f^2+f+1 of H divisible by p.
   - With d >= 3 this applies to t^2+t+1, so t^2+t+1 divides H to multiplicity >= 2. Hence deg rad(H) <= deg H - 2.
   - deg H = 2 sum_{j=1}^{d-1} j 2^{j-1} = 2+2m(d-2), which I checked numerically. R | rad(H) gives (VI.7).
   - Computed deg R = 6, 26, 60, 176, 396 for d = 3..7, against bounds 8, 32, 96, 256, 640. H ≡ 0 mod R was checked.
     Primes of N-relevant differences have degree <= 2(d-1), so enumerating the primes up to that degree is complete.
5. **T = E + S^2 and the trace prescription.**
   - Tr_{F_{2^e}/F2}(theta) = (e/2) mod 2, the same for both embeddings because they differ by 1 and Tr(1) = 0.
   - Tr(S^2) = Tr(S) and the trace is onto F2, so S can be chosen by CRT with deg S < deg R <= L/2. Then
     2 deg S <= L-2, so deg T = L. Checked.
   - Astra's script writes T = t + S^2 with deg S = L/2. This is the same as t^L + t + S'^2 and is harmless.
6. **Pairwise coprimality.** Suppose alpha is a common root of T+theta+v and T+theta+v', with v != v'.
   - Then (v+v')(alpha) = 0 and theta = T(alpha)+v(alpha) lies in F2(alpha). So the minimal polynomial p of alpha has
     even degree and divides the nonzero element v+v' of V, hence p | R.
   - The trace of T(alpha)+theta is 1, while the trace of v(alpha) is 0: a contradiction.
   - This also covers a factor against its conjugate, since conj(T+theta+v) = T+theta+(v+1) and 1 is in V.
   - **Squarefreeness of each factor.** The derivative is T' (v is a square), which is 1 for even L and t^(L-1) for
     odd L. In the odd case the only possible repeated root is 0, and F(0) = theta + (element of F2) is nonzero.
   - Correct as written.
7. **Squarefree N.** N = prod_{v in V\{0,1}}(T+theta+v) is a product of pairwise coprime squarefree factors over F4,
   so it is squarefree over F4[t] and hence over F2[t]. One bookkeeping point: as an F2[t] product, N is the product
   of the m-1 block norms N_v = (T+v)^2+(T+v)+1 over 0 != v in E_0. It is **not** the product of the norms of all
   2m-2 factors, which would square it. My first draft made exactly this error and got a non-squarefree "N". The
   final code asserts that the F4 product equals (N,0).
8. **(V.1) still applies.**
   - The identity z_H - z_H' = P_J(1) P_J(Y)/Y uses only F2-additivity of P_J and u+u' ≡ 1 mod J, so it is valid for
     any F2-subspace V containing 1 and any Y.
   - J = H ∩ H' has dim d-2 and does not contain 1, so all nonzero elements of J have degree in {4,...,4(d-1)}.
     Echelon leading degrees therefore satisfy a_i <= 4(i+1), and deg(1+j) = deg j.
   - Sum: deg P_J(1) <= 4 sum_{i=1}^{d-2} 2^{i-1}(i+1) = 4(d-2)2^{d-2} = 2m(d-2) <= L/2.
   - **The bound is attained**: computed max deg P_J(1) = 8, 32, 96, 256, 640 = 2m(d-2), via J spanned by the images
     of t^2,...,t^(d-1). So within this V the threshold L >= 4m(d-2) is exact, and the least degree within the
     construction is 8m(m-1)(log2 m - 1).
   - Each pair has the exact gap deg(z_H - z_H') = (m/2-1)L + deg P_J(1). I checked this for every pair.
9. **Degrees.** There are 2m-2 factors of degree L, so deg N = 2(m-1)L and D = (m-1)L. At L = 4m(d-2) this gives
   deg N = 8m(m-1)(log2 m - 1).
10. **Arbitrary M.** Take a subset of the m = 2^r >= M marks. The norm and the pair bounds are inherited. The record's
    model (round29/ff.md §1) asks only for distinct points with one norm and the pair bound. The subset may share a
    common factor, which is allowed in the model, so no gcd cost is discarded. This gives O(M^2 log M).

## Reproduction (own code, sqf.py)

For each case the code builds every mark by multiplication in F4[t]. It checks:
- all m norms equal N;
- gcd(N, N') = 1;
- the m marks are distinct;
- each of the m(m-1)/2 pair differences has degree exactly (m/2-1)L + deg P_J(1), and twice that is <= D;
- H H^T = mI for the sign matrix with the all-ones column adjoined;
- the m-1 block norms N_v are pairwise coprime;
- each N_v is squarefree.

Each case was run with two random CRT seeds.

| m | d | L | deg R | #primes | deg N | D | max gap | pairs |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 3 | 16, 17, 18, 23 | 6 | 2 | 96 (at L=16) | 48 | 24 | 6 |
| 8 | 4 | 64, 65, 66, 71 | 26 | 6 | 896 | 448 | 224 | 28 |
| 16 | 5 | 192, 193, 194, 199 | 60 | 11 | 5760 | 2880 | 1440 | 120 |
| 32 | 6 | 512, 513, 514, 519 | 176 | 25 | 31744 | 15872 | 7936 | 496 |
| 64 | 7 | 1280 | 396 | 47 | 161280 | 80640 | 40320 | 2016 |

All checks pass. Odd D values such as L=17, D=51 give a max gap of 25, which is <= floor(D/2).

**Control (the trace condition matters).** With S random of degree < deg R and no trace prescription, N was
squarefree in 81/200 trials at m=4, 23/200 at m=8 and 5/200 at m=16. The CRT trace step is what makes the construction
work uniformly.

## Sharpness paragraph

- **Hadamard profile.** E_0 = span of the non-constant basis. Each conjugate pair {T+theta+v, T+theta+v+1}, for
  0 != v in E_0, is one block. Exactly one of v and v+1 lies in H because l(1) = 1. Restriction l -> l|_{E_0} is a
  bijection from the m forms to E_0^*, so the rows are the full character table of F2^(d-1) with the zero column
  removed. Checked: H H^T = mI.
- **Blocks.** Each f_v is monic of degree L and coprime to its conjugate. The block norms N_v are pairwise coprime
  (checked), so the prime supports are disjoint. No common factor remains, since each column is nonconstant and f_v,
  fbar_v are coprime. So the family satisfies the hypotheses of VI.2.
- **Bounds.** VI.2 gives D >= m(m-1)(log2(m-1)-3)/6, and the construction gives D = 4m(m-1)(log2 m - 1). That is
  Theta(m^2 log m) along powers of two, conditional on VI.2.
- **Spot check of VI.2's intermediate inequalities on this family.**
  - depth(r_v) = L >= D/m, because f_v + fbar_v = 1.
  - deg(A_vB_w + A_wB_v) = deg(v+w) <= 4(d-1) <= 2L/m = 8(d-2) for d >= 3, which is what depth(r_v/r_w) >= 2D/m
    requires. This is consistent.
- **Scope.** Theta(m^2 log m) is a statement about the complete Hadamard prime profile (one construction class). It is
  not a statement about the model.

## Novelty

- The record lists "is sup_D M_2(D,0) still infinite for squarefree N?" as open (round31/README_V.md item 3 and the
  caveats around lines 76, 278 and 387). It also records that the Part V N is never squarefree (deg gcd(N,N') = 8 and
  36).
- VI.5 answers this open question: **yes, infinite**, with explicit parents of degree 8m(m-1)(log2 m - 1), which is
  4 times the non-squarefree Part V degree.
- The device is new to the record: replace V by a space of squares with trace zero at all even-degree primes, and
  choose T by a CRT trace condition at the split primes dividing differences.

## Significance (blunt)

- **Modest.** It closes question (c) and extends the filter F3+ to squarefree parents. No argument valid verbatim over
  F2[t] with this norm, even one that uses squarefree N essentially, can prove a uniform bound.
- **Growth rate.** It says nothing about the growth rate. The degree is Theta(m^2 log m), about (8/3)m times the
  growth scale 3M log2 M.
- **Integers.** It does not transfer to Z[i]: it is a characteristic-2 Artin-Schreier translation mechanism, killed
  by Lemma P in odd characteristic.
- **Sharpness pairing.** It shows VI.2 is the right order for complete Hadamard cores, which is a statement about one
  construction class.

## What remains open

1. **Squarefree at growth scale.** Is the least deg N of a squarefree binary parent carrying M points at c=0 of order
   M log M, M^{1+o(1)}, or M^2 log M? Known: 3M log2 M - O(M) (the record's growth bound, valid for all N)
   <= W_sqf(M) <= 8m(m-1)(log2 m - 1) with m = 2^ceil(log2 M). The same gap is open without squarefreeness, with
   upper bound 2m(m-1)(log2 m - 1).
2. **Least squarefree degree for fixed m.** The exact least squarefree degree is unknown even for m = 4: the
   construction gives deg N = 96, against the non-squarefree minimum D = 12 (deg N = 24). The exhaustive M_2(D,0)
   tables in the record do not record whether their extremal N are squarefree.
3. **Within this V the threshold is exact**, because max deg P_J(1) = 2m(d-2) is attained. Any improvement of the
   factor 4 over Part V needs a different V. The squaring is used both for the derivative argument and for the trace
   condition.
4. **Theta(m^2 log m) for Hadamard cores** depends on VI.2 being correct (another block).

Scripts: /tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref6_squarefree/{sqf.py, run2.py, run2.out, d7.out}
