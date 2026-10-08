# A uniform count for mixed Fibonacci classes with a common slope

Mixed affine Fibonacci classes with a common positive slope have a
uniform finite orientation count, even with arbitrary odd cyclotomic
layers, arbitrary intercepts, and arbitrary layer multiplicities.
The proof uses the integral width budget to separate the classes at
large harmonics, then recovers the missing small harmonics modulo
three primes. This does **not** bound the full rational kernel of the
coupled equations; that kernel can have unbounded dimension.

This is a fixed-template theorem. It does not reduce arbitrary lattice
circles to Fibonacci templates. Constants below are deliberately crude.

## 1. Statement and normalization

Use the normalization of
[the mixed-class phase equations](mixed_affine_cyclotomic_rank_scope.md).
The classes are

\[
d_c(n)=An+B_c,
\qquad A>0\text{ even},\quad B_c\text{ distinct odd integers},
\]

on a fixed parity subsequence of \(n\). Negative intercepts are allowed;
all indices are eventually positive. Let \(a_{i,c,e}\) be the integer
orientation counts, where \(e\) ranges over finitely many odd positive
layer indices. Normalize by subtracting the minimum in each coordinate,
and let \(W_{c,e}>0\) be the actual remaining widths. Put

\[
D=\sum_{c,e}\phi(e)W_{c,e}>0,
\qquad \mathcal L=AD,
\qquad D_0=2^{32}.                                      \tag{1}
\]

Assume the rows have common total base parity and distinct phase
polynomials. A family with at most one row is trivial; otherwise let
\(Ah\) be its smallest pair contact frequency.
Here \(h\) is odd. Assume the endpoint budget

\[
D\le4h.                                                \tag{2}
\]

For actual fixed templates, equal-norm Gaussian prefactors and aligned
limiting phases give the total base parity condition; the endpoint
arc condition gives (2). No parity condition is imposed separately
on each class.

**Theorem.** If \(D\ge D_0\), the affine dimension of the row set is at
most \(1125\), and the number of rows is at most \(2251\). For every
\(D>0\), the number of rows is at most

\[
2D_0+1=8589934593.                                     \tag{3}
\]

These constants are independent of the common slope, intercepts,
layer orders, multiplicities, and number of classes.

## 2. Integral separation at large harmonics

Fix two rows and write

\[
b_{c,e}=a_{i,c,e}-a_{j,c,e},\qquad
S_{c,k}=\sum_e b_{c,e}c_e(k),
\]

where \(c_e(k)\) is the Ramanujan sum. Then

\[
\sum_c|S_{c,k}|\le\sum_{c,e}\phi(e)|b_{c,e}|\le D.       \tag{4}
\]

Let \(\varphi=(1+\sqrt5)/2\), \(q=-\varphi^{-2}\), choose one
intercept \(B_0\), and write \(B_c=B_0+2s_c\). The \(s_c\) are
distinct integers. The phase constants satisfy
\(\gamma_c/\gamma_0=q^{s_c}\), so vanishing at an odd harmonic
\(k<h\) is exactly

\[
P_k(q^k)=0,\qquad
P_k(X)=\sum_c S_{c,k}X^{s_c}\in\mathbb Z[X,X^{-1}].      \tag{5}
\]

The quadratic conjugate of \(q^k\) has modulus \(\varphi^{2k}\).
If a nonzero integer Laurent polynomial \(P\) vanishes there, multiply
it by a monomial and let its leading coefficient be \(a_m\ne0\).
At any root \(r\) with \(|r|>1\), the leading-term estimate gives

\[
|a_m||r|^m\le
(\|P\|_1-|a_m|)|r|^{m-1},
\quad\text{hence}\quad \|P\|_1\ge |r|+1.               \tag{6}
\]

Thus (4)--(6) force \(P_k=0\) whenever \(\varphi^{2k}>D\).
Distinct exponents then give \(S_{c,k}=0\) for every class. This
argument is independent of the span or signs of the intercepts.

Set

\[
K=\left\lceil\frac{\log D}{2\log\varphi}\right\rceil+1.
                                                               \tag{7}
\]

Every odd \(k\) with \(K\le k<h\) therefore vanishes class by class.

## 3. Recovering all smaller harmonics

For every prime \(p\) and every positive integer \(d\) with
\(p\nmid d\), the divisor formula gives

\[
c_e(dp)\equiv c_e(d)\pmod p.                            \tag{8}
\]

Indeed \(c_e(k)=\sum_{r\mid(e,k)}r\mu(e/r)\), and the new divisors
appearing when \(d\) is replaced by \(dp\) are multiples of \(p\).
This includes all possible depths of \(p\) in \(e\).

We verify one absolute threshold. Write \(x=\log D\). If
\(D\ge2^{32}\), then

\[
K\le3x,\qquad \sqrt D>192x.                             \tag{9}
\]

For the first inequality use \(2\log\varphi>\log2>1/2\) and
\(2x+2\le3x\). For the second, it holds at \(2^{32}\) since
\(\log(2^{32})<32\) and \(2^{16}>192\cdot32\); the function
\(\sqrt D/\log D\) increases thereafter.

Fix an odd \(d<K\), and put

\[
n=\left\lfloor\frac{h-1}{8d}\right\rfloor.
\]

Using (2) and (9),

\[
n\ge\frac{D}{32K}-2
 \ge\frac{D}{96x}-2
 >2\sqrt D-2\ge\sqrt D>K.                              \tag{10}
\]

Bertrand's postulate supplies three distinct primes

\[
p_1\in(n,2n),\quad p_2\in(2n,4n),\quad p_3\in(4n,8n).
\]

They are odd, exceed \(d\), and satisfy

\[
K<dp_j<h,\qquad p_1p_2p_3>8n^3>D.                       \tag{11}
\]

The large-harmonic conclusion gives \(S_{c,dp_j}=0\).
Equation (8) shows that their product divides \(S_{c,d}\), while
\(|S_{c,d}|\le D\). Therefore \(S_{c,d}=0\).

We have proved, for every actual row difference,

\[
S_{c,d}=0\quad\text{for every class }c
\text{ and every odd }d<h.                              \tag{12}
\]

These equations follow from bounded integral differences. They are
not asserted for arbitrary rational vectors in the coupled kernel.

**Upper-budget version.** The same argument applies with any integer
\(B\ge D_0\) in place of \(D\), provided
\(\sum_{c,e}\phi(e)|b_{c,e}|\le B\le4h\). In particular \(B\)
need not equal a width degree. Define \(K\) from \(B\), use
\(\|P_k\|_1\le B\), and repeat (9)--(11) with \(B\). The conclusion
is again (12). This form allows a common-slope frequency group to be
bounded using a larger global budget.

## 4. The colored rank bound with ordinary totient cost

Since

\[
S_{c,d}=\sum_{r\mid d}r B_{c,r},\qquad
B_{c,r}=\sum_{e:r\mid e}\mu(e/r)b_{c,e},
\]

triangular inversion changes (12) into the classwise truncated
Möbius equations \(B_{c,r}=0\) for odd \(r<h\).

Apply the prime-power compression from
[the uniform rank theorem](cyclotomic_uniform_rank.md) to each class
separately. Here use the ordinary cost \(\sum_e\phi(e)\), charging
the base index by \(1\), since its width need not be even in a single
class. The same proof applies: removing \(p^a s\) and inserting all
missing lower powers costs at most
\(\phi(s)p^{a-1}\), whereas the removed cost is
\(\phi(s)(p-1)p^{a-1}\). The nullity-preserving linear map is
unchanged. No additional base surcharge is needed.

It produces divisor downsets \(F_c\) with total cost

\[
U=\sum_c\sum_{e\in F_c}\phi(e)\le
\sum_{c,e:W_{c,e}>0}\phi(e)\le D\le4h.                  \tag{13}
\]

Their direct-sum nullity is exactly the number of colored high orders
\((c,n)\) with \(n\in F_c\), \(n\ge h\). View their cyclic root
subgroups \(H_n\) in disjoint copies indexed by \(c\). The union has
size \(U\). Same-color intersections have size \(\gcd(n,n')\);
different-color intersections are empty.

Join two high orders of the same color if
\(\gcd(n,n')\ge2h/15\). Each order is at most \(4h\). Writing
\(n=ga,n'=gb\) for their gcd \(g\), the positive odd integers
\(a,b\) are at most \(30\). For fixed \(n\), the pair \((a,b)\)
determines \(n'\), so there are at most \(15^2=225\) possibilities
including the vertex itself. The maximum degree is at most \(224\).

If there were at least \(1126\) high orders, a greedy selection would
give six independent vertices. Their root subgroups would have union
size, by the first two terms of inclusion-exclusion, strictly greater
than

\[
6h-\binom62\frac{2h}{15}=4h,
\]

contradicting (13). Thus the direct-sum kernel containing every row
difference has dimension at most \(1125\).

## 5. From affine dimension to count, and scope

The global sigma-reciprocity argument in
[the mixed-class note, Section 4](mixed_affine_cyclotomic_rank_scope.md)
gives the weighted pair separation

\[
\delta_{ij}=\sum_{c,e}\phi(e)|a_{i,c,e}-a_{j,c,e}|
\ge2h\ge D/2.                                         \tag{14}
\]

It uses common **total** base parity. Applying an even-degree
reciprocity argument to each class separately would be incorrect
when that class has odd base difference.

For completeness, embed each row by coordinates

\[
f_i(c,e)=\sqrt{\phi(e)W_{c,e}}
 \left(\frac{2a_{i,c,e}}{W_{c,e}}-1\right).
\]

The inequality \(uv\le1-|u-v|\) for \(u,v\in[-1,1]\) yields
\(\langle f_i,f_j\rangle\le D-2\delta_{ij}\le0\).
The affine obtuse-set bound gives at most \(2r+1\) distinct rows in
affine dimension \(r\); see
[the affine rank reduction](cyclotomic_affine_rank_reduction.md).
Section 4 gives \(r\le1125\) when \(D\ge D_0\). If \(D<D_0\),
the number of active coordinates is at most \(D\), so the same
argument gives at most \(2D+1\), proving (3).

The earlier fixed-shift base-only theorem gives the much sharper
bound of eight in that narrower setting. The present result permits
all odd layers, and thus raw factors with rates \(eA\), while the
class representatives retain the common slope \(A\). For actual
circle sequences, the index after which the normalized phase and
radius comparisons apply may depend on the fixed template. No
template-independent entry index or arbitrary-circle bound follows.

The independent exact arithmetic checks are in
[check_mixed_affine_prime_amplification.py](check_mixed_affine_prime_amplification.py).
They test the Ramanujan congruences and threshold arithmetic; the
proof above supplies the unbounded-parameter assertions.
