# A structural endpoint bound for fixed shifted-Pell products

This note proves a uniform count within a specified family class. It
does not prove the endpoint conjecture for arbitrary lattice circles.
In particular it rules out a counterexample to the intrinsic
eight-point bonus obtained from one fixed template of shifted Pell
factors, even allowing arbitrary fixed multiplicities and fixed
Gaussian prefactors.

Let
\[
\lambda=9+4\sqrt5,\quad F=1+2i,\quad
x_r+y_r\sqrt5=\lambda^r,\quad H_r=x_r+Fy_r.
\tag{1}
\]
Take finitely many distinct fixed integer shifts \(s\), nonnegative
fixed integers \(e_s\), and a fixed finite collection of rows
\[
Z_i(n)=c_i\prod_s
H_{n+s}^{\,a_{i,s}}\overline{H_{n+s}}^{\,e_s-a_{i,s}},
\quad 0\leq a_{i,s}\leq e_s,
\quad c_i\in\mathbf Z[i]\setminus\{0\}.
\tag{2}
\]
The \(c_i\) have one common norm. All shifts, exponents, row vectors,
and prefactors are independent of \(n\). We restrict to sufficiently
large \(n\) so that all \(n+s\) are positive. The rows have equal
norm for each \(n\).

**Theorem.** Suppose the collection (2) has at least two distinct
rows and, for an unbounded sequence of \(n\), all its primitive
normalized points lie on one arc of length at most \(C\sqrt R\),
with fixed \(C>0\). Suppose also that the primitive radii tend to
infinity. Put
\[
E_s=\max_i a_{i,s}-\min_i a_{i,s},\qquad K=\sum_s E_s.
\tag{3}
\]
Then \(K\leq4\), and there are at most eight distinct rows.
If there are eight rows, there are exactly four active factors, each
of width one, and the row vectors are one complete parity class of
the four-dimensional binary cube.

For such an eight-row endpoint family, its intrinsic conductor
statistics satisfy
\[
D=O(1),\qquad W_4=W-O(1),\qquad
D+W_4-2W=-W+O(1)\longrightarrow-\infty.
\tag{4}
\]
Constants here may depend on the fixed template. Thus this class
cannot produce a growing violation of the proposed eight-point
bonus at \(C=1/2\), or at any other fixed endpoint constant.

## 1. Primitive normalization does not hide a growing factor

First divide all rows by the exact common block factor
\[
\prod_s H_{n+s}^{\,\min_i a_{i,s}}
\overline{H_{n+s}}^{\,e_s-\max_i a_{i,s}}.
\]
Relabel the exponents so that every active coordinate has minimum
zero, maximum \(E_s\), and range \(0,\ldots,E_s\). The resulting
unnormalized radius is asymptotic to a positive constant times
\(\lambda^{Kn}\).

Distinct shifted blocks have uniformly bounded norm gcds. One direct
way to see this uses Fibonacci numbers \(f_0=0,f_1=1\):
\[
N(H_r)=f_{12r+1}.
\tag{5}
\]
For example, (5) follows immediately from Binet's formula and
\(\lambda=((1+\sqrt5)/2)^6\).
The elementary Euclidean identity
\(\gcd(f_u,f_v)=f_{\gcd(u,v)}\) follows from the Fibonacci addition
formula and \(\gcd(f_u,f_{u-1})=1\). Consequently, for \(s\ne t\),
\[
\gcd(N(H_{n+s}),N(H_{n+t}))
\mid f_{12|s-t|}.
\tag{6}
\]
Each individual \(H_r\) is conjugate-coprime: its coordinates are
relatively prime and its norm is odd.

Here is why (6) suffices even when an individual bad prime divides
one block to an unbounded exponent. Outside a fixed finite set of
primes, the block norm supports are disjoint and avoid the
prefactors. An active block occurs in both extreme orientations
among the rows, so such a prime cannot divide their common gcd.
At a prime in that finite set, all but at most one block have
bounded valuation, by (6). Choose a row using the opposite extreme
orientation of that possibly large block. Its valuation at the
given Gaussian prime is bounded by the fixed prefactor and the
remaining blocks. This bounds the common gcd valuation as well.
Fixed multiplicities only multiply these bounds by fixed constants.

Thus the remaining common Gaussian gcd has bounded norm. The
primitive radius satisfies
\[
R_n\asymp\lambda^{Kn}.
\tag{7}
\]
Removing this bounded gcd preserves all angular separations.

## 2. The exact angle expansion

Put
\[
\rho=\frac{\sqrt5-1}{2},\qquad q=\lambda^{-2},\qquad
h=\frac{1+F/\sqrt5}{2},\qquad \theta=\arg h.
\]
The exact Pell expression is
\[
H_{n+s}=h\lambda^{n+s}
 \bigl(1-i\rho q^{n+s}\bigr).
\tag{8}
\]
Therefore
\[
\arg H_{n+s}
=\theta-\arctan(\rho q^{n+s}).
\tag{9}
\]
Any two endpoint rows must have the same limiting argument modulo
\(2\pi\). Choose continuous argument lifts near that common limit.
For their exponent difference \(b_s=a_{i,s}-a_{j,s}\), let
\[
P(X)=\sum_s b_s X^{s-s_{\min}}\in\mathbf Z[X].
\]
It is nonzero: identical exponent vectors and equal limiting phases
would give identical fixed prefactors, hence identical rows.
Moreover
\[
\|P\|_1=\sum_s|b_s|\leq K.
\tag{10}
\]

After the limiting arguments cancel, (9) gives an odd-power series
whose coefficient at \(q^{\ell n}\), for odd \(\ell\geq1\), is
\[
-2\,\frac{(-1)^{(\ell-1)/2}\rho^\ell}{\ell}\,
q^{\ell s_{\min}}P(q^\ell).
\tag{11}
\]
There is a first odd \(\nu\) with \(P(q^\nu)\ne0\), because a
nonzero polynomial cannot vanish at the infinitely many distinct
numbers \(q,q^3,q^5,\ldots\). Their angular separation is therefore
asymptotic to a positive constant times \(\lambda^{-2\nu n}\).

The endpoint angle bound and (7) now imply
\[
K\leq4\nu.
\tag{12}
\]

## 3. Higher-order phase cancellation costs too many factors

If \(\nu\geq3\), then \(P(q^{\nu-2})=0\). Conjugating in
\(\mathbf Q(\sqrt5)\) shows that
\[
r=\lambda^{2(\nu-2)}>1
\]
is also a root of \(P\). Every nonzero integer polynomial with
a positive root \(r>1\) satisfies
\[
r\leq\|P\|_1-1.
\tag{13}
\]
Indeed, if its leading coefficient is \(b_d\), then
\[
|b_d|r^d
\leq(\|P\|_1-|b_d|)r^{d-1},
\]
and \(|b_d|\geq1\).

Equations (10), (12), and (13) would give
\[
\lambda^{2(\nu-2)}+1\leq K\leq4\nu,
\]
which is impossible for \(\nu\geq3\), since \(\lambda^2>300\).
Thus \(\nu=1\), and (12) proves
\[
K\leq4.
\tag{14}
\]
This argument handles arbitrary fixed multiplicities. It does not
assume uniqueness of binary subset sums.

## 4. Rational prefactors force a parity class

One has the exact identity
\[
\frac{h}{\overline h}=\frac{F}{\sqrt5}.
\]
For two rows with equal limiting argument, the ratio of fixed
prefactors consequently satisfies
\[
\frac{c_i}{c_j}
=\left(\frac{F}{\sqrt5}\right)^{-\sum_s(a_{i,s}-a_{j,s})}.
\tag{15}
\]
The left side lies in \(\mathbf Q(i)\). An odd power of
\(F/\sqrt5\) does not lie in \(\mathbf Q(i)\), whereas every even
power does. Hence all row exponent sums have the same parity.
For any one exponent vector, equality of the limiting phase fixes
its prefactor uniquely, because all prefactors have the same norm.

For \(K\leq3\) this gives at most four vectors. For \(K=4\), the
possible positive width partitions and their largest parity classes
are:

| Widths \(E_s\) | Largest parity class |
|---|---:|
| \(4\) | \(3\) |
| \(3,1\) | \(4\) |
| \(2,2\) | \(5\) |
| \(2,1,1\) | \(6\) |
| \(1,1,1,1\) | \(8\) |

Thus eight distinct rows force four binary active blocks and one
full parity class. Each block occurs in each orientation exactly
four times.

## 5. The intrinsic eight-point bonus in this class

On eight such rows, every growing block has a four-versus-four
orientation cut. At a prime where block supports are disjoint and
the prefactors have zero valuation, every threshold layer is
therefore balanced.

At a prime in the fixed exceptional set from Section 1, all but at
most one block have bounded valuation. Its eight valuation entries
have the form
\[
a_i=e\,v_i+\beta_i,\qquad
v_i\in\{0,1\},\qquad \sum_i v_i=4,
\]
where all \(\beta_i\) and all other conductor contributions are
bounded independently of \(n\). All but a bounded number of the
threshold layers are again four-versus-four, even if \(e\) grows.
The bounded common gcd removes only boundedly many layers.

Thus the total weight of unbalanced intrinsic layers is bounded.
Their squared defects are at most \(16\) times their weights, giving
\[
D=O(1),\qquad W-W_4=O(1).
\]
Since \(W=\log R_n^2\to\infty\), this proves (4). The conclusion
uses actual primitive conductor layers, not a formal layer model.

The template constants are essential. This result allows fixed
shifts and fixed multiplicities of any size, but does not cover
shifts, prefactors, or exponent ranges growing with \(n\). It also
does not turn its template-dependent bounds into a universal bonus
constant over all Gaussian configurations.

## 6. A nonresonant extension with different Pell rates

There is also a binary-factor extension. Replace the active factors
in (2) by
\[
H_{k_jn+s_j},\qquad k_j\in\mathbf Z_{>0},
\]
with fixed integer shifts and one binary orientation per factor.
Keep the fixed equal-norm prefactors. Suppose that every two active
factors satisfy
\[
k_i(12s_j+1)-k_j(12s_i+1)\ne0.
\tag{16}
\]
Then the norm gcds are again bounded: by (5), their Fibonacci
indices are \(12k_jn+12s_j+1\), and the gcd of two such indices
divides the fixed nonzero integer in (16).

After removing inactive factors, put
\[
S=\sum_j k_j,\qquad k_{\min}=\min_j k_j.
\]
The remaining primitive radius is comparable to \(\lambda^{Sn}\).
Some pair of rows differs on the group of factors having rate
\(k_{\min}\). Their leading angular coefficient is a nonzero
constant times
\[
\sum_{j:k_j=k_{\min}}b_jq^{s_j},
\qquad b_j\in\{-1,0,1\},
\]
with some \(b_j\ne0\). The shifts in this rate group are distinct.
After removing the lowest power, the coefficient of the constant
term has absolute value one, while all later powers have total
absolute value at most \(q/(1-q)<1\). Hence the sum is nonzero.
Higher-rate blocks cannot cancel this leading term.

The pair's angular separation is therefore comparable to
\(\lambda^{-2k_{\min}n}\). The endpoint condition forces
\[
S\leq4k_{\min}.
\tag{17}
\]
There are consequently at most four active binary factors.
The limiting-phase parity argument from Section 4 still applies,
and gives at most eight points. If eight occur, all four rates
must equal \(k_{\min}\), reducing to the balanced binary case of
Sections 4--5.

Condition (16) cannot simply be omitted from this normalization
argument. For example,
\[
N(H_n)=f_{12n+1},\qquad
N(H_{13n+1})=f_{13(12n+1)},
\]
so the first norm divides the second and their overlap grows.
These factors have different affine indices and would not be
removed by merely combining identical factors. Such resonant-rate
templates, as well as growing shifts or prefactors, remain outside
the proved class.

## 7. The finer Fibonacci Pell orbit has the same eight-row bound

The choice of step \(\lambda=9+4\sqrt5\) in (1) is not essential for
fixed equal-rate shifted templates. Put
\[
\varphi=\frac{1+\sqrt5}{2},\qquad
\lambda_0=\varphi^2=\frac{3+\sqrt5}{2},\qquad
J_r=f_{2r+1}+i f_{2r}\quad(r\geq0).
\tag{18}
\]
In (2), replace every \(H_{n+s}\) by \(J_{n+s}\), keeping exactly the
same hypotheses on fixed shifts, fixed multiplicities, equal-norm
Gaussian prefactors, endpoint arcs, and unbounded primitive radii.
Then the conclusion of the theorem still holds: \(K\leq4\) and at
most eight distinct rows; equality at eight forces four binary
active factors and one complete parity class. The intrinsic bonus
conclusion (4) also holds.

Here are the normalization and phase details, including the ramified
prime. Consecutive Fibonacci numbers are coprime, so any common
Gaussian divisor of \(J_r\) and \(\overline{J_r}\) divides \(2\): it
divides both \(2f_{2r+1}\) and \(2if_{2r}\), and an integer Bezout
identity finishes the claim. Thus their common divisor has bounded
norm even when the norm of \(J_r\) is even. Moreover
\[
N(J_r)=f_{4r+1},\qquad
\gcd\bigl(N(J_{n+s}),N(J_{n+t})\bigr)
 \mid f_{4|s-t|}\quad(s\ne t).
\tag{19}
\]
The same primewise argument as in Section 1 therefore removes only
a bounded common Gaussian gcd, and the primitive radius is
\(R_n\asymp\lambda_0^{Kn}\).

With \(\rho=\varphi^{-1}\), \(q_0=\lambda_0^{-2}\), and
\(h_0=(\varphi+i)/\sqrt5\), Binet's formula gives the exact identity
\[
J_{n+s}=h_0\lambda_0^{n+s}
       (1-i\rho q_0^{n+s}),\qquad
\arg J_{n+s}=\arg h_0-\arctan(\rho q_0^{n+s}),
\qquad \frac{h_0}{\overline{h_0}}=\frac{1+2i}{\sqrt5}.
\tag{20}
\]
For two distinct rows, the difference polynomial \(P\) from Section
2 has \(\|P\|_1\leq K\). If its first nonvanishing odd phase term is
\(\nu\), the endpoint bound still gives \(K\leq4\nu\). For
\(\nu\geq5\), conjugating the root \(q_0^{\nu-2}\) of \(P\) gives the
positive root \(\lambda_0^{2(\nu-2)}\); (13) would force
\[
\lambda_0^{2(\nu-2)}+1\leq K\leq4\nu,
\]
which is impossible for every odd \(\nu\geq5\).

The potentially resonant case is \(\nu=3\), for which \(K\leq12\)
and \(P(q_0)=0\). Since
\[
q_0^2-7q_0+1=0,
\]
Gauss's lemma writes \(P=(X^2-7X+1)Q\) with \(Q\in\mathbf Z[X]\).
The \(\ell_1\) convolution inequality gives
\[
\|P\|_1\geq(7-1-1)\|Q\|_1=5\|Q\|_1,
\]
so \(\|Q\|_1\leq2\). A single monomial of coefficient \(\pm2\)
costs \(18\). Two unit monomials cost at least \(14\): adjacent
equal signs give coefficients \((1,-6,-6,1)\), the minimum case;
opposite signs or wider separation give at least \(16\). Hence the
only possibility with \(\|P\|_1\leq12\) is
\(P=\pm X^a(X^2-7X+1)\), whose coefficient sum is \(\pm5\), odd.
But the last identity in (20) and the equal-norm Gaussian prefactors
force the exponent-sum difference \(P(1)\) to be even, exactly as in
Section 4. This contradiction rules out \(\nu=3\). Therefore every
distinct pair has \(\nu=1\), and the endpoint condition gives
\(K\leq4\), even for a template with only two rows. Sections 4--5
then give the asserted count and conductor conclusions.

Finally, \(J_{3r}=H_r\) exactly, so this result extends the original
fixed-rate Pell suborbit. It remains a bound for fixed templates, not
a radius-uniform count for arbitrary Gaussian configurations.
The [exact supplemental checker](check_finer_pell_parity_obstruction.py)
tests the Fibonacci recurrence, norm and shifted-gcd identities, and
the low-coefficient convolution cases.

## 8. A uniform count across fixed real-quadratic Pell units

The eight-row count is uniform over a wider fixed-template class.
Let \(D>0\) be a nonsquare integer, let \(F=a+ib\in\mathbf Z[i]\)
satisfy \(N(F)=a^2+b^2=D\), and let
\[
\lambda=x_1+y_1\sqrt D>1
\]
be a norm-one algebraic integer in \(\mathbf Q(\sqrt D)\), with
integer trace \(T=\lambda+\lambda^{-1}\geq3\). Write
\(\lambda^r=x_r+y_r\sqrt D\). Choose one fixed positive integer
\(c\) such that \(c x_r,c y_r\in\mathbf Z\) for every \(r\geq0\);
this is possible from the integer-trace recurrence. Put
\[
J_r=c(x_r+F y_r)\in\mathbf Z[i].                         \tag{21}
\]

Use these \(J_{n+s}\) in the fixed equal-rate template (2), with
distinct fixed shifts, fixed nonnegative multiplicities and exponent
vectors, and fixed Gaussian prefactors of one common norm. Assume
at least two distinct rows, the same fixed-endpoint arc hypothesis
along an unbounded sequence, and primitive radii tending to infinity.
Then, with \(E_s,K\) defined by (3),
\[
K\leq4,\qquad \#\{\text{distinct rows}\}\leq8.             \tag{22}
\]
The constant eight is independent of \(D,F,\lambda,c\), and the
fixed template. If eight rows occur, their active widths are four
ones and their exponent vectors form a complete parity class. The
intrinsic conclusion (4) follows with template-dependent constants.

Here is a bounded-gcd argument that does not rely on a special
Fibonacci identity. Put \(X_r=cx_r\), \(Y_r=cy_r\). They are integers
and \(X_r^2-DY_r^2=c^2\), so
\(\gcd(X_r,Y_r)\leq c\). The real and imaginary coordinates of
\(J_r\) are \(X_r+aY_r\) and \(bY_r\); their ordinary gcd is at most
\(|b|c\). As in Section 7, any common Gaussian divisor of \(J_r\)
and \(\overline{J_r}\) divides twice that ordinary gcd, so its norm
is bounded by \(4b^2c^2\). Notice that \(b\ne0\), since \(D\) is
nonsquare.

Set \(M_r=N(J_r)\) and \(A=T^2-2\). The Pell identity gives
\[
M_r=c^2(x_{2r}+a y_{2r}),\qquad
M_{r+2}=A M_{r+1}-M_r,\qquad
M_rM_{r+2}-M_{r+1}^2=c^4b^2y_2^2>0.                 \tag{23}
\]
The last quantity is a fixed nonzero integer, say \(\Delta\).
Consequently \(\gcd(M_r,M_{r+1})\leq\sqrt\Delta\). For a fixed
shift difference \(d>0\), the recurrence writes
\(M_{r+d}=U_dM_{r+1}-U_{d-1}M_r\), where the integers \(U_d\)
depend only on \(d,A\). Hence
\[
\gcd(M_r,M_{r+d})\leq |U_d|\sqrt\Delta.              \tag{24}
\]
This bounds every distinct shifted norm overlap. The primewise
extreme-orientation argument of Section 1, now including the bounded
within-block conjugate gcd, leaves only a bounded common Gaussian
factor after normalization. Thus the primitive radius is
\(R_n\asymp\lambda^{Kn}\).

The exact phase calculation is just as stable. Define
\[
h=\frac{1+F/\sqrt D}{2},\qquad
\rho=\frac{b}{\sqrt D+a}\ne0,\qquad q=\lambda^{-2}.
\]
Since \(N(F)=D\), direct division gives
\[
J_r=ch\lambda^r(1-i\rho q^r),\qquad
\arg J_r=\arg h-\arctan(\rho q^r),\qquad
\frac{h}{\overline h}=\frac{F}{\sqrt D}.           \tag{25}
\]
The ratio on the right has even powers in \(\mathbf Q(i)\) and odd
powers outside \(\mathbf Q(i)\). Therefore equal limiting phases
again force all row exponent sums to have the same parity.

For a distinct pair, let \(P\) and its first nonzero odd phase order
\(\nu\) be as in Section 2. Equations (23)--(25) give
\(\|P\|_1\leq K\leq4\nu\). If \(\nu\geq3\), the conjugate of the
earlier root \(q^{\nu-2}\) and (13) yield
\[
\lambda^{2(\nu-2)}+1\leq K\leq4\nu.              \tag{26}
\]
For \(T\geq4\), already \(\lambda^2+1>12\), and (26) is
impossible for every odd \(\nu\geq3\). The smallest trace is
\(T=3\), giving \(\lambda=(3+\sqrt5)/2\). In that case
\(\nu\geq5\) is excluded by (26), while \(\nu=3\) has
\(K\leq12\) and \(q^2-7q+1=0\). The exact convolution and parity
argument in Section 7 rules out this last case for every \(F\) of
norm \(D\): its only possible low-weight polynomial has odd
coefficient sum, contradicting (25). Thus every distinct pair has
\(\nu=1\), proving \(K\leq4\) even when there are only two rows.
Sections 4--5 finish (22) and the intrinsic assertion.

This corollary keeps the Pell unit, field, shifts, multiplicities,
and prefactors fixed as \(n\) grows. It does not cover varying units,
resonant unequal rates, or arbitrary same-norm lattice points.

## 9. Norm-minus-one units: a uniform fixed-template bound

Keep \(D,F,c\) and the fixed equal-rate shifted template of Section 8,
but let \(\lambda=x_1+y_1\sqrt D>1\) be a norm-minus-one algebraic
integer. Its integer trace \(T=\lambda-\lambda^{-1}\) is at least one,
so \(\lambda\geq\varphi\). Again choose the fixed positive integer
\(c\) to clear all \(x_r,y_r\), and put
\(J_r=c(x_r+F y_r)\). Under the same assumptions of at least two
distinct rows, fixed endpoint arc, and unbounded primitive radii,
the active width and number of rows satisfy
\[
K\leq12,\qquad \#\{\text{distinct rows}\}\leq2^{12}=4096. \tag{27}
\]
The bound is independent of the field, unit, and fixed template.
Together with Section 8, every template of this stated type built
from a fixed real-quadratic unit of norm \(+1\) or \(-1\) has a
radius-independent count. This is still a fixed-template result.

The normalization argument remains valid with the sign change.
Writing \(X_r=cx_r\), \(Y_r=cy_r\), one now has
\(X_r^2-DY_r^2=(-1)^r c^2\). Thus
\(\gcd(X_r,Y_r)\leq c\), and the Gaussian gcd of \(J_r\) with
\(\overline{J_r}\) is again bounded by \(2|b|c\). With
\(M_r=N(J_r)\), the exact formulas (23)--(24) hold with
\[
A=\lambda^2+\lambda^{-2}=T^2+2
\]
and the same positive Wronskian
\(\Delta=c^4b^2y_2^2\). Hence distinct fixed shifts have bounded
norm overlap, the common Gaussian normalization costs a bounded
factor, and \(R_n\asymp\lambda^{Kn}\).

The phase has a crucial alternating sign. With \(h,\rho\) as in
Section 8 and \(q=-\lambda^{-2}\),
\[
J_r=ch\lambda^r(1-i\rho q^r),\qquad
\arg J_r=\arg h-\arctan(\rho q^r),\qquad
\frac{h}{\overline h}=\frac{F}{\sqrt D}.              \tag{28}
\]
The last ratio still forces the exponent sums of rows sharing the
limiting angle to have the same parity. The endpoint comparison is
unchanged: for a nonzero pair-difference polynomial \(P\), if
\(\nu\) is its first nonvanishing odd phase order, then
\(\|P\|_1\leq K\leq4\nu\).

The conjugate of \(q^j\), for odd \(j\), is
\(-\lambda^{2j}\). If \(\nu\geq7\), the earlier root
\(q^{\nu-2}\) and the absolute-value form of (13) give
\[
\lambda^{2(\nu-2)}+1\leq K\leq4\nu,
\]
impossible already at \(\nu=7\), since
\(\varphi^{10}+1>28\), and increasingly so for larger odd \(\nu\).
The case \(\nu=5\) also fails by an elementary coefficient bound.
Put \(A=T^2+2\geq3\) and \(B=A^3-3A\geq18\). The distinct
minimal polynomials of \(q\) and \(q^3\) are
\[
f_1(X)=X^2+AX+1,\qquad f_3(X)=X^2+BX+1.
\]
Both divide \(P\), so \(P=f_1f_3Q\) in \(\mathbf Z[X]\).
For \(f_C=X^2+CX+1\), the reverse triangle inequality gives
\(\|f_C Q\|_1\geq(C-2)\|Q\|_1\). If \(A\geq4\), then
\[
\|P\|_1\geq(A-2)(B-2)>20=4\cdot5,
\]
contradicting the endpoint bound. If \(A=3\), then \(B=18\)
and \(\|P\|_1\leq20\) would imply
\(\|f_3Q\|_1\leq20\) and \(16\|Q\|_1\leq20\), hence
\(\|Q\|_1=1\). But then \(P\) is a signed monomial multiple
of \(f_1f_3\), whose coefficient norm is
\((3+2)(18+2)=100\), again a contradiction. Thus
\(\nu\leq3\) and \(K\leq12\). The number of integer exponent
vectors in the active box is at most
\(\prod_s(E_s+1)\leq2^{\sum_s E_s}=2^K\), proving (27).
As in Section 2, equal exponent vectors and equal limiting angles
force equal fixed prefactors, so each distinct row uses a distinct
vector.

Unlike the norm-one case, a third-order resonance can obey the
prefactor parity condition. For the smallest trace \(T=1\), the
minimal polynomial \(f_1=1+3X+X^2\) gives
\[
f_1(X)(1-X)=1+2X-2X^2-X^3,
\]
whose coefficient norm is six and whose coefficient sum is zero.
The width bound \(K=12\) is attained by an actual four-row template.
Take the unscaled Fibonacci factors
\(J_r=f_{r+1}+if_r\), shifts \(s=0,\ldots,4\), widths
\[
(E_0,E_1,E_2,E_3,E_4)=(1,3,4,3,1),
\]
all prefactors equal to one, and exponent vectors
\[
(0,0,2,3,1),\quad(0,1,4,1,0),\quad
(1,2,0,2,1),\quad(1,3,2,0,0).                    \tag{29}
\]
Every vector has sum six and polynomial value \(q^2\) at
\(q=-\varphi^{-2}\), using \(q^2+3q+1=0\). Every pair therefore
cancels the first phase term. Its difference cannot also cancel
the third: then \(f_1f_3\) would divide its difference polynomial
\(P\), forcing \(\|P\|_1\geq16>12\) by the convolution bounds
above. The primitive radius is
\(\asymp\varphi^{12n}\), while each pair angle is
\(\asymp\varphi^{-6n}\). These four actual points therefore lie on
an arc of length \(O(\sqrt R)\) for some fixed endpoint constant.
The cleared blocks of the theorem are twice these \(J_r\), a
fixed common scalar that does not change primitive normalized rows.
Their even-index suborbit is the norm-one family of Section 7.
The example attains the width bound; it does not approach the
4096-row count bound.

Only the global Gaussian gcd is bounded in this example. Individual
coordinate contents grow: the first and fourth displayed rows contain
the real factor `Norm(J_(n+2))^2`, while the second and third contain
`Norm(J_(n+1))Norm(J_(n+3))`. Individual coordinate primitivity is not
an assumption of the circle problem or of the fixed-template theorem.

The [exact supplemental checker](check_finer_pell_parity_obstruction.py)
also tests eight norm-\(\pm1\) unit fixtures, their recurrence and
Wronskian identities, the two minimal polynomials, and (29).
