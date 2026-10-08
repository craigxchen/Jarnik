# A uniform cardinality bound for fixed-shift quadratic-unit templates

## The theorem and its scope

Let \(K\) be a real quadratic field, let \(\tau\) be its nontrivial
automorphism, and extend \(\tau\) to \(L=K(i)\) by fixing \(i\).
Let \(\lambda>1\) be a quadratic algebraic-integer unit, with

\[
\sigma=N_{K/\mathbf Q}(\lambda)\in\{1,-1\},\qquad
\lambda'=\tau(\lambda)=\sigma\lambda^{-1}.
\]

Fix \(a\in L^\times\) with

\[
\boxed{\operatorname{Im}\bigl(\tau(a)/a\bigr)\ne0.}
\tag{1}
\]

For integers \(r\), put

\[
H_r=a\lambda^r+\tau(a)(\lambda')^r.
\tag{2}
\]

These values lie in \(\mathbf Q(i)\). Their integrality in
\(\mathbf Z[i]\) is assumed on the indices used below.
Fix finitely many distinct integer shifts \(s\), integers \(e_s\geq0\),
and finitely many rows

\[
z_j(n)=c_j\prod_s H_{n+s}^{\,a_{j,s}}
                    \overline{H_{n+s}}^{\,e_s-a_{j,s}},
\quad
0\leq a_{j,s}\leq e_s,
\tag{3}
\]

where \(c_j\in\mathbf Z[i]\setminus\{0\}\) have one common modulus.
All shifts, exponents, and prefactors are fixed as \(n\) varies.
The rows must be distinct for the indices in question.

Define the active widths and their sum by

\[
w_s=\max_j a_{j,s}-\min_j a_{j,s},\qquad
K_0=\sum_s w_s.
\tag{4}
\]

Divide the evaluated rows by their common Gaussian gcd, and let \(R_n\)
be the resulting primitive radius. Suppose \(R_n\to\infty\) along an
unbounded sequence of indices, and along that sequence the primitive
points occupy an arc of length at most \(C\sqrt{R_n}\), for one fixed
\(C>0\).

Then

\[
\boxed{K_0\leq12,\qquad \#\{z_j\}\leq8.}
\tag{5}
\]

Both bounds are independent of the quadratic field, the unit, the shifts,
the exponents, and the fixed Gaussian prefactors.
Seven or more rows force \(K_0\leq4\); eight rows force four binary
active coordinates and one full parity class. Sections 6--7 supply
these sharper cardinality conclusions after the initial width bound.
If \(\lambda>\sqrt{11}\), the stronger conclusions are

\[
\boxed{K_0\leq4,\qquad \#\{z_j\}\leq8.}
\tag{6}
\]

In particular (6) applies to \(\lambda=9+4\sqrt5\).
No pairwise coprimality hypothesis on the evaluated blocks is needed:
the bounded common-factor loss used in the proof is established below.

For `C < 2 sqrt(2)`, the cardinality bound improves to **six**.
The [multiplicative-rectangle corollary](multiplicative_rectangle_separation.md)
combines Section 7.4 below with exact cancellation of fixed prefactors
in complementary exponent pairs. The bounds of eight stated above
remain valid for arbitrary fixed `C`.

When `lambda>sqrt(11)` and `10 C^2<=8`, the count further improves
to **five**. The [small-width circuit argument](small_width_affine_circuits.md)
combines the width-four conclusion with product rigidity through ten
factors. In particular this applies at `C=1/2` to `9+4 sqrt(5)`.

This is a theorem about fixed templates whose endpoint property persists
along unbounded indices. It gives no uniform radius threshold beyond which
an arbitrary isolated configuration must belong to such a template, or
beyond which a varying template has entered its asymptotic regime.
It does not prove the uniform endpoint theorem for arbitrary circles.

## 1. Fixed shifts have bounded common factors

For two different oriented blocks from the finite list

\[
H_{n+s},\quad \overline{H_{n+s}},
\]

their common Gaussian gcd has bounded norm, independently of \(n\).
Here the two orientations for one shift count as different blocks.

To prove this, write the two functions as

\[
f(n)=A\lambda^n+B(\lambda')^n,\qquad
g(n)=D\lambda^n+E(\lambda')^n
\]

with fixed coefficients in \(L\). If \(AE-BD\ne0\), a common Gaussian
divisor \(\gamma\) of the integral values \(f(n),g(n)\) divides a fixed
nonzero algebraic integer in \(L\), up to multiplication by a unit.
Indeed, clear the coefficient denominators by one fixed rational integer
\(m\), and eliminate one of \(\lambda^n,(\lambda')^n\) from
\(mf(n),mg(n)\). This gives divisibility in the ring of integers of \(L\)
by

\[
m^2(AE-BD)\lambda^n.
\]

The power \(\lambda^n\) is a unit. Taking absolute field norms bounds
\(N_{\mathbf Q(i)/\mathbf Q}(\gamma)^2\) by a fixed positive integer.

The required determinants do not vanish. For the same orientation and
different shifts, their nonzero scalar factor is

\[
\lambda^s(\lambda')^t-\lambda^t(\lambda')^s,
\]

which is nonzero because \(|\lambda/\lambda'|>1\).
For \(H_{n+s}\) and \(\overline{H_{n+t}}\), vanishing would imply

\[
\frac{a\,\tau(\overline a)}
     {\overline a\,\tau(a)}
=\left(\frac{\lambda}{\lambda'}\right)^{t-s}.
\tag{7}
\]

The left side has complex modulus one. Thus \(t=s\), and then (7)
would make \(\tau(a)/a\) real, contrary to (1).
The remaining conjugated cases follow by conjugation.

Because there are only finitely many oriented pairs, there is one fixed
nonzero Gaussian integer \(D_0\) divisible by every such pairwise gcd for
all permitted indices. One can obtain \(D_0\) by taking a common multiple
of the finitely many possible gcds of bounded norm; no uniform bound on
\(D_0\) across different templates is asserted.

Now divide every row in (3) by the exact common product

\[
G_{\rm obvious}(n)=
\prod_s H_{n+s}^{\,\min_j a_{j,s}}
 \overline{H_{n+s}}^{\,e_s-\max_j a_{j,s}}.
\]

The remaining rows have exponents from \(0\) to \(w_s\), and each
orientation is absent from at least one row.
Their remaining common gcd divides

\[
\left(\prod_j c_j\right)D_0^{K_0}
\tag{8}
\]

up to a unit. To check (8) at a Gaussian prime, select an oriented block
with maximal valuation and a row in which that block is absent.
Every other block in that row has valuation at most the valuation of
\(D_0\), since its gcd with the maximal block divides \(D_0\).
The total block multiplicity in that row is \(K_0\).

Finally

\[
|H_{n+s}|\sim |a|\lambda^{n+s}.
\]

Equation (8) therefore proves

\[
R_n\asymp\lambda^{K_0n}.
\tag{9}
\]

The constants may depend on the fixed template. In particular, the
assumed unbounded primitive radius implies \(K_0\geq1\).

## 2. A phase expansion that covers both unit norms

If \(\sigma=-1\), restrict the unbounded index sequence to one parity.
This leaves an unbounded subsequence. Define

\[
q=\sigma\lambda^{-2},\qquad
\eta=\sigma^n\,\frac{\tau(a)}a,\qquad
x=\lambda^{-2n}.
\tag{10}
\]

On the selected parity, \(\eta\) is fixed and nonreal. In the norm-one
case it is fixed without passing to a subsequence.
Exactly,

\[
H_{n+s}=a\lambda^{n+s}(1+\eta q^s x).
\tag{11}
\]

Thus negative shifts cause no problem: \(q^s\) is a fixed nonzero real
number, and the logarithmic expansion is valid for sufficiently small
positive \(x\).

For distinct rows \(j,k\), put

\[
\delta_s=a_{j,s}-a_{k,s},\qquad
P(X)=\sum_s\delta_sX^s\in\mathbf Z[X,X^{-1}].
\tag{12}
\]

The widths give

\[
\|P\|_1:=\sum_s|\delta_s|\leq K_0.
\tag{13}
\]

If \(P=0\), the row quotient is the fixed number \(c_j/c_k\).
The endpoint angular width tends to zero by (9), so this quotient must
equal one; the rows would be identical. Thus \(P\ne0\).

The limiting phase of their quotient must likewise equal one. After
canceling this limit, its analytic logarithm at \(x=0\) is

\[
L_P(x)=
\sum_{\ell\geq1}
\frac{(-1)^{\ell+1}}{\ell}
  \bigl(\eta^\ell-\overline\eta^{\,\ell}\bigr)
  P(q^\ell)x^\ell.
\tag{14}
\]

All coefficients are purely imaginary. Let \(\nu\) be the first order
with a nonzero coefficient.
This order exists: a nonzero Laurent polynomial has finitely many
nonzero roots, whereas the \(q^\ell\) are infinitely many distinct
numbers; and the nonreal number \(\eta\) has infinitely many powers with
nonzero imaginary part.

The first nonzero coefficient in (14) gives a nonzero phase difference
of order \(\lambda^{-2\nu n}\). The endpoint assumption and (9) give an
upper bound of order \(\lambda^{-K_0n/2}\). Hence

\[
\boxed{K_0\leq4\nu.}
\tag{15}
\]

The conclusion uses only a comparison of exponents. Its starting index
and its implicit constants can depend on the fixed coefficients of
\(P\), the shifts, and the prefactors.

## 3. Two vanished moments force too much coefficient mass

For a nonzero integer polynomial
\(Q(X)=b\prod_\ell(X-r_\ell)\), its Mahler measure is

\[
M(Q)=|b|\prod_\ell\max(1,|r_\ell|).
\]

We use only the elementary inequality

\[
M(Q)\leq\|Q\|_1.
\tag{16}
\]

For completeness, the mean of \(\log|e^{it}-r|\) over the unit circle
is \(\log\max(1,|r|)\). Applying this identity to the linear factors
expresses \(\log M(Q)\) as the mean of \(\log|Q(e^{it})|\), which is at
most \(\log\|Q\|_1\). Roots on the unit circle follow by a limit.
For a Laurent polynomial, multiply by a monomial; both quantities in
(16) remain unchanged.

The conjugate of \(q\) in \(K\) is

\[
\tau(q)=\sigma\lambda^2.
\tag{17}
\]

If \(P(q^\ell)=0\), its integer coefficients also give
\(P(\tau(q)^\ell)=0\). This conjugate root has modulus
\(\lambda^{2\ell}>1\).
For different positive orders, these roots are distinct.
Moreover \(q^\ell\) is irrational: it is a rational number of field norm
one only if it is \(1\) or \(-1\), excluded by \(|q^\ell|<1\).
Thus different orders correspond to distinct quadratic factors.

Consequently if two distinct positive orders \(u,v\) vanish, then

\[
\|P\|_1\geq M(P)\geq\lambda^{2(u+v)}.
\tag{18}
\]

Next, the numbers
\(\eta^\ell-\overline\eta^{\,\ell}\) cannot vanish at two consecutive
positive orders. If they did, the two corresponding powers of
\(\eta/\overline\eta\) would both equal one, forcing
\(\eta=\overline\eta\).
The coefficient at order one is nonzero by (1).

Suppose now that \(\nu\geq4\). Since it is the first nonzero order in
(14), we have \(P(q)=0\). Among the two orders \(\nu-2,\nu-1\), at
least one, say \(j\), has
\(\eta^j-\overline\eta^{\,j}\ne0\); hence \(P(q^j)=0\).
Here \(j\geq\nu-2\geq2\), so (18) gives

\[
\|P\|_1\geq\lambda^{2(1+j)}
             \geq\lambda^{2(\nu-1)}.
\tag{19}
\]

Every real quadratic algebraic-integer unit \(\lambda>1\) satisfies
\(\lambda\geq\varphi=(1+\sqrt5)/2\).
Indeed, in norm minus one its positive integer trace is
\(\lambda-\lambda^{-1}\geq1\); in norm one its integer trace is
\(\lambda+\lambda^{-1}\geq3\).
Also

\[
\varphi^{\,2(\nu-1)}>4\nu\qquad(\nu\geq4).
\tag{20}
\]

The base case is \(\varphi^6=8\varphi+5>16\), and induction follows by
multiplying by \(\varphi^2>2\).

Equations (13), (15), (19), and (20) contradict one another. Therefore

\[
\boxed{\nu\leq3,\qquad K_0\leq12.}
\tag{21}
\]

If \(\nu=3\), then \(\eta\) must be purely imaginary. Otherwise the
coefficient at order two is nonzero, so \(P(q)=P(q^2)=0\), and (18)
would give \(\|P\|_1\geq\varphi^6>12\), again contradicting (13)--(15).
Thus a nonreal phase ratio that is not purely imaginary gives the
additional bound \(K_0\leq8\).

## 4. The sharper bound for larger units

A nonzero integer polynomial with a root \(r\) satisfying \(|r|>1\)
obeys

\[
|r|\leq\|P\|_1-1.
\tag{22}
\]

In fact the leading-term inequality gives the stronger bound
\(|r|\leq\|P\|_1/|b|-1\), where \(b\) is its leading coefficient.
This also applies to a Laurent polynomial after multiplying by a
monomial.

If \(\nu=2\) or \(3\), then \(P(q)=0\), and the conjugate root has
modulus \(\lambda^2\). Equations (13), (15), and (22) imply

\[
\lambda^2+1\leq K_0\leq4\nu\leq12.
\]

Thus \(\lambda>\sqrt{11}\) forces \(\nu=1\), and hence \(K_0\leq4\).

## 5. Counting the rows

After subtracting the minimum exponent in each coordinate, there are
at most \(\prod_s(w_s+1)\leq2^{K_0}\) possible exponent vectors.
One further factor of two follows from the leading phase.

Put \(u=a/\overline a\). Condition (1) implies
\(u\notin\mathbf Q(i)\); otherwise \(\tau(u)=u\) would make
\(\tau(a)/a\) real.
Equality of the limiting phases of two rows implies

\[
\frac{c_j}{c_k}\,u^{\sum_s(a_{j,s}-a_{k,s})}=1.
\tag{23}
\]

Their total exponent sums therefore cannot differ by \(1\) or \(-1\).

Encode an exponent \(b_s\in\{0,\ldots,w_s\}\) by a binary string of
length \(w_s\), consisting of \(b_s\) ones followed by zeros.
Concatenating these strings injects the row vectors into the binary
cube of dimension \(K_0\). Two adjacent cube vertices have total
exponent sums differing by one, so they cannot both represent rows.
Pair all cube vertices by flipping one fixed coordinate. Each pair
contains at most one row, giving

\[
\#\{z_j\}\leq2^{K_0-1}.
\tag{24}
\]

Together with (21), this initially gives at most \(2048\) rows.
It also proves (6). Sections 6--7 improve the unrestricted-unit
cardinality bound to eight.

## 6. Eight rows reduce to the smallest norm-minus-one unit

Suppose there are eight rows and \(K_0>4\). Then every pair has
\(\nu\geq2\), so every pair polynomial satisfies

\[
P(q)=0.
\tag{25}
\]

We first check that their coefficient lengths are even. Write
\(K=\mathbf Q(\sqrt D)\), with \(D>1\) squarefree.
If some nonzero power of \(u=a/\overline a\) lies in \(\mathbf Q(i)\),
then \(\zeta=\tau(u)/u\) is a root of unity satisfying
\(\tau(\zeta)=\overline\zeta=\zeta^{-1}\).
Consequently \(\zeta=x+iy\sqrt D\), with \(x,y\in\mathbf Q\).
The rational algebraic integer \(2x=\zeta+\zeta^{-1}\) belongs to
\(\{-2,-1,0,1,2\}\), and \(x^2+Dy^2=1\).
The option \(x=0\) is impossible for nonsquare \(D\); the options
\(x=\pm1/2\) require \(D=3\).
For \(D\ne3\), therefore, \(\zeta=\pm1\).
The positive sign is excluded by (1), leaving \(\zeta=-1\).
Equation (23) then forces all total exponent sums to have the same parity.
If no nonzero Gaussian-rational power of \(u\) exists, (23) instead
forces their total exponent sums to be equal, which also suffices.
Thus, for \(D\ne3\),

\[
\|P\|_1\equiv P(1)\equiv0\pmod2
\quad\hbox{for every row pair.}
\tag{26}
\]

The exceptional field \(D=3\) does not arise when \(K_0>4\).
Its real quadratic units have norm one and smallest value greater than
one equal to \(2+\sqrt3>\sqrt{11}\); Section 4 already gives \(K_0\leq4\).
For example, the norm-minus-one equation in
\(\mathbf Z[\sqrt3]\) is impossible modulo \(3\), and the
norm-one equation has least positive first coordinate \(2\).

By (22) and (25), every pair satisfies

\[
\|P\|_1\geq\lambda^2+1.
\tag{27}
\]

For norm one, \(\lambda\geq\varphi^2\), so the integer on the left
is at least \(8\). For norm minus one with trace at least \(2\),
\(\lambda\geq1+\sqrt2\), and (27) first gives \(\|P\|_1\geq7\);
the parity in (26) improves this to \(8\).

On the other hand, summing the coefficient distances over all 28 row
pairs gives

\[
\sum_{j<k}\sum_s|a_{j,s}-a_{k,s}|\leq16K_0\leq192.
\tag{28}
\]

Indeed, each of the \(w_s\) integer threshold layers in one coordinate
separates at most \(4\cdot4=16\) pairs.
The lower bound \(28\cdot8=224\) contradicts (28).

We have proved the following stronger reduction:

\[
\boxed{\text{Eight rows with }K_0>4
\ \Longrightarrow\
\sigma=-1,\quad\lambda=\varphi,\quad
q=-\varphi^{-2}.}
\tag{29}
\]

The last equality for the unit follows because its remaining possible
positive integer trace is \(1\).

### The exact remaining finite-width question

Here \(q\) has minimal polynomial

\[
M(X)=X^2+3X+1.
\tag{30}
\]

The active row polynomials have identical evaluations at \(q\), and
their exponent sums have the same parity. Equivalently, each pair
difference is \(P=M Q\) for an integer Laurent polynomial \(Q\), with
\(P(1)\) even.
The triangle inequality gives the useful coercive estimate

\[
\|P\|_1=\|X^2Q+3XQ+Q\|_1\geq\|Q\|_1.
\tag{31}
\]

Every nonzero such pair difference has \(\|P\|_1\geq6\).
To prove this, lengths below six, being even, would be two or four.
Expand a putative relation into that many signed monomials
\(\epsilon q^s\). Each is a nonzero real number of quadratic norm one.
For four terms, group the relation as \(A+B=C+D\).
Conjugation gives the corresponding equality of reciprocal sums.
If the common sum is nonzero, equality of reciprocal sums implies
\(AB=CD\), so the two unordered pairs coincide.
If the common sum is zero, each pair cancels internally.
In either case there is cancellation of identical powers with opposite
signs, since distinct powers of \(q\) have distinct absolute values.
That contradicts having used the reduced coefficient list of a nonzero
polynomial. The two-term case is immediate.

This further shows that eight rows require

\[
168=28\cdot6\leq16K_0,\qquad K_0\in\{11,12\}.
\tag{32}
\]

Their phase ratio \(\eta\) must also be purely imaginary.
Otherwise Section 3 gives \(K_0\leq8\), contradicting (32).
Thus every distinct row pair has first nonzero phase order exactly
\(\nu=3\).

The remaining finite-width question concerns integer vectors in arbitrary
finitely supported boxes with the same \(q\)-evaluation and exponent-sum
parity. The next section resolves the needed bound by an exact
classification of the shortest nonzero difference. It does not use an
enumeration restricted to consecutive shifts.

## 7. The golden-ratio evaluation fiber has at most six rows

**Finite-width theorem.** Let finitely supported integer vectors
\(b_j=(b_{j,s})_{s\in\mathbf Z}\) lie in boxes of total width at most
\(12\). Suppose their evaluations
\(\sum_s b_{j,s}q^s\) are equal, where \(q=-\varphi^{-2}\), and their
coefficient sums have the same parity. Then there are at most six
distinct vectors.

### 7.1. Classification of the shortest trade

Every pair difference is divisible by \(M=X^2+3X+1\), and the preceding
argument proves its nonzero coefficient length is at least six.
The exact equality case is

\[
\boxed{\|P\|_1=6,\quad M\mid P
\quad\Longrightarrow\quad
P=\pm X^s T,\qquad
T=M(1-X)=1+2X-2X^2-X^3.}
\tag{33}
\]

Here \(P\) is an integer Laurent polynomial, and \(s\) is an integer.

To prove (33), multiply by a monomial so that
\(P=p_0+\cdots+p_dX^d\) has \(p_0p_d\ne0\).
Put \(\beta=\varphi^2\). Its root \(-\beta\) and the leading-coefficient
version of (22) give

\[
|p_d|\leq\frac6{\beta+1}<2.
\]

Applying the same argument to the reciprocal polynomial gives
\(|p_0|=|p_d|=1\).

The neighboring coefficient satisfies \(|p_{d-1}|\geq2\).
If it is zero, the leading term at \(-\beta\) would give
\(\beta^2\leq5\), whereas \(\beta^2>5\).
If its absolute value is one, the first two terms have combined
absolute value at least
\(\beta^{d-1}(\beta-1)\). The remaining coefficient mass is four,
which would imply \(\beta(\beta-1)\leq4\), again false.
The reciprocal argument gives \(|p_1|\geq2\).

Degree two is impossible: a degree-two multiple of \(M\) with leading
coefficient \(\pm1\) is \(\pm M\), of length five.
For \(d\geq3\), the four endpoint positions are distinct and use all six
units of coefficient mass. Thus the interior coefficients vanish, and
the endpoint magnitudes are \(1,2,2,1\).
The coefficients \(p_{d-1}\) and \(p_d\) have the same sign.
Otherwise their contribution at \(-\beta\) has magnitude
\(\beta^{d-1}(\beta+2)\), which cannot be canceled by the remaining
mass three in lower degrees. Reciprocally, \(p_0,p_1\) have the same sign.
Consequently

\[
P=\epsilon(1+2X)+\delta X^{d-1}(2+X),
\qquad \epsilon,\delta\in\{-1,1\}.
\]

Both \(1+2q\) and \(2+q\) are positive, and \(M(q)=0\) implies

\[
\frac{1+2q}{2+q}=q^2.
\]

Taking absolute values in \(P(q)=0\) gives
\(|q|^{d-1}=|q|^2\), so \(d=3\); then \(\delta=-\epsilon\).
This proves (33).

### 7.2. Shortest-trade edges form an integer-grid graph

Make a graph on the row vectors by joining precisely the pairs at
coefficient distance six. In each connected component, choose a base
row. By (33), every other row differs from it by \(T R_j\), for an
integer Laurent polynomial \(R_j\). This polynomial is uniquely
determined, because multiplication by the nonzero polynomial \(T\)
is injective.
An edge becomes

\[
R_j-R_k=\pm X^s.
\]

Thus the component is an induced graph on finitely many vertices of
an integer coordinate grid. Different components may be embedded in
one larger grid by adding a separating coordinate with widely spaced
constant values.

For at most eight vertices, the following bounds hold for the number
of edges in any integer-grid graph:

\[
\begin{array}{c|rrrrrrrr}
n&1&2&3&4&5&6&7&8\\ \hline
f(n)&0&1&2&4&5&7&9&12.
\end{array}
\tag{34}
\]

This small table is an induction, not an experimental enumeration.
Split the vertices by an integer coordinate threshold into nonempty
parts of sizes \(a,b\). Edges crossing that threshold form a matching,
so their number is at most \(\min(a,b)\).
The entries in (34) satisfy

\[
f(a)+f(b)+\min(a,b)\leq f(a+b)\qquad(a+b\leq8),
\]

which is checked directly for the finitely many displayed integer
splits. This proves (34) in every dimension.

### 7.3. The distance count excludes seven rows

If there were seven rows, let \(h\) be the number of distance-six pairs.
Equation (34) gives \(h\leq9\).
All other pair distances are even and at least eight, so

\[
\sum_{j<k}\|b_j-b_k\|_1
\geq6h+8(21-h)\geq150.
\tag{35}
\]

At a coordinate of width \(w\), each integer threshold separates at most
\(\lfloor7^2/4\rfloor=12\) pairs. Hence the same sum is at most
\(12\cdot12=144\), a contradiction.
This proves the finite-width theorem.
For the original eight-row question, the analogous contradiction is
\(224-2\cdot12=200>16\cdot12=192\).

### 7.4. Consequence for all quadratic-unit templates

The template (3) has at most eight persistent endpoint rows.
More precisely,

\[
\boxed{\#\{z_j\}\geq7\quad\Longrightarrow\quad K_0\leq4.}
\tag{36}
\]

To verify this, if \(K_0>4\), the parity proof in Section 6 still
applies. Outside the minimal norm-minus-one unit \(\lambda=\varphi\),
all pair distances are at least eight.
Seven rows would then have total distance at least \(21\cdot8=168\),
whereas their coordinate-width bound is at most \(12K_0\leq144\).
For \(\lambda=\varphi\), Section 7.3 gives the contradiction instead.
The row count (24) now proves

\[
\boxed{\#\{z_j\}\leq8}
\tag{37}
\]

for every fixed template in the theorem, without a restriction on the
quadratic unit's size.

If equality holds, the four active units of width must be four binary
coordinates. Indeed, the other partitions of total width four,
\((4),(3,1),(2,2),(2,1,1)\), allow at most \(3,4,5,6\) mutually
nonadjacent vectors respectively. These bounds follow by alternating
along a path through the grid vertices, whose sizes are respectively
\(5,8,9,12\). Such a path is obtained by the usual back-and-forth
traversal of successive rows of each rectangular grid.
Smaller total width is excluded by (24).

The eight rows in the binary four-cube form one entire parity class.
They are a largest independent set. Since the cube is regular, every
edge must join this set to its complement; the complement is therefore
also independent. The cube is connected, so this partition is its
unique bipartition, up to interchange.
Every growing active threshold layer consequently has four rows on
each side.

## 8. Almost all conductor layers of an eight-row template are balanced

The bounded oriented-gcd argument makes the last conclusion precise even
when evaluated block supports overlap at fixed primes.
Retain the Gaussian integer \(D_0\) from Section 1 and put
\(C_0=\prod_j c_j\).
For the primitive eight-row configuration, let \(W=\log R_n^2\), let
\(W_4\) be the weight of its four-versus-four split-prime threshold
layers, and let

\[
D=\sum_{p,\ell}(r_{p,\ell}-4)^2\log p.
\]

Then there is an explicit constant depending only on the fixed template,

\[
B_{\rm fix}=2\log N_{\mathbf Q(i)/\mathbf Q}(C_0D_0^4),
\]

such that throughout the endpoint family

\[
\boxed{0\leq W-W_4\leq B_{\rm fix},\qquad
       0\leq D\leq16B_{\rm fix}.}
\tag{38}
\]

Here and below, finitely many initial indices at which a block vanishes
are omitted.

To check this, fix a Gaussian prime \(\pi\) above a rational split prime
\(p\), and write

\[
b_\pi=\max_jv_\pi(c_j)+4v_\pi(D_0).
\]

Among the eight oriented block factors for the four binary coordinates,
choose one with maximal \(\pi\)-valuation \(v\).
Every other oriented factor has valuation at most \(v_\pi(D_0)\).
The four rows omitting the maximal factor have valuations in
\([0,b_\pi]\), and the four rows including it have valuations in
\([v,v+b_\pi]\). Subtracting the valuation of the common Gaussian gcd
merely shifts both intervals.
At most \(2b_\pi\) integer threshold layers can therefore be unbalanced:
outside those two bounded bands, a threshold separates exactly the four
high rows from the four low rows. This remains valid if the bands
overlap, because their entire valuation range then has length at most
\(2b_\pi\).

Summing over one orientation above each split prime gives
\[
W-W_4\leq
2\sum_\pi b_\pi\log p
\leq2\log N(C_0D_0^4)=B_{\rm fix}.
\]
Non-split factors disappear in the primitive normalization.
The defect coefficient of any unbalanced eight-row layer is at most
sixteen, proving the second inequality in (38).

Thus an eight-row persistent template has
\(D=O_{\rm template}(1)\) and \(W_4=W-O_{\rm template}(1)\).
This statement allows bounded overlap between the growing blocks and
the fixed prefactors; disjoint norm supports were not assumed.
The additive constants remain dependent on the fixed template.

## What does not transfer to arbitrary circles

The fixed-template assumptions supply one exponent polynomial \(P\),
one fixed phase expansion, and fixed bounds for all common-factor
losses. They are essential to comparing powers of \(\lambda^n\).
The proof excludes an oversized fixed template from meeting the endpoint
condition along unbounded indices. It does not give a uniform bound on
the largest exceptional index as the template changes.

In particular, the theorem does not allow replacing fixed shifts by
shifts depending on \(n\), allowing the private Gaussian factors to grow,
or independently varying several quadratic units.
No extraction of arbitrary endpoint lattice configurations into the
template (3) is established.

The argument is a prose proof. No Lean formalization of this theorem is
claimed.
