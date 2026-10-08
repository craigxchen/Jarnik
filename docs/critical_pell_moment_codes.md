# Eleven-factor slack and a twelve-factor Hadamard obstruction

This note concerns the first-moment stage of the positive, shared-ray,
quadratic Pell model. It proves two finite obstructions; it does not
settle all twelve-factor codes or the general eight-point bonus.

Let \(S\) be an \(8\)-by-\(m\) sign matrix with distinct rows.
Assume its rows have one common parity and that positive distinct
coefficients \(r_1,\ldots,r_m\) satisfy
\[
Sr=\alpha\mathbf1,\qquad S(r^{-1})=\beta\mathbf1,
\tag{1}
\]
where reciprocal is taken coordinatewise.
In the Pell application, \(r_j\in\mathbf Q(\sqrt5)\) have norm one,
so the second equality is the conjugate of the first.

## 1. The eleven-factor case has a strict conductor margin

Two rows of (1) cannot have Hamming distance two or four.
After canceling their common terms, a one-term side cannot equal a
sum of at least two positive terms both for the coefficients and
their reciprocals. A two-versus-two relation would give the same
sum and reciprocal sum for two disjoint pairs, hence the same
product, and would make those unordered pairs identical.
The common row parity makes all distances even, so every pair
has Hamming distance at least six.

For a column with \(u_j\) plus signs, write
\[
\Delta=\sum_{j=1}^m(u_j-4)^2,\qquad
b=\#\{j:u_j=4\}.
\]
Counting separated pairs gives
\[
\Delta=16m-\sum_{i<k}d(S_i,S_k).
\tag{2}
\]
At \(m=11\), this yields
\[
\Delta\leq176-28\cdot6=8,\qquad
\Delta+b-2m\leq8+11-22=-3.
\tag{3}
\]
Thus the coefficient-one bonus has a strict growing margin.
More precisely, if the eleven block weights are
\(2\log t+O(1)\), fixed private factors contribute only \(O(1)\),
and the growing supports give the displayed column cuts, then
\[
D+W_4-2W
=2(\Delta+b-22)\log t+O(1)
\leq-6\log t+O(1).
\tag{4}
\]

For \(m=12\), equality in the minimum-distance bound means that
all pair distances are six. Its eight sign rows are orthogonal,
and (2) gives \(\Delta=24\). The coefficient of the bonus excess
then equals \(b\). This explains why a twelve-factor equality code
with balanced columns is a natural test case.

## 2. No reciprocal moment solution on eight rows of a Hadamard matrix

**Theorem.** Let \(S\) consist of any eight rows of any real
Hadamard matrix of order twelve. Row and column signs may be
chosen arbitrarily. There is no vector \(r\in\mathbf R_{>0}^{12}\)
with pairwise distinct entries satisfying (1).

In fact, the proof only uses that the signed vector introduced
below has pairwise distinct absolute values. No algebraicity
or norm-one hypothesis is needed.

*Proof.* Let \(K\) be the four omitted rows. Multiply every column
of the complete Hadamard matrix by the entry of one fixed omitted
row. That omitted row becomes the all-ones row. At the same time
replace \(r\) by the signed vector \(x\) obtained by the same
column signs. Then
\[
Sx=\alpha\mathbf1,\qquad S(x^{-1})=\beta\mathbf1,
\tag{5}
\]
and the \(|x_j|\) remain positive and pairwise distinct. We use
the same letters \(S,K\) for the transformed matrices.

Let the other three omitted rows have column patterns
\(\varepsilon=(\varepsilon_1,\varepsilon_2,\varepsilon_3)\).
Their individual sums and their pairwise correlations vanish.
If \(T\) is their triple correlation, the number of columns with
pattern \(\varepsilon\) is therefore
\[
n_\varepsilon=\frac{12+T\varepsilon_1\varepsilon_2
\varepsilon_3}{8}.
\tag{6}
\]
Nonnegative integrality forces \(T\in\{-12,-4,4,12\}\).
Thus the columns split either into four groups of size two and
four singleton groups, or into four groups of size three.

Consider two columns \(j,k\) in the same group. All four omitted
row entries agree in these columns, so
\(e_j-e_k\) is orthogonal to the row space of \(K\), and hence
lies in the row space of \(S\). Because \(SS^\mathsf T=12I\),
equation (5) gives
\[
x_j-x_k=\frac{\alpha u_{jk}}{12},\qquad
x_j^{-1}-x_k^{-1}=\frac{\beta u_{jk}}{12},
\quad
u_{jk}=\sum_{i=1}^8(S_{ij}-S_{ik}).
\tag{7}
\]
Two columns of a Hadamard matrix of order twelve differ in
exactly six entries. All these differences occur in the eight
selected rows, since the omitted entries agree. Their six
nonzero differences are each \(2\) or \(-2\), so
\[
u_{jk}\in\{0,\ \pm4,\ \pm8,\ \pm12\}.
\tag{8}
\]

Distinct absolute values exclude \(u_{jk}=0\), \(\alpha=0\),
and \(\beta=0\). Combining the two equations in (7) then gives
\[
x_jx_k=-\frac{\alpha}{\beta}=:p
\tag{9}
\]
for every pair in a common omitted-pattern group.
If a group has three columns, applying (9) to two pairs sharing
one column immediately makes the other two entries equal,
a contradiction.

It remains to consider four disjoint groups of size two.
There are only three possible nonzero values of \(|u_{jk}|\),
so two of these pairs have the same \(|u_{jk}|\).
Equation (7) gives the same squared difference for those pairs,
and (9) gives the same product \(p\).
The unordered squares of the two entries are determined by
\[
x_j^2+x_k^2=(x_j-x_k)^2+2p,\qquad
x_j^2x_k^2=p^2.
\]
Thus the two disjoint pairs have the same unordered absolute
values, again contradicting the distinctness of the \(|x_j|\).
This proves the theorem. \(\square\)

## 3. Scope of the critical obstruction

The theorem excludes every twelve-factor distance-six equality
code obtained by selecting eight rows from a complete Hadamard
matrix of order twelve, including signings with no constant
column. It therefore excludes positive distinct norm-one Pell
coefficients for all such codes.

A normalized Hadamard submatrix with a constant column already
has a simpler limitation: that growing factor is common to all
rows, so primitive normalization removes it and returns to at
most eleven active factors. The theorem above does not rely on
the presence of such a column.

The elementary proof in Section 2 assumes an actual completion;
it does not claim a completion theorem for arbitrary
eight-by-twelve partial Hadamard matrices. Section 4 removes
that hypothesis using the rank-four complement. Codes with
some distances larger than six still require other arguments.

Distinct coefficient values are a substantive hypothesis.
Expanding a repeated block into several unary columns produces
repeated coefficients, so this proof does not apply to that
multiplicity model. The fixed-shift multiplicity obstruction in
pell_product_endpoint_obstruction.md uses a different polynomial
argument and remains logically separate.

## 4. A completion-free theorem for every orthogonal eight-row code

The Hadamard completion hypothesis can be removed.

**Theorem.** Let \(S\in\{-1,1\}^{8\times12}\) satisfy
\(SS^\mathsf T=12I_8\). No positive vector with twelve distinct
entries satisfies (1).

### 4.1. The rank-four complement gives a root system

Set
\[
G=\frac{12I_{12}-S^\mathsf TS}{4}.
\tag{10}
\]
It is positive semidefinite, has diagonal entries one, and has
four nonzero eigenvalues, all equal to three. Thus there is a
matrix \(V\in\mathbf R^{4\times12}\), with unit columns \(v_j\),
such that
\[
V^\mathsf TV=G,\qquad VV^\mathsf T=3I_4,\qquad VS^\mathsf T=0.
\tag{11}
\]
A dot product of two sign columns of \(S\) is an even integer.
Positivity of the corresponding two-by-two principal minor of
\(G\) bounds its absolute value by four. Hence
\[
\langle v_j,v_k\rangle\in
\{0,\ \pm\tfrac12,\ \pm1\}\quad(j\ne k).
\tag{12}
\]

The vectors \(w_j=\sqrt2\,v_j\) generate a positive definite
even integral lattice of rank four. To check discreteness,
choose four independent \(w_j\)'s as a real basis. Every
lattice vector has integral pairings with that basis, so its
coordinates lie in the inverse of a fixed integral Gram
matrix applied to \(\mathbf Z^4\), a discrete set. Evenness
follows from the squared generator norms two and integral
mutual pairings.

The squared-norm-two vectors of this lattice form a finite
simply laced root system: reflection in such a vector \(w\)
is \(z\mapsto z-\langle z,w\rangle w\), which preserves the
lattice and the norm. These roots span the lattice because
the original generators have squared norm two.
The classification of simply laced root systems therefore
allows only irreducible components
\[
A_1,\ A_2,\ A_3,\ A_4,\ D_4
\]
in rank at most four. This is the rank-four specialization of
[Etingof, *Lie Groups and Lie Algebras*, Theorem 23.7 and
Section 23.9, pages 123 and 128](https://ocw.mit.edu/courses/18-755-lie-groups-and-lie-algebras-ii-spring-2024/mit18_755_s24_lec_full.pdf).
The standard root realizations used below have squared root
norm two.

### 4.2. The frame equation classifies parallel-column groups

For a component of type \(A_d\), realize its unit root lines as
\[
(e_i-e_j)/\sqrt2,\qquad 1\leq i<j\leq d+1,
\]
in the subspace whose coordinate sum is zero.
Let \(n_{ij}\) count the input columns on a given line, including
both orientations. The restriction of the frame equation in
(11) to this component makes the weighted graph Laplacian
\[
\sum_{i<j}n_{ij}(e_i-e_j)(e_i-e_j)^\mathsf T
=6\left(I-\frac{J}{d+1}\right).
\tag{13}
\]
Comparing off-diagonal entries gives
\(n_{ij}=6/(d+1)\) for every pair.
Thus \(A_3,A_4\) are impossible. An \(A_1\) component gives
three parallel input columns, and an \(A_2\) component gives
three parallel pairs. In rank four, the cases built only
from these components consequently contain a triple or,
for \(A_2\oplus A_2\), six disjoint parallel pairs.

For type \(D_4\), the unit root lines are
\[
(e_i+e_j)/\sqrt2,\quad(e_i-e_j)/\sqrt2,
\qquad1\leq i<j\leq4.
\]
The off-diagonal frame equations force the two lines belonging
to each pair \(i,j\) to have the same multiplicity, say \(n_{ij}\).
The diagonal equations then say
\[
\sum_{j\ne i}n_{ij}=3\quad(1\leq i\leq4).
\]
They imply equal multiplicities on opposite edges:
\[
n_{12}=n_{34}=a,\quad
n_{13}=n_{24}=b,\quad
n_{14}=n_{23}=c,\qquad a+b+c=3.
\tag{14}
\]
Up to permutation, the possibilities are:

- \((3,0,0)\): four parallel triples;
- \((2,1,0)\): four parallel pairs and four singleton lines;
- \((1,1,1)\): every one of the twelve \(D_4\) root lines once.

### 4.3. The twelve-distinct-lines case cannot have eight sign rows

For the final case in (14), orient and label the columns so that
the two columns belonging to an edge \(i<j\) are
\((e_i+e_j)/\sqrt2\) and \((e_i-e_j)/\sqrt2\).
Consider any sign vector \(s\in\{-1,1\}^{12}\) in the kernel
of \(V\), as every row of \(S\) must be by (11).

If the entries of \(s\) on this edge are \(x,y\), its contribution
is \(x+y\) at vertex \(i\) and \(x-y\) at vertex \(j\), ignoring
the common factor \(1/\sqrt2\). Exactly one endpoint receives
a nonzero contribution, equal to \(2\) or \(-2\). Assign the
edge to that endpoint. The zero-sum condition at each vertex
makes its number of assigned edges either zero or two. Since
there are six edges, the four counts must be \(0,2,2,2\).
Thus every kernel sign vector has one distinguished vertex
receiving no edges.

For a fixed distinguished vertex, all its three edges are
assigned to the other endpoints. On the remaining triangle,
each vertex must receive one edge, giving one of two directed
cycles. At each of those three vertices, the two assigned
contributions have opposite signs.

Two kernel sign vectors with this same distinguished vertex
cannot be orthogonal. If their triangle cycles agree, their
dot product is four times a sum of three signs. If their cycles
are opposite, the three triangle-edge contributions to the
dot product vanish, and the result is two times a sum of three
signs from the three remaining edges. Neither expression is
zero.

There are only four choices of distinguished vertex, so at
most four kernel sign vectors can be mutually orthogonal.
This contradicts the eight rows of \(S\). The final case of
(14) is therefore impossible.

### 4.4. Parallel columns finish the reciprocal obstruction

In every remaining case there is a parallel triple or at least
four disjoint parallel pairs among the \(v_j\). Flip column signs
within each group so the vectors agree exactly. Apply the same
sign changes to \(S\) and to \(r\), obtaining a real vector \(x\)
with distinct positive absolute values and the two moment
equations (5).

For a pair in a group, \(V(e_j-e_k)=0\); thus \(e_j-e_k\) is in
the row space of \(S\). Also \(G_{jk}=1\), so the corresponding
sign columns have dot product \(-4\) and differ in six entries.
Equations (7)--(9) and their three possible nonzero absolute
differences apply unchanged. A triple forces equal entries;
four disjoint pairs force a repeated unordered pair of absolute
values. This proves the completion-free theorem.

The proof uses the rank-four Gram complement and the cited root
classification, never a Hadamard completion assertion. Its final
finite combinatorial step was also checked by enumerating the
64 sign vectors in the full \(D_4\)-line kernel; their maximum
orthogonal subset has size four. The argument in Section 4.3
is the proof, and does not depend on that enumeration.

## 5. What remains after the equality-code exclusion

All orthogonal eight-by-twelve sign matrices are now excluded
from the positive distinct reciprocal-moment problem. Thus the
distance-six equality case cannot provide the proposed critical
Pell family.

The twelve-factor cases with some distances larger than six
remain outside this obstruction, as do unary expansions having
repeated coefficient values. No bound for those cases, or for
arbitrary Gaussian endpoint configurations, has been inferred.

## 6. A nonorthogonal critical code and its separate obstruction

The remaining combinatorial regime is not empty. The following sign
matrix has all columns active, one common row parity, twenty-seven
pair distances equal to six, and one pair distance equal to eight:
\[
S=\begin{pmatrix}
+&-&-&+&-&-&-&+&+&+&-&+\\
+&+&-&-&+&-&-&-&+&+&+&-\\
+&-&+&-&-&+&-&-&-&+&+&+\\
-&+&-&+&-&-&-&-&-&-&+&+\\
+&+&+&-&+&-&-&+&-&-&-&+\\
+&+&+&+&-&+&-&-&+&-&-&-\\
+&-&+&+&+&-&+&-&-&+&-&-\\
+&-&-&+&+&+&-&+&-&-&+&-
\end{pmatrix}.
\tag{15}
\]
Its column plus counts are
\[
(7,4,4,5,4,3,1,3,3,4,4,4),
\]
so
\[
\Delta=22,\qquad b=6,\qquad \Delta+b-24=4>0.
\tag{16}
\]
Thus minimum distance and parity alone do not make the bonus
coefficient nonpositive once the orthogonal case is excluded.

This particular matrix nevertheless has a short additional
reciprocal-moment obstruction. Number rows and columns from zero.
Directly from (15),
\[
-S_0-S_2+S_3+S_4+S_5-S_6
=(-2,6,0,0,0,0,-2,0,0,-6,0,0).
\tag{17}
\]
The coefficients of the row combination sum to zero. Consequently
the two constant-row-sum equations would imply
\[
3(r_1-r_9)=r_0+r_6>0,\qquad
3(r_1^{-1}-r_9^{-1})=r_0^{-1}+r_6^{-1}>0.
\tag{18}
\]
The first inequality makes \(r_1>r_9\), contradicting the second.
This certificate excludes positive coefficients even if repetitions
are allowed.

The matrix can be constructed from the first eight nonconstant
rows of the normalized Paley Hadamard matrix of order twelve,
then reversing row three in columns zero and six. Its exact
distance multiset, column counts, and relation (17) were checked
by integer enumeration. Formula (18) is the arithmetic obstruction;
the matrix itself is not an example of admissible Pell coefficients.
Other nonorthogonal twelve-factor codes remain unresolved here.

## 7. The same nonorthogonal pattern is excluded under every column signing

The obstruction for the matrix in (15) is stronger than its particular
positive-coordinate certificate (18). It admits no nonzero real vector
with pairwise distinct absolute values satisfying both constant-row-sum
equations. Consequently every independent signing of its twelve columns
is excluded for positive distinct reciprocal coefficients.

### A signed four-term reciprocal lemma

Suppose nonzero real numbers satisfy

\[
a+b=c+d,\qquad a^{-1}+b^{-1}=c^{-1}+d^{-1}.
\tag{19}
\]

Their four absolute values cannot be pairwise distinct. If the common
sum is zero, then \(a=-b\) already repeats an absolute value. Otherwise
the reciprocal equation, written as
\((a+b)/(ab)=(c+d)/(cd)\), gives \(ab=cd\). The two unordered pairs
are then roots of the same quadratic and coincide.

More generally, let \(A\) be the matrix of row differences of any sign
matrix. If its real row span contains a nonzero vector proportional to
\(e_i+e_j-e_k-e_l\), with four distinct indices, then
\(Ax=A(x^{-1})=0\) is incompatible with pairwise distinct absolute
values. Arbitrary choices of the four signs give the same obstruction:
absorb those signs into the four numbers before applying (19).
This is an exact linear certificate, not a numerical feasibility test.

### An explicit certificate for all 4,096 column signings

For the matrix (15), the following identity uses zero-based row indices:

\[
-S_0-S_1-S_2+2S_4+2S_5-S_7
=6(e_1+e_2-e_9-e_{10}).
\tag{20}
\]

Its row coefficients sum to zero. Hence constant row sums for both
\(x\) and \(x^{-1}\) imply

\[
x_1+x_2=x_9+x_{10},\qquad
x_1^{-1}+x_2^{-1}=x_9^{-1}+x_{10}^{-1}.
\tag{21}
\]

The signed reciprocal lemma proves the stated collision. Identity (20)
is verified directly by integer addition of the displayed rows.

The algebraic signed reciprocal lemma is also checked in
[ReciprocalCollision.lean](../GaussianChain/ReciprocalCollision.lean),
as `pair_collision` and `abs_collision`. The displayed matrix certificate
and its application to Pell factors remain in the prose argument.

For a diagonal sign matrix \(E\), a positive solution for \(SE\) would
give \(x=Er\), with \(x^{-1}=E(r^{-1})\). Its absolute values would
be the distinct positive entries of \(r\), contradicting (21).
Row and column permutations also preserve this exclusion. Independent
row sign changes are not asserted to preserve the zero-sum coefficient
condition in (20), and are not covered by this extension.

The initial bounded search through small row combinations excluded
4,088 column signings by positivity and weighted reciprocal inequalities.
The four-term certificate closes the other eight without relying on
that enumeration. It proves a whole nonorthogonal pattern exclusion,
but not a classification of all nonorthogonal twelve-factor codes.

## 8. A cubic-moment obstruction after two arbitrary column deletions

There is a separate real-algebraic variant relevant to rational linear
Gaussian blocks. Let `H` be any Hadamard matrix of order twelve, with
arbitrary row and column signs. Select eight rows and delete any two
columns, obtaining `S` of size `8 x 10`. Then there is no real vector
`a` with ten nonzero, pairwise distinct absolute values such that

```text
S a = alpha 1,             S(a^3) = beta 1,             (22)
```

where the cube is coordinatewise. This statement does not require a
constant deleted column. Its proof below uses the actual Hadamard
completion; a general length-ten distance-five code is not assumed
to have such a completion.

### The two zero coordinates and a Bessel bound

Restore the deleted columns and put zero in both corresponding
coordinates of `a`, giving `x` of length twelve. Denote the selected
eight-row matrix by `T` and the omitted four-row matrix by `K`.
As in Section 2, orient columns to make one omitted row constant,
and absorb these signs into `x`. Odd powers respect that change,
so `T x=alpha 1` and `T(x^3)=beta 1` still hold. There are exactly
two zero coordinates; all ten nonzero absolute values remain distinct.
The omitted-pattern groups are four pairs and four singletons, or
four triples, by (6).

For columns `j,k` in one pattern group, (7) now reads

```text
x_j-x_k = alpha u_jk/12,
x_j^3-x_k^3 = beta u_jk/12,
u_jk in {0, +/-4, +/-8, +/-12}.                         (23)
```

One has `alpha!=0`: otherwise every repeated-pattern group has
equal entries, and at least one such group avoids both zeros,
contradicting the distinct nonzero absolute values. Whenever
`x_j!=x_k`, cancellation in (23) gives

```text
q := beta/alpha = x_j^2+x_j x_k+x_k^2.                  (24)
```

The values of `q` and `d^2=(x_j-x_k)^2` determine the unordered
absolute values of the pair, since

```text
x_j x_k=(q-d^2)/3,
x_j^2+x_k^2=(2q+d^2)/3.                                (25)
```

In the four-pair case, unequal pairs must therefore have distinct
nonzero values of `|u_jk|`. If both zeros are not in one pair, all
four pairs are unequal, contradicting the three choices `4,8,12`.
If the zeros occupy one pair, that pair has `u=0`, while the other
three must use all three nonzero choices. Let `v_p` be the difference
of the two selected sign columns in pair `p`. These four vectors
are orthogonal and have squared norm `24`: their omitted entries
cancel, and full Hadamard column orthogonality gives both assertions.
Consequently Bessel's inequality against the eight-dimensional
all-ones vector gives

```text
sum_p u_p^2 <= 24 ||1||^2 = 192.
```

But the required values give `0^2+4^2+8^2+12^2=224`, a contradiction.
This is the additional step needed when deleted coordinates can
belong to the same pattern pair.

### The triple groups

If the zeros lie in different triples, a triple containing one zero
has two nonzero entries `x,y`; (24) gives `x^2=q=y^2`, a forbidden
absolute-value collision. Otherwise the zeros lie in one triple.
Its third entry `t` obeys `q=t^2`. Put `h=|alpha|/3>0`.
Equation (23) gives `|t|/h in {1,2,3}`, so
`q/h^2 in {1,4,9}`.

Every other triple has three distinct nonzero entries. Equality of
(24) for its three pairs forces their sum to be zero. Their positive
consecutive gaps, in units of `h`, are either `(1,1)`, `(1,2)`, or
`(2,1)`, since all pair distances belong to `{h,2h,3h}`. The first
choice gives `(-h,0,h)` and is impossible. Either other choice gives
`q/h^2=7/3`, contradicting `{1,4,9}`. This proves (22)'s exclusion.

### What this proves for fixed rational linear blocks

For example, fix nonzero integers `a_j` with distinct absolute values
and take the ten Gaussian blocks `H_j(t)=t+i a_j`, with literal
product rows prescribed by `S`. Assume all ten retained columns are
nonconstant. The row units must agree for an arc whose angular width
tends to zero; equivalently they may already have been chosen to
align the limiting phases in the formulation `a_j+i t`.

The common Gaussian divisor of the product rows is bounded in `t`.
Indeed their polynomials in `Q(i)[t]` have no common root: the roots
`+/-i a_j` are all distinct, and each block occurs in both orientations
among the rows. A polynomial Bezout identity, with denominators cleared,
then bounds their common Gaussian divisor by one fixed nonzero Gaussian
integer. Thus the primitive radius is comparable to `t^10`.
The exact phase expansion is

```text
theta_i(t) = (S_i a)/t - (S_i a^3)/(3t^3) + O(t^-5).
```

An endpoint arc would require both displayed coefficients to be
independent of the row, contrary to (22). In fact its primitive
normalized arc constant tends to infinity at least as `c t^2`, with
`c>0` depending on the fixed template. This is an asymptotic family
exclusion, with no uniform threshold over moving coefficients. A
constant retained column must be removed and the active radius degree
recomputed before applying this particular degree-ten conclusion.

There is a limitation relevant to the intrinsic bonus: every eight-row,
length-ten code of minimum distance five already has a strict bonus
margin. Appending its parity bit gives length eleven and minimum
distance six, so its original defect satisfies `D<=176-168=8`.
With all columns active, write `n_j` for the number of cuts of minority
size `j`. Then `D=9n_1+4n_2+n_3`, while
`D+b-20=8n_1+3n_2-10`. Thus `n_1=0`, `n_2<=2`, and the excess is
at most `-4`. The cubic theorem excludes endpoint realization of the
stated Hadamard family; it does not remove a new obstruction to the
general eight-point bonus. The broader degree calculation is recorded
in [the rational moment continuation](odd_moment_separation.md).

[The exact cubic checker](check_critical_hadamard_cubic_moments.py)
audits every eight-row restriction of the order-twelve Paley matrix,
all selected-row signings, the projection and Bessel identities,
and the finite triple-gap calculation. These checks supplement the
proof for arbitrary Hadamard matrices of order twelve.

## 9. Even column counts yield short integer moment certificates

The following argument allows integer multiplicities and does not require
a Hadamard completion. It gives a family obstruction, not a classification
of all eight-point endpoint profiles.

### A weighted odd-moment certificate

Let `S` be an eight-row sign matrix, and suppose a vector `a` has nonzero,
pairwise distinct absolute values and satisfies

```text
S(a^(2r-1)) = alpha_r 1,             1<=r<=s.
```

There cannot be a zero-sum integer row vector `lambda` and a nonzero
rational scale `q` for which

```text
c=lambda^t S/q is integral,       0<sum_j |c_j|<=2s.     (26)
```

Indeed these hypotheses give `sum_j c_j a_j^(2r-1)=0`. Replace each
nonzero coefficient by `|c_j|` copies of `sign(c_j)a_j`. The resulting
multiset has size `h<=2s` and all its odd power sums through degree
`2s-1` vanish. Newton's identities imply that every odd elementary
symmetric function of degree at most `h` vanishes: in the identity
for an odd degree, each term has either a lower odd elementary
symmetric function or an odd power sum. If `h` is odd its product
therefore vanishes. If `h` is even its monic root polynomial is even,
so it has an opposite pair of roots. The first alternative contradicts
nonzero coefficients `a_j`; the second contradicts their distinct
absolute values, because all copies belonging to one coordinate have
the same sign. Repetition introduced by integer multiplicities causes
no exception.

### Orthogonal rows with even column counts

Let `T` be an `8 x L` sign matrix with `T T^t=L I`, where `4|L`.
Assume every column has an even number of plus signs. Delete any two
columns to obtain `S`. There is a balanced retained column `b`.
Indeed the full defect is `||T^t 1||^2/4=2L`, and every nonbalanced
column has defect at least four. Hence at least `L/2` columns are
balanced; `L>=8` by row rank, so at least two remain after deletion.
Set `lambda=T_b`, whose entries sum to zero. The dot product of any
two even-weight sign columns of length eight is divisible by four.
Thus

```text
c_j=(lambda^t S_j)/4 is integral,   c_b=2,
sum_j c_j^2 <= L/2,
0<sum_j |c_j| <= sum_j c_j^2-2 <= L/2-2.              (27)
```

The square bound is Parseval for the orthogonal rows of `T`:
`||T^t lambda||^2=8L`, followed by deletion of two nonnegative terms.
Every integer satisfies `|c_j|<=c_j^2`; the retained coordinate two
saves exactly two more units. Nonvanishing is immediate from that
coordinate. Consequently (26) excludes the corresponding odd-moment
system whenever `2s>=L/2-2`.

For fixed distinct-absolute rational linear blocks with all `L-2`
retained columns active, endpoint geometry forces
`s=floor((L-1)/4)=L/4-1` odd moments, exactly sufficient for (27).
The polynomial-gcd argument in Section 8 justifies
the primitive radius degree; no unbounded common factor is discarded.
This is an exclusion of this specified fixed-template family for
arbitrarily large `L`. It does not place arbitrary endpoint tuples
in that family or give a general point-count bound.

### A parity-padded extension with a small total Gram defect

There is also an exact variant before full row orthogonality. Let `S`
have `m` columns, let `s=floor((m+1)/4)`, and suppose its minimum
Hamming distance is at least `2s+1`. Append its parity column
`lambda`, then `h=4s+3-m` constant columns, so

```text
T=[S,lambda,1,...,1],       L=4s+4,       1<=h<=4.
```

The parity extension has even distances at least `2s+2`; therefore

```text
T T^t=L I-4A,
A_ii=0, A_ij nonnegative integers, K=sum_(i<j) A_ij.
```

Assume every original column has even plus-count and `lambda` is
balanced. Then `c=lambda^t S/4` is integral, the entries of `lambda`
sum to zero, and the constant columns have zero projection. Exactly

```text
sum_j c_j^2=L/2-4-lambda^t A lambda/4,
L/2-4-K/2 <= sum_j c_j^2 <= L/2-4+K/2.                (28)
```

Here `|lambda^t A lambda|<=2K`. If `K<=4` and `L-K>8`, equations
(28) and `|c_j|<=c_j^2` give the nonzero certificate (26), since
`2s=L/2-2`. More generally one may use the exact middle expression
in (28) whenever it is positive and at most `2s`. All parity,
column-count, and defect hypotheses are part of this corollary;
none has been inferred for arbitrary Gaussian configurations.

### The positive degree-twenty-two Paley profiles

An exact enumeration of all `binom(24,8)=735471` eight-row subsets
of the normalized order-twenty-four Paley matrix finds exactly `759`
subsets with a positive degree-twenty-two bonus after deleting the
constant column and one other column. Each has one constant column,
eight minority-size-two columns, and fifteen balanced columns.
The second deleted column must be one of the balanced columns.
Thus there are exactly `759*15=11385` such deleted codes, each with

```text
minority-size counts (0,8,0,14),
D=32, balanced count=14, D+14-44=2,
minimum pair distance=11.
```

These are genuine binary codes with a positive intrinsic-bonus
coefficient, not Gaussian endpoint examples. Formula (27) applies
to every one of them. More specifically the finite audit finds
an even shorter certificate by choosing `lambda` to be the deleted
balanced column: its remaining squared coefficient sum is
`(24*8-8^2)/16=8`, and every one of these `11385` certificates has
eight nonzero entries, all equal to `+1` or `-1`.
Odd moments through degree seven therefore exclude
all these fixed rational linear-block families; the endpoint
degree twenty-two condition supplies moments through degree nine.

[The stored exact fixture](paley24_even_column_moment_fixture.json)
contains the matrix, all `759` positive row sets, and one displayed
certificate. [The checker](check_hadamard_even_column_moment_certificate.py)
reconstructs the complete enumeration, checks all `11385` certificates,
and verifies the signed energy identity in (28) on actual parity-padded
sign matrices. No assertion that every critical degree-twenty-two code
is a Paley restriction is made.
