# Eight endpoint points in a shared-ray Pell product family

This note proves an asymptotic classification for a restricted class of
Pell block products. It is not a reduction of arbitrary Pell-linear
blocks to this class, and it is not a proof of the general mixed
eight-point bonus.

**Theorem.** Consider at most seven growing Pell-linear Gaussian
blocks whose leading coefficients lie on one common complex ray,
with arbitrary fixed Gaussian private factors. Assume the blocks
are individually conjugate-coprime and their rational norms have
disjoint supports. If eight fixed, distinct product rows lie in an
arc of length \(C\sqrt R\), for one fixed \(C\), along an unbounded
sequence of Pell indices, then there are exactly four growing
blocks and the eight rows form one full parity class of the
four-dimensional binary cube. In particular every growing block
has a balanced \(4\)-versus-\(4\) cut.

This excludes a mixed critical eight-row counterexample in this
class: any unbalanced conductor in such a persistent endpoint
family comes from its fixed private factors.

## 1. Exact hypotheses

Let \(K=\mathbf Q(\sqrt5)\), let \(\tau\) be its nontrivial
automorphism extended to \(K(i)\) by \(\tau(i)=i\), and put
\[
\lambda=9+4\sqrt5,\qquad t_n=\lambda^n.
\]
Fix a nonzero \(h\in K(i)\) and nonzero real
\(q_1,\ldots,q_m\in K\), where \(1\leq m\leq7\). The blocks are
\[
H_j(n)=h q_j t_n+\tau(h)\tau(q_j)t_n^{-1}\in\mathbf Z[i].
\tag{1}
\]
Integrality is required only on the chosen infinite progression of
indices. The common-ray hypothesis is precisely the factorization
of all leading coefficients as \(h q_j\), with real \(q_j\).
Arbitrary coefficients in \(K(i)\) need not have this form.

Assume, on an unbounded sequence of these indices,
\[
\gcd(H_j(n),\overline{H_j(n)})=1,\qquad
\gcd\bigl(N(H_j(n)),N(H_k(n))\bigr)=1\quad(j\ne k).
\tag{2}
\]
Gaussian gcd equal to \(1\) here means a unit. Each block grows
like a positive constant times \(t_n\).

Fix a set \(\mathcal C\subset\{-1,1\}^m\) of row patterns, and for
each row choose a fixed nonzero \(E_s\in\mathbf Z[i]\), all of the
same modulus \(\rho\). Define
\[
z_s(n)=E_s
 \prod_{j:s_j=1}H_j(n)
 \prod_{j:s_j=-1}\overline{H_j(n)}.
\tag{3}
\]
The common radius satisfies
\[
R_n=\rho\prod_j|H_j(n)|\asymp t_n^m.
\tag{4}
\]
The points are assumed distinct for the indices under discussion.
Allowing more than one fixed private factor for the same row
pattern cannot enlarge a cluster whose angular width tends to
zero: their quotient would be a fixed unit-modulus number
different from \(1\). Thus using a set of distinct patterns loses
nothing.

Suppose that for an unbounded sequence of indices all these
points lie in one arc with angular width
\[
\Delta_n\leq C R_n^{-1/2}=O(t_n^{-m/2}).
\tag{5}
\]
No assumption on the Gaussian units of the \(E_s\) is needed for
the classification.

## 2. Leading phases and first moments

Multiplying an individual block by \(-1\) changes every row by
the same sign, so take \(q_j>0\) in the chosen real embedding.
Write
\[
u=\frac h{\overline h},\qquad
\kappa=\frac{\tau(h)}h,\qquad
r_j=\frac{\tau(q_j)}{q_j}\in K^\times.
\tag{6}
\]
Then
\[
\tau(r_j)=r_j^{-1},\qquad
H_j(n)=h q_j t_n(1+\kappa r_j t_n^{-2}).
\]
For sufficiently large \(n\), continuous phase lifts give
\[
\arg H_j(n)=\arg h+
 \operatorname{Im}(\kappa)r_jt_n^{-2}+O(t_n^{-4}).
\tag{7}
\]
The constants may depend on the fixed family.

Let \(k(s)\) be the number of plus signs in a row \(s\).
Equality of the limiting phases of two rows implies
\[
\frac{E_s}{E_v}\,u^{\,k(s)-k(v)}=1.
\tag{8}
\]
If \(m>4\) and \(\operatorname{Im}\kappa\ne0\), the bound (5)
is \(o(t_n^{-2})\). Comparing the next coefficients in (7)
therefore gives the exact first-moment equality
\[
\sum_js_jr_j=\sum_jv_jr_j
\quad\text{for every }s,v\in\mathcal C.
\tag{9}
\]

The \(r_j\) are distinct. Otherwise \(r_j=r_k\) implies
\(q_j/q_k\in\mathbf Q\), hence \(H_j/H_k\) is a fixed nonzero
rational number. Their norm ratio is then a fixed positive
rational number. Two relatively prime positive integer norms
with a fixed rational ratio are bounded, contradicting (4).

## 3. The two possible leading-phase regimes

Suppose first that \(u^\ell\in\mathbf Q(i)\) for some nonzero
integer \(\ell\). Then
\[
v=\frac{\tau(u)}u
\]
is a root of unity, and \(\tau(v)=v^{-1}=\overline v\).
Consequently \(v=a+i b\sqrt5\) for rational \(a,b\).
Because \(v\) is a root of unity, \(2a=v+v^{-1}\) is a rational
algebraic integer, hence an integer. The identity
\(a^2+5b^2=1\) leaves only \(a=\pm1,b=0\):
the alternatives \(a=0,\pm1/2\) would require \(1/5\) or
\(3/20\) to be a rational square. Thus
\[
\tau(u)=u\quad\text{or}\quad \tau(u)=-u.
\tag{10}
\]

The first possibility cannot occur under (1)--(2). It would make
\(\kappa\) real, so every \(H_j(n)\) would lie on one fixed
complex line, with \(H_j/\overline H_j=u\in\mathbf Q(i)\).
Because the line contains a nonzero Gaussian integer, it has a
primitive Gaussian direction \(A+iB\). Every \(H_j(n)\) is an
integer multiple of that direction. As its magnitude tends to
infinity, its coordinate gcd tends to infinity, contrary to
conjugate coprimality.

We are left with two regimes:

1. **Quadratic phase regime:** \(\tau(u)=-u\). Here
   \(\kappa/\overline\kappa=-1\), so \(\kappa\) is nonzero and
   purely imaginary. Moreover \(u^d\in\mathbf Q(i)\) holds
   exactly when \(d\) is even. Equation (8) forces every
   \(k(s)\) to have the same parity.
2. **No Gaussian-rational power:** no nonzero power of \(u\)
   lies in \(\mathbf Q(i)\). Here
   \(\operatorname{Im}\kappa\ne0\), and (8) forces all
   \(k(s)\) to be equal.

In the first regime one also has
\[
r_j\ne-r_k\qquad(j\ne k).
\tag{11}
\]
Indeed such an equality would give \(q_j/q_k=c\sqrt5\) for a
nonzero rational \(c\). Directly from (1),
\[
H_j(n)=Q\,\overline{H_k(n)},\qquad
Q=c\sqrt5\,u\in\mathbf Q(i),
\tag{12}
\]
where the membership follows from \(\tau(u)=-u\).
Again their norm ratio would be fixed rational, contradicting
the unbounded pairwise-coprime norms in (2).

## 4. Two elementary subset-sum bounds

The following facts concern distinct real numbers \(r_j\in K\)
with \(\tau(r_j)=1/r_j\).

### 4.1. Positive coefficients

Suppose all \(r_j>0\). A relation between disjoint subsets,
\[
\sum_{j\in A}r_j=\sum_{j\in B}r_j,
\tag{13}
\]
also holds for their reciprocals.
Neither side can have just one term if the other has at least
two: the one term is larger than each term on the other side,
while its reciprocal is smaller than each of their reciprocals.
A \(2\)-versus-\(2\) relation is also impossible for disjoint
subsets: equality of sums and reciprocal sums gives equality
of products, so the two unordered pairs coincide.

Thus any two distinct binary sign rows with equal weighted
sums have Hamming distance at least \(5\). Three such rows
would have total pairwise Hamming distance at least \(15\).
Each coordinate contributes at most \(2\) to this total, so
three rows need at least \(8\) coordinates. For \(m\leq7\),
there can be at most two rows.

### 4.2. A fixed number of terms

Without assuming positivity, fix a cardinality \(k\).
For \(m\leq7\), any fiber of the map
\[
A\longmapsto\sum_{j\in A}r_j,\qquad |A|=k,
\tag{14}
\]
has at most three members.

Taking complements reduces the claim to \(k\leq3\).
The cases \(k=0,1\) are immediate. For pairs, a nonzero sum
and its reciprocal sum determine the product, and hence the
pair uniquely. Pairs of sum zero consist of opposite elements;
there are at most \(\lfloor m/2\rfloor\leq3\) such pairs.

For triples, two distinct members of a fiber cannot share two
indices. If they share exactly one index \(j\), the two
remaining pairs have the same sum. Pair uniqueness forces
both pairs to have sum zero. The common triple sum is therefore
\(r_j\), and \(j\) is uniquely determined by that sum.
Every further triple in the fiber either contains \(j\), or
must be disjoint from both of the first two triples. The latter
is impossible: their union already uses five of at most seven
indices. All triples in the fiber then consist of \(j\) and an
opposite pair, giving at most three triples.
If no two triples intersect, there are at most two.

## 5. Proof of the classification

For \(m\leq3\), the quadratic phase regime permits at most one
parity class of patterns, of size \(2^{m-1}\leq4\).
The no-power regime permits at most one constant-weight layer,
also of size at most \(3\). Eight rows are impossible.

For \(m=4\), the no-power regime has at most
\(\binom42=6\) rows. The quadratic phase regime permits at most
the eight patterns in one parity class. If there are eight,
that class is complete. Each of the four coordinates is plus
in exactly four rows and minus in four rows.

Finally let \(5\leq m\leq7\). Equation (9) holds.
In the quadratic phase regime, (11) permits replacing each
negative \(r_j\) by \(-r_j\) and simultaneously reversing that
coordinate in every row. The resulting positive coefficients
remain distinct, still have norm one, and preserve all
equalities of weighted sign sums. Section 4.1 bounds the
number of rows by two.
In the no-power regime, all rows have the same weight.
Equation (9) is equality of their plus-subset sums, so
Section 4.2 bounds the number of rows by three.
This proves the theorem.

There is a further exclusion in the quadratic phase regime alone:
eight endpoint rows are impossible for every \(5\leq m\leq10\).
Their weights have one common parity, so their pairwise Hamming
distances are even. Section 4.1 consequently improves the minimum
distance from \(5\) to \(6\). The sum over the \(28\) pairs of
eight rows would be at least \(28\cdot6=168\). Each binary
coordinate contributes at most \(4\cdot4=16\), giving at most
\(16m\leq160\), a contradiction. The fixed-weight no-power
regime for \(m=8,9,10\) is not covered by this extra argument.

## 6. Consequence for the conductor bonus, and its limits

Suppose additionally that the fixed factors \(E_s\) are represented
by one fixed inherited Gaussian conductor of rational norm \(N_0\),
whose rational prime support is disjoint from every growing block.
This is the usual fixed-private-factor product model.
Let \(W_0=\log N_0\), and write \(W_{\rm grow}\) for the total
log norm of the variable blocks.

In an eight-row endpoint family covered by the theorem, every
variable layer has size \(4\). Thus all its defect is fixed:
\[
D=D_0,\qquad
W_4=W_{\rm grow}+W_{0,4},\qquad
W=W_{\rm grow}+W_0.
\]
For each fixed conductor layer the coefficient
\[
(r-4)^2+\mathbf1_{r=4}-2
\]
is at most \(14\). Therefore the desired coefficient-one bonus
holds throughout that product family with the explicit fixed bound
\[
\boxed{D+W_4\leq2W+14\log N_0.}
\tag{15}
\]
If \(N_0\) is bounded in advance, so is this additive constant.

This conclusion concerns families whose endpoint property
persists along unbounded Pell indices. It does not bound isolated
endpoint configurations before the leading-phase asymptotics apply,
and it does not give an additive constant independent of unrestricted
private height.

The same-ray assumption is substantive. General Pell-linear
blocks can approach different rays, including powers of one ray;
their limiting phases impose different weighted relations.
For variable degrees \(m=8,9,10\), Section 5 additionally excludes
eight rows in the quadratic phase regime. The no-power regime at
these degrees, degrees \(m\geq11\), variable private factors,
and different-ray systems remain outside this classification theorem.
Both explicit four-point counterexamples in the neighboring
notes satisfy the shared-ray hypothesis with \(m=3\), which is
consistent with the allowed four-row parity class.

[The higher odd-moment continuation](odd_moment_separation.md)
proves further coefficient-count ranges and excludes mixed critical
eight-row families of every degree in the narrower pure-Pell-shift
subclass.
