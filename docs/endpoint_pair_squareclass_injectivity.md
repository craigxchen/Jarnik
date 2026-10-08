# Pair-squareclass injectivity through the direct constant C = 2

The argument is valid for a fixed Gaussian-unit class of an actual
short-arc cluster. It uses angular information absent from the allocation
countermodel in `pair_norm_squareclass_code_audit.md`.

**Relation to the existing results.** Section 3 of
[squareclass_baseline_and_power_groups.md](squareclass_baseline_and_power_groups.md)
already proves four-wise augmented-label separation for C<1, using the
same three balanced sign partitions and handling zero words. The present
proof gives a direct fixed-unit version for C<=2. It extends the range of
that elementary baseline; it does not provide a new qualitative packing
bound for arbitrary fixed C, since the earlier result already allows
subdivision. No improvement of the general growth rate follows.

Suppose distinct Gaussian integers \(z_i\) have common squared modulus
\(N\), arguments \(\theta_i\) in an interval of width
\(\Delta\le C N^{-1/4}\), and \(0<C\le2\). In one fixed unit class write
the split-prime allocation exponents as \(e_{i,p}\in[0,n_p]\). Here \(N\) is the squared radius of the full source cluster, not a radius
obtained by dividing out the common Gaussian gcd of a selected quartet.
Common Gaussian factors cancel from every ratio below; at each prime the
selected quartet still has allocation range at most the full exponent
\(n_p\). Their presence therefore only weakens the upper bounds expressed
using the original \(N\).

For a pair \(i,j\), the primitive Gaussian quotient factor has norm
\(d_{ij}=\prod_p p^{|e_{i,p}-e_{j,p}|}\). Its squareclass is the parity
vector \(\rho_i+\rho_j\), where \(\rho_i=(e_{i,p}\bmod2)_p\).
We prove that these squareclasses are distinct over all unordered pairs.

## Four distinct indices

Suppose two disjoint edges have the same squareclass. Relabel their four
vertices \(0,1,2,3\). At every prime, \(e_0+e_1+e_2+e_3\) is even, so
\[
 s_1=(e_0+e_1-e_2-e_3)/2,\quad
 s_2=(e_0+e_2-e_1-e_3)/2,\quad
 s_3=(e_0+e_3-e_1-e_2)/2
\]
are integers. Build a primitive Gaussian integer \(U_j\) by choosing, at
each prime, the Gaussian orientation prescribed by the sign of \(s_j\)
and exponent \(|s_j|\). Its norm is \(D_j=\prod_p p^{|s_{j,p}|}\).

For the partition \(ab\mid cd\) associated to \(j\), exact factorization
gives
\[
 (U_j/\overline U_j)^2=z_a z_b/(z_c z_d).
\]
Consequently
\(\arg U_j=(\theta_a+\theta_b-\theta_c-\theta_d)/4\) modulo \(\pi/2\).
Multiplication by a Gaussian unit permits choosing exactly this argument,
which has absolute value at most \(\Delta/2\). This choice is legitimate:
multiplication by \(i\) changes the phase \(U_j/\bar U_j\) by a sign,
and leaves its square unchanged.

None of these three chosen Gaussian integers is real. Otherwise its squared
phase would be one, giving a nontrivial equality of two disjoint source
pair products. The established [multiplicative-rectangle separation theorem](multiplicative_rectangle_separation.md)
excludes such an equality whenever \(\Delta N^{1/4}<2\sqrt2\), which holds
here. Thus \(|\Im U_j|\ge1\), and
\[
 D_j\ge\frac1{\sin^2(\Delta/2)}
 >\frac4{\Delta^2}\ge\frac4{C^2}\sqrt N\ge\sqrt N.
\]
The strict inequality follows from \(\sin x<x\), and remains strict at
\(C=2\).

On the other hand, at each prime,
\[
 |s_1|+|s_2|+|s_3|
 \le\frac32(\max_i e_i-\min_i e_i)
 \le\frac32 n_p.
\]
One proof is to subtract the minimum, divide by the positive range, and
observe that the left side is convex on \([0,1]^4\); its maximum occurs at
a vertex. If that vertex has respectively zero, one, two, three, or four
coordinates equal to one, the half-sum totals are respectively
\(0,3/2,1,3/2,0\). Rescaling gives the asserted bound. The expression is invariant under
subtracting a common constant. Zero range is immediate.
Multiplying over primes gives
\(D_1D_2D_3\le N^{3/2}\), contradicting the three strict lower bounds.

## Shared endpoints and zero labels

The existing nonsquare pair-norm conclusion already excludes adjacent-edge
collisions: equal labels on \(ij\) and \(ik\) would make \(d_{jk}\) a square.
For completeness the following direct half-difference proof establishes the
same fact, including the boundary constant.

If two edge squareclasses with a shared endpoint agree, the other two rows
have identical parity vectors. More generally suppose \(\rho_i=\rho_j\)
for distinct rows. The integer half-differences
\((e_{i,p}-e_{j,p})/2\) define a Gaussian integer \(U\) with
\[
 N(U)\le\sqrt N,\qquad (U/\bar U)^2=z_i/z_j.
\]
A Gaussian unit again gives \(\arg U=(\theta_i-\theta_j)/4\), of absolute
value at most \(\Delta/4\). This is nonzero because the rows are distinct
and the arc has width less than \(\pi\). Therefore
\[
 N(U)\ge\frac1{\sin^2(\Delta/4)}
 >\frac{16}{\Delta^2}\ge4\sqrt N,
\]
a contradiction. Hence labels are nonzero and a repeated edge label cannot
share an endpoint. The preceding four-row argument handles the remaining
possibility.

Thus, for \(m\) points in the retained unit class and \(\omega\) split
primes, the actual angular hypotheses give
\[
 \binom m2\le2^\omega-1.
\]
The bound also holds with \(\omega\) replaced by the binary dimension of
the allocation labels, so the retained labels form a Sidon set in their
binary ambient group. A general cluster may first be partitioned into at
most four unit classes.
This conclusion depends on angular separation and does not follow from
ordinary norm injectivity or allocation constraints alone. It also does not
bound \(m\) uniformly when \(\omega\) is unrestricted.

An independent finite check tested the half-sum inequality on all 5,000
quadruples in \(\{0,\ldots,9\}^4\) having even total sum; every case passed.
The convexity proof establishes it at arbitrary real values and hence at
all allocation depths.

The persistent [checker](check_angular_pair_squareclass_injectivity.py)
reproduces those 5,000 inequalities and checks 360 exact Gaussian root
identities, including all four unit choices and configurations with a
nontrivial common source factor. The analytic argument establishes the
strict inequality at C=2 independently of these finite checks.
