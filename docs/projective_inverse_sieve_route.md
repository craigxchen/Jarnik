# Projective inverse sieve: normalization and interpolation costs

The uniform endpoint bound remains unproved. This note proves a finite
projective interpolation bound and an exact limitation of adaptive
multiplicity choices. **The actual radial endpoint configuration has
height R, not R^(1/2); no transfer into the abstract fourth-root-height
regime is claimed.** These are prose proofs.

## 1. The literal two-residue model

Suppose N=R^2 is squarefree and odd. Each circle point z_i=x_i+i y_i
gives a primitive integer vector (x_i,y_i), since a common integer
divisor would have its square dividing N. At a prime p|N, the equation

    x_i^2+y_i^2=0 mod p

therefore places [x_i:y_i] in the two isotropic classes
[iota_p:1], [-iota_p:1], with iota_p^2=-1 mod p.
This is exact for squarefree norms. Higher powers require individual
coordinate contents and valuation bookkeeping, which are not removed
by merely writing projective coordinates.

An integral matrix A preserves these two distinct classes at primes
not dividing det A. Thus, after a rational projective change represented
by a primitive integral matrix, the conductor retained by this argument is

    Q_good = product_(p|N, p not dividing det A) p.        (1)

The coordinate distinction is important: [x_i:y_i] is the **radial**
projective direction. The rational half-angle parameter used in
`fresh_algebraic_parallel.md` parametrizes the circle itself. The
map from half-angle parameter to radial direction has degree two.
The two-residue condition above is not automatically a two-residue
condition on those half-angle parameters.

## 2. Three-anchor normalization loses every conductor prime

If a rational projective map sends three distinct radial directions
to 0,1,infinity, its primitive integral matrix A satisfies

    N divides det A.                                      (2)

Indeed, if its reduction were invertible at any p|N, it could not send
at most two anchor residue classes to the three distinct classes
0,1,infinity. Thus all conductor primes are bad for this normalization.
No positive good-prime conductor can silently be retained afterward.

Two-anchor normalization has an exact smaller cost. For primitive
vectors v_a,v_b, write D_ij=det(v_i,v_j) and use the map

    v -> [det(v_a,v):det(v_b,v)].

Its determinant is D_ab up to sign. A conductor prime survives exactly
when the two anchors reduce to different classes. The two surviving
classes become 0 and infinity. If a third point maps to a reduced
fraction [u:v], then

    Q_good divides u v.                                   (3)

Here u,v are nonzero, so max(|u|,|v|)>=sqrt(Q_good).
In the uniform bit-pattern model, the retained logarithmic conductor
is one half of the original. Fixing the two residue values therefore
does not give a free gain in a height-versus-conductor comparison.

This is not an endpoint height bound. Literal radial determinants
can have size R L, and the contents divided out of the transformed
coordinate pairs must still be tracked.

## 3. A finite interpolation theorem with projective denominators

Let S consist of distinct rational projective points represented by
primitive integer pairs (x,y) with max(|x|,|y|)<=H. Suppose that at
each prime dividing a squarefree integer Q>1, S reduces into a
prescribed union of two distinct projective classes. The classes
may vary arbitrarily with the prime.

For integers d>=1 and r>=0 with 2r<=d+1, there is a nonzero
homogeneous integer binary form F of degree d such that

    max absolute coefficient of F <= Q^(r(r+1)/(d+1)),
    Q^r divides F(x,y) for every point of S.                (4)

Consequently,

    (d+1) H^d < Q^(r-r(r+1)/(d+1))                         (5)

implies |S|<=d.

For the proof, work over Z_p and use an invertible local coordinate
change taking the two classes to X=0 and Y=0. In the expansion

    F(X,Y)=sum_(j=0..d) c_j X^j Y^(d-j),

impose p^(r-j)|c_j for j<r, and
p^(r-(d-j))|c_j when d-j<r. The ranges are disjoint, and these
conditions ensure p^r-divisibility on either whole residue class.
Their coefficient-lattice index is exactly

    p^[2(1+...+r)] = p^(r(r+1)).                           (6)

The induced binary-form coordinate change is invertible over Z_p,
so the index is chart-independent. The global lattice has index
Q^(r(r+1)). Pigeonholing coefficient vectors in an integer cube
provides a nonzero vector with the coefficient bound in (4).
Under (5), the integer F(x,y) has absolute value less than Q^r
and is divisible by Q^r, so it vanishes. A nonzero homogeneous
binary form of degree d has at most d rational projective zeros.

The height exponent is

    h_(d,r)=r(d-r)/(d(d+1)),
    max_r h_(d,r)=floor(d^2/4)/(d(d+1)) < 1/4.             (7)

It tends to 1/4 from below. Counting the coefficient-lattice index
as Q^(2r) would create a false improvement: several derivative
congruences modulo p do not make all evaluations divisible by p^r.

## 4. Multiplicities adapted to every cut still supply no saving

Suppose there are k candidate rows. At a prime with s rows in the
first class, allow arbitrary nonnegative multiplicities a,b with
a+b<=d+1. They can depend on the prime and its entire cut.

The index exponent at that prime is

    E(a,b)=[a(a+1)+b(b+1)]/2,                              (8)

while mean evaluation divisibility across the rows is
(s/k)a+(1-s/k)b. Define the convex function

    f_d(u)=max_(a,b>=0 integers, a+b<=d+1)
           [u a+(1-u)b-E(a,b)/(d+1)].                     (9)

For uniform logarithmic weights on all 2^k oriented bit patterns,
every d<k satisfies the exact bound

    average_cut f_d(|cut|/k) <= d/4.                      (10)

To prove it, let X_1,...,X_n be independent fair bits. Their mean
is the average of the n leave-one-out means. Jensen's inequality
therefore shows that

    E f_d((X_1+...+X_n)/n)

is nonincreasing in n. Reduce to n=d+1. At u=s/n the quadratics
in (9) are maximized by a=s or s-1, and b=n-s or n-s-1,
with appropriate nonnegative choices. The degree constraint is
satisfied. Thus

    f_d(s/n)=[s(s-1)+(n-s)(n-s-1)]/(2n).

Its fair-binomial expectation is (n-1)/4=d/4, proving (10).

The smallest row divisibility is at most its mean. Thus the
generic coefficient bound from the lattice index, even with all
prime-varying and cut-dependent multiplicities, cannot force
simultaneous vanishing at a height exponent greater than 1/4.
This goes beyond optimizing one constant multiplicity.

### Intrinsic normalization and nearly uniform weights

Remove the two constant cuts and give all remaining cuts equal
weight. Since f_d(0)=f_d(1)=d/2, (10) becomes

    average_cut f_d(|cut|/k) <= d beta_k,
    beta_k=(2^(k-2)-1)/(2^k-2).                           (11)

Equality holds when d=k-1. This is exactly the pair-determinant
averaging threshold: a fixed pair agrees with probability
2 beta_k under a uniform nonconstant cut.

If each cut weight differs from the uniform weight by relative
error at most eta, the bound changes by at most eta d/2, since
0<=f_d<=d/2. Hence there is no fixed positive gain as the profiles
tend to uniformity.

Exact arithmetic checks for k=4,6,8,16,32 and every 1<=d<k
agree with (10); the argument above proves it for all k,d.

## 5. Why a small coefficient vector is not yet an endpoint criterion

For the literal radial configuration, Q=N and H=R=N^(1/2).
This is outside the fourth-root regime of Section 3. Moreover,
the coefficient lattice has obvious exceptionally short vectors:
powers of the circle form

    F(X,Y)=(X^2+Y^2)^r

have small coefficients and satisfy

    F(x_i,y_i)=N^r,

which is exactly their forced divisor, not zero. Thus a generic
request for a shortest vector smaller than the covolume bound
could be satisfied by these forms without resolving anything.

A useful new criterion would have to control evaluations in an
arc-adapted norm, with the subspace generated by the circle
equation accounted for, or supply a conductor-preserving
transformation that actually reaches the needed height regime.
This note proves neither. The adaptive-multiplicity obstruction
does not rule out such a stronger arithmetic or geometric argument.

In particular, Sections 3 and 4 are a rigorous inverse-sieve
scope calculation, not a reduction of endpoint uniformity to
an already formulated shortest-vector theorem.

## 6. Available inverse-sieve theorems do not fill these gaps

Zywina's [larger sieve for rational points, Theorem 3.1](https://davidzywina.github.io/papers/Quantitative-HIT.pdf)
already handles projective denominators. With two residue classes
over Q, its decisive denominator is log Q/2-log(2H^2), giving
the same fourth-root threshold. Replacing rational messages by
bounded integers is unnecessary for that theorem.

Menconi, Paredes, and Sasyk's [projective inverse theorem, Theorem 1.6](https://arxiv.org/html/1907.02049)
requires a weighted dense set of small primes and 0<=k<d-1.
For a projective curve d=1, there is no eligible nonnegative k;
the paper explicitly leaves the boundary case open. Also,
logarithmic conductor weight does not imply the required density
weighted by log p/p. Its affine alternatives retain a positive-power
small-set exception, which does not control arbitrarily slowly
growing endpoint clusters. These are limitations of the cited
statements, not an assertion that no other method could work.

The new concrete conclusions are the complete loss of good
conductor under three-anchor normalization and the exact
no-saving calculation for adaptive projective interpolation.
Uniform boundedness remains open in this work.
