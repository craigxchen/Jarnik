# Angular order and the joint inert-prime sieve

This note tests whether retaining angular order and simultaneous
congruences can improve the growth rate in
[inert_prime_cofactor_bound.md](inert_prime_cofactor_bound.md).
There is an exact sharpness theorem: genuine ordered Gaussian circle
points can attain the separate inert-prime collision lower bounds
simultaneously, up to `O(M^2)` in the total logarithmic cost. The
construction satisfies the actual conic and every Plucker identity.
Its conductor profile is asymptotically singleton, so it does not
settle the question with near-uniform cut weights. We also compute
the optimal conductor budget for nonnegative pair reweightings and
an explicit angular-order upper bound. None changes the established
`log R/log log R` growth rate.

## 1. Simultaneous near-minimal collisions on an ordered actual circle

For integers `M>=2,n>=1`, set

```text
H_j=2(n+j)+i,                  0<=j<M,
z_j=H_j product_(l!=j) bar(H_l),
R_n=product_j |H_j|.
```

These are distinct Gaussian points of common modulus `R_n`.
Their arguments differ by twice the arguments of the `H_j`, so
their angular order is exactly the reverse of the index order.
This is the genuine family already used for fixed-size residue
examples in [residue_prime_incidence.md](residue_prime_incidence.md).
The following statement is uniform in both `M` and `n`.

For every inert rational prime `p=3 mod 4` and every `a>=1`,

```text
z_i=z_j mod p^a  if and only if  i=j mod p^a.       (1)
```

Indeed all `H_j` and their conjugates are units modulo `p^a`.
Writing `B=product_l bar(H_l)`, one has `z_j=B H_j/bar(H_j)`.
After canceling units, the numerator of a pair difference is

```text
H_i bar(H_j)-H_j bar(H_i)=4i(j-i).
```

Since `p` is odd, this proves (1) to every depth. It follows that
the number of colliding pairs at depth `a` is exactly

```text
E(M,p^a),
E(M,q)=q binom(floor(M/q),2)+(M mod q)floor(M/q).  (2)
```

Thus the occupancies are balanced simultaneously at every inert
prime power, and those residue classes occur periodically in the
geometric ordering. Any common Gaussian divisor of the whole tuple
is a unit at these primes, so primitive normalization preserves (1).

Define the actual inert part of the primitive cofactor product by

```text
V(M)=sum_(p=3 mod 4)sum_(a>=1) E(M,p^a) log p.
```

The universal conic-occupancy bound uses `(p+1)p^(a-1)` classes:

```text
F(M)=sum_(p=3 mod 4)sum_(a>=1)
       E(M,(p+1)p^(a-1)) log p.
```

Then the following entirely elementary bound holds for every `M>=2`:

```text
0 <= V(M)-F(M) <= 2M^2.                          (3)
```

To prove it, use the exact identity

```text
E(M,q)=sum_(b>=1) max(M-bq,0).
```

For `q=p^a`, increasing `q` to `q+q/p` changes only the positive
summands indexed by `b<M/q`. Their total decrease is at most

```text
(q/p)sum_(b<=M/q)b <= M^2/(2pq)+M/(2p).
```

Summing over depths and primes therefore gives

```text
V(M)-F(M)
 <= (M^2/2)sum_(p=3 mod 4) log p/[p(p-1)]
    +(M log M/2)sum_(p<=M, p=3 mod 4)1/p
 <= (M^2/2)sum_(r>=3) log r/[r(r-1)]
    +(M/2)(log M)^2
 < 2M^2.                                        (4)
```

For the final bound, the displayed convergent series is less than
two by comparison with `2 sum_(r>=3)log(r)/r^2`, and
`(log M)^2/M<=4/e^2`. No prime-distribution theorem is used in
(3). The weighted Mertens estimate used in the inert-prime note
gives `F(M)=(1/4+o(1))M^2 log M`: its matching upper bound follows
from `E(M,q)<=M^2/(2q)`, with the depths beyond the first contributing
`O(M^2)`. Hence the simultaneous cost has the same leading term.

This rules out a stronger leading collision cost based only on
angular order, all inert-prime congruence partitions at all depths,
and the algebraic equations of actual circle points. Those
hypotheses all hold in this family. It does not rule out an argument
that additionally uses the leading conductor weights or the
endpoint size of the residues.

### The radius and conductor costs remain essential

When `M` is fixed and `n` tends to infinity, the angular span is

```text
Delta_n=2[arctan(1/(2n))-arctan(1/(2(n+M-1)))]
       ~ (M-1)/n^2,
R_n~(2n)^M.
```

Hence `Delta_n sqrt(R_n)` grows like a positive constant times
`n^(M/2-2)` for `M>=5`. The bounded common-factor removal described
in the earlier note does not change that exponent. This is not an
endpoint counterexample.

The cut-profile failure is equally explicit. A Gaussian prime
occurring in only one factor among the `H_j` and their conjugates
creates a singleton cut. Repeated factors divide one of

```text
H_i-H_j=2(i-j),
H_i-bar(H_j)=2(i-j)+2i.
```

The sum of all pair-overlap logarithmic norm bounds is
`O(M^2 log(M+1))`. Thus, if `log n` is chosen much larger than
`M log(M+1)`, the conductor has singleton cuts outside a negligible
fraction of its total height `2M log n+O(M)`. The example specifically
does not combine near-minimal inert collisions with the extracted
near-uniform conductor distribution.

## 2. The optimal cut budget for nonnegative pair weights

Let `w_ij>=0` be arbitrary weights on pairs of `M` labels, allowed
to depend on their angular ordering or observed residue classes.
Put

```text
A=sum_(i<j) w_ij,
C_max=max_S sum_(i in S, j outside S) w_ij.
```

A uniformly chosen balanced cut separates a fixed pair with
probability `M/[2(M-1)]` when `M` is even, and `(M+1)/(2M)` when
`M` is odd. Averaging therefore proves

```text
C_max-A/2 >= A/[2(M-1)]    (M even),
C_max-A/2 >= A/(2M)        (M odd).               (5)
```

Equal weights on all pairs attain equality. This is an exact
finite optimization statement, not an assumption that actual
conductor cuts are uniformly balanced.

For primitive Gaussian chord cofactors, the pair identity in the
inert-prime note is

```text
log|c_ij| <= log(|z_i-z_j|/sqrt(R))-W/4+d_ij/2,
W=log(R^2).
```

Consequently the usual cut upper bound for any such weighting is

```text
sum_(i<j) w_ij log|c_ij|
 <=sum_(i<j) w_ij log(|z_i-z_j|/sqrt(R))
    +(W/2)(C_max-A/2).                           (6)
```

The coefficient of `W` per unit pair weight in this procedure
cannot improve on the equal-weight complete graph, by (5).
Selecting consecutive pairs or pairs colliding at chosen primes
does not circumvent that fact. Such a selection could still be
useful if it supplied a stronger arithmetic lower bound or a new
constraint on the actual conductor cuts; (5) concerns only the
nonnegative weighted cut maximum used in (6).

## 3. Retaining angular order in the complete product

There is a clean upper bound for the complete geometric product.
Let the arc length be `L<=C sqrt(R)` and normalize the ordered
arc parameters to `s_i in [0,1]`. Chord length is at most arc
separation, so

```text
product_(i<j)|z_i-z_j| <= L^K product_(i<j)|s_i-s_j|,
K=binom(M,2).
```

For `j>=1`, the monic polynomial
`2^(1-2j) T_j(2s-1)` has supremum at most `2^(1-2j)` on `[0,1]`;
here `T_j` is the elementary Chebyshev polynomial. Replace the
monomial columns of the Vandermonde determinant by these monic
polynomial columns. The determinant is unchanged, and the column
norm form of Hadamard's inequality gives

```text
product_(i<j)|s_i-s_j|
 <= M^(M/2) 2^[-(M-1)^2].                        (7)
```

Combining (7) with the same conductor product bound gives

```text
F(M) <= b_M log R+K log C
        -(M-1)^2 log 2+(M/2)log M,               (8)
b_M=M/4 (M even),      b_M=(M-1)/4 (M odd).
```

Thus angular ordering supplies a real improvement in the additive
geometric term. Its order is `M^2`, while the inert-prime term has
order `M^2 log M` and the conductor allowance has order `M log R`.
It does not improve the growth rate. The near-minimal simultaneous
collision theorem (3) explains why adding compatibility of those
same inert partitions does not automatically create the missing
extra factor.

## 4. What remains after this test

The tested new route combines three pieces of information that
could, in principle, have interacted: angular order, simultaneous
prime-power collisions, and conic/Plucker compatibility. The first
two quantitative conclusions are (3) and (5), and the explicit
combined estimate is (8). They leave a precise remaining issue:
whether the **near-uniform conductor weights** force the actual
ordered inert collision pattern, or the primitive residue values,
to be much more expensive than in the singleton-profile family.
No such implication is proved here. A bound uniform in `R` remains
unproved.

Exact arithmetic checked (1) and (2) in 288 ordered Gaussian
configurations, making 374,400 pair/depth comparisons at inert
primes `3,7,11,19` and depths one through three. A separate finite
evaluation of (3) through `M=2000` had maximum
`(V(M)-F(M))/M^2` approximately `0.114105`. These checks supplement
the unrestricted proofs; they are not a formalization.
