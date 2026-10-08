# Composite inert CRT collisions do not add valuation mass

This note tests whether simultaneous residue compatibility at products of
inert prime powers strengthens the collision lower bound in
`inert_prime_cofactor_bound.md`.  There is an exact accounting identity:
the available logarithmic divisibility mass is supported on prime powers.
Composite congruences are intersections of those events and carry no new
von Mangoldt weight.  An actual equal-norm Gaussian family realizes all
these intersections simultaneously with the same
`(1/4+o(1))M^2 log M` cost as the separate-prime bound.

This closes a residue-only composite-modulus enhancement.  It does not
exclude an argument that also forces new compatibility from the endpoint
arc scale or the near-uniform conductor profile.

## 1. Exact joint collision accounting

Use the primitive cofactor notation

```text
c_ij=(z_i-z_j)/gcd_G(z_i,z_j).
```

After removing the common Gaussian gcd, every rational prime
`p=3 mod 4` is absent from the common norm.  Therefore

```text
p^a | c_ij  iff  z_i=z_j mod p^a.                  (1)
```

For an integer `m` supported on inert primes, put

```text
C_m=#{i<j:z_i=z_j mod m}.
```

CRT gives the exact intersection identity

```text
1_(m|c_ij)=product_(p^a || m) 1_(p^a|c_ij).         (2)
```

Let `n_ij` be the inert rational-prime part of the largest rational
integer dividing both coordinates of `c_ij`.  Applying
`log n=sum_(d|n) Lambda(d)` pairwise gives

```text
sum_(i<j) log n_ij
 =sum_(m>=2, inert-supported) Lambda(m) C_m
 =sum_(p=3 mod 4)sum_(a>=1) C_(p^a) log p.          (3)
```

The middle sum is over all inert-supported integers, but `Lambda(m)=0`
when `m` contains two distinct primes.  Thus a collision modulo `pq` is
already the intersection of the collisions modulo `p` and `q`; assigning
it an additional weight `log(pq)` counts both existing valuations again.
For one modulus the inequality

```text
log(m) 1_(m|c_ij)
 <=sum_(p^a|m) log(p) 1_(p^a|c_ij)
```

is valid, but summing it over overlapping composite moduli is not valid
without allocating the prime-power weights between them.  Any weighted
composite-modulus larger-sieve argument must obey this same finite
valuation budget.

## 2. The CRT conic classes

If `gcd(m,N)=1`, the norm conic modulo
`m=product p^(a_p)` is the Cartesian product of its local fibers.  For
inert primes its exact number of classes is

```text
Q(m)=product_(p^a || m) (p+1)p^(a-1).               (4)
```

Hence the direct composite occupancy bound is

```text
C_m>=E(M,Q(m)),                                     (5)
E(M,q)=q binom(floor(M/q),2)+(M mod q)floor(M/q).
```

Formula (5) is compatible with the separate local bounds, but cannot be
added to them with a fresh `log m` weight.  It controls an intersection in
(2), not another divisor of the cofactor.

There is also no abstract obstruction to balanced joint CRT occupancy.  On
a finite collection of local fibers, repeat every element of their
Cartesian product equally often.  Every projection and every joint
projection is then balanced.  This relaxed observation alone would not be
enough for the circle problem, because arbitrary CRT classes need not lift
to small points on one integral circle.  The next section supplies actual
equal-norm Gaussian points instead.

## 3. An actual compatible family at every composite modulus

For `M>=2` and `n>=1`, set

```text
H_j=2(n+j)+i,
z_j=H_j product_(ell!=j) conjugate(H_ell),     0<=j<M.       (6)
```

All `z_j` have the same norm and satisfy every conic and Plucker identity.
Let `B=product_ell conjugate(H_ell)`.  Modulo any inert odd integer `m`,
all `H_j`, their conjugates, and `B` are units, and

```text
z_j=B H_j/conjugate(H_j),
H_i conjugate(H_j)-H_j conjugate(H_i)=4i(j-i).
```

It follows simultaneously for every inert-supported `m` that

```text
z_i=z_j mod m  iff  i=j mod m.                      (7)
```

In particular

```text
C_m=E(M,m).                                         (8)
```

The statement remains true for the primitive pair cofactors.  No inert
prime divides the common norm, so the Gaussian pair gcd is a unit at every
such prime and cannot change (7).

At a prime-power level the universal conic lower bound uses
`Q(p^a)=(p+1)p^(a-1)`, whereas (8) uses `p^a`.  Define

```text
F(M)=sum_(p=3 mod 4,a>=1) E(M,(p+1)p^(a-1)) log p,
V(M)=sum_(p=3 mod 4,a>=1) E(M,p^a) log p.
```

The exact identity `E(M,q)=sum_(b>=1)max(M-bq,0)` gives

```text
0<=V(M)-F(M)<=2M^2.                                 (9)
```

The proof is the bounded calculation already recorded in
`ordered_inert_sieve_route.md`: increasing `p^a` by the factor `1+1/p`
and summing the resulting triangular differences costs `O(M^2)`, uniformly
in `M`.  Since weighted Mertens gives

```text
F(M)=(1/4+o(1))M^2 log M,
```

the exact functions `F(M)` and `V(M)` have the same leading term as
`M` tends to infinity.  This is separate from the fixed-`M`, `n`-varying
circle construction and does not assert a joint uniform error term in
`(M,n)`.  The residue partitions of (6) at different primes and powers are
fully compatible by (7), including every composite intersection.

Consequently any lower bound valid for all actual equal-norm conic tuples
and based only on these joint congruence partitions and primitive cofactor
divisibility cannot force more than order `M^2 log M`, or improve its
leading coefficient through compatibility alone.  Combined with the
geometric upper bound of size `M log R`, it cannot improve the established
`log R/log log R` growth rate.

## 4. Scope of the countermodel

The family (6) is the actual circle construction already used in
`ordered_inert_sieve_route.md`, here reused to audit composite CRT
accounting rather than presented as a new circle countermodel.  For fixed
`M` and `n` tending to infinity, its angular span is asymptotic to
`(M-1)/n^2`, while its radius is asymptotic to `(2n)^M`.

There is also a uniform obstruction to using this family as a growing-point
endpoint example.  The angular gap between labels zero and one is exactly

```text
2 atan(2/(4n(n+1)+1)) >= 2/(9n^2).
```

Indeed `(2n+1)^2<=9n^2` and `atan(x)>=x/2` for `0<=x<=1`.  If `R`
denotes the common radius, then `sqrt(R)>=(2n)^(M/2)`.  Hence for
`M>=4` the normalized angular span satisfies

```text
Delta sqrt(R) >= (2/9) 2^(M/2) n^(M/2-2).           (10)
```

This diverges with `M`, uniformly for every `n>=1`.  Thus the family is not
a bounded `C sqrt(R)` endpoint-arc family when its number of points grows,
and its large conductor factors are mostly singleton cuts.  It does not
rule out a theorem coupling joint residues to endpoint geometry or to the
near-uniform cut profile.

What it does rule out is a gain obtained merely by replacing separate
prime-power occupancies with their simultaneous CRT occupancies.  Equation
(3) shows why such a calculation double-counts, and (6)--(9) show that the
compatible pattern already occurs on a genuine integral norm conic with
the primitive cofactor convention retained.

`check_composite_inert_crt_collision_budget.py` verifies (1)--(8), the
composite conic class count, the von Mangoldt accounting, and exact Gaussian
cofactors on bounded instances.
