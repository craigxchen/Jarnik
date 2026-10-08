# Inert-prime collisions and a bound for all unit classes together

The new arithmetic input in this note concerns prime divisibility of
the actual primitive Gaussian chord quotients. It gives a sharper
general bound,

```text
M <= (1+epsilon) log R/log log R
```

for each fixed `C>0,epsilon>0` and all sufficiently large `R`, uniformly
in the position of an arc of length at most `C sqrt(R)`. This improves
the previous prose leading constant four, but **does not improve the
growth rate and does not complete the uniform-bound objective**.

Unlike the common-unit primitive imaginary-part formula, this proof
handles all Gaussian unit classes at once. It is a prose proof; no
Lean formalization of its prime-distribution input is claimed.

## 1. Remove only the actual common Gaussian factor

Let `z_1,...,z_M` be distinct Gaussian integers of common modulus `R`.
If their common Gaussian gcd is `g`, divide the whole configuration
by `g`. The new radius and containing arc length satisfy

```text
R'=R/|g| <= R,
L'=L/|g| <= (C/sqrt(|g|)) sqrt(R') <= C sqrt(R').
```

Thus an upper bound in terms of `R'` and `C` implies the corresponding
bound in terms of `R`. In the normalized configuration the common gcd
is a unit. Its norm `N=R^2` has no inert prime factor and is odd:
an inert prime or the ramified Gaussian prime would otherwise divide
every point to the same positive order. Hence

```text
N=product_(p=1 mod 4) p^(e_p),
z_j=epsilon_j product_p pi_p^(a_jp) conjugate(pi_p)^(e_p-a_jp),
0<=a_jp<=e_p.
```

There is no restriction on the units `epsilon_j`.

For each pair choose any Gaussian gcd `g_ij` and put

```text
c_ij=(z_i-z_j)/g_ij in Z[i] minus {0},
P=product_(i<j)|c_ij|,
K=binom(M,2).
```

Changing the associate of a gcd changes no modulus or divisibility
used below.

## 2. Geometry and conductor cuts give the upper bound

Define the weighted allocation distance

```text
d_ij=sum_p |a_ip-a_jp| log p.
```

Unique factorization gives

```text
log|g_ij|=log R-d_ij/2.
```

The chord bound `|z_i-z_j|<=C sqrt(R)` therefore implies

```text
log|c_ij| <= log C-(log R)/2+d_ij/2.                (1)
```

At each threshold of one prime allocation, at most
`floor(M^2/4)` pairs lie on opposite sides. Summing over all
thresholds gives

```text
sum_(i<j) d_ij <= floor(M^2/4) log N.
```

Consequently

```text
log P <= b_M log R+K log C,                         (2)
b_M=floor(M^2/4)-M(M-1)/4
   =M/4       if M is even,
   =(M-1)/4   if M is odd.
```

This step needs neither common point units nor a primitive imaginary
part; the nonzero Gaussian quotient itself is the integral object.

## 3. Exact finite divisibility at inert primes

Fix a rational prime `p=3 mod 4`. Since `p` does not divide `N`, it
does not divide any `g_ij`. Therefore for every `a>=1`,

```text
p^a divides c_ij in Z[i]
    iff z_i=z_j mod p^a.                            (3)
```

The conic `x^2+y^2=N` has exactly

```text
q_(p,a)=(p+1)p^(a-1)
```

solutions modulo `p^a`. At level one this is the norm fiber in
`F_(p^2)^*`, of size `p+1`; equivalently it follows by parametrizing
the nonsingular conic from the reduction of one of the given points.
At every higher level its nonzero gradient gives exactly `p` lifts
per point.

For `M` objects placed in at most `q` classes, the minimum number of
equal-class unordered pairs is

```text
E(M,q)=q binom(b,2)+r b,
b=floor(M/q),       r=M-qb.                         (4)
```

Indeed, if two occupancies differ by two or more, transferring one
object from the larger to the smaller strictly decreases their pair
count. Thus the minimum has `r` occupancies `b+1` and the others `b`.

Applying (4) at every prime-power depth in (3) yields the entirely
finite arithmetic function

```text
F(M)=sum_(p=3 mod 4) sum_(a>=1)
       E(M,(p+1)p^(a-1)) log p
```

and the lower bound

```text
log P >= F(M).                                      (5)
```

Only terms with `(p+1)p^(a-1)<M` can be nonzero. Divisibility by
an inert rational prime costs modulus `p`, not `sqrt(p)`, so the
weight in (5) is `log p`. Combining with (2), every endpoint
configuration satisfies

```text
F(M) <= b_M log R+binom(M,2) log C.                 (6)
```

The original, unnormalized radius can also be used on the right of
(6), by Section 1. This is a finite, computable bound without an
asymptotic prime theorem.

## 4. Consequence for the general asymptotic bound

The first-depth contribution already suffices. Cauchy--Schwarz on
the occupancies, or (4), gives

```text
E(M,p+1) >= M^2/(2(p+1))-M/2.
```

Sum over inert primes `p<=M`. The unconditional weighted Mertens
estimate in the residue class `3 mod 4` gives

```text
sum_(p<=M, p=3 mod 4) log p/(p+1)
    =(1/2+o(1)) log M.                              (7)
```

Replacing `p+1` by `p` changes a bounded sum. The weighted
prime-in-progression asymptotic is recorded in Section 5 of
Kenneth S. Williams,
[Mertens' Theorem for Arithmetic Progressions](https://people.math.carleton.ca/~williams/papers/pdf/057.pdf),
Journal of Number Theory 6 (1974), 353--359. It also follows by
partial summation from the prime number theorem in this fixed
progression. The elementary upper bound
`sum_(p<=M) log p=O(M)` controls the subtracted term.

Thus

```text
F(M) >= (1/4+o(1)) M^2 log M.                       (8)
```

Using `b_M<=M/4` in (6), for fixed `C` and `M` tending to infinity,

```text
(1+o(1)) M log M <= log R+O_C(M).                  (9)
```

For clarity, the `o(1)` on the left need not be positive. The precise
consequence is that for every `epsilon>0` there is a threshold
depending only on `C,epsilon` above which

```text
M <= (1+epsilon) log R/log log R.                  (10)
```

This follows by monotonicity of `x log x`, or by substituting
`(1+epsilon) log R/log log R` into (9); bounded `M` is immediate.
No factor four for the unit classes is lost anywhere in the proof.

Two independent audits checked the Gaussian normalization, the exact
finite inequality, and this asymptotic deduction. Exact calculations
on 1,035 configurations on circles of integer norm at most 1,000,
including mixed unit classes, verified the gcd/cofactor identities
and 7,859 outside-prime occupancy checks. These checks supplement
the proofs, and do not constitute a Lean formalization.

## 5. What this does and does not provide for uniformity

The primes outside the conductor impose a genuine residue cost of
order `M^2 log M`. At a fixed extracted tuple size, however, that
cost is independent of the radius. It does not imply
`log P >= delta log R` for any fixed positive `delta`.

In particular, configurations whose cardinalities grow arbitrarily
slowly in `R` are not excluded by (6). To reach the requested bound
one must obtain additional information about the actual residue
values, their interaction with the conductor, or some other
structure. The constant-one corollary is an incidental improvement,
not a redefinition of the active goal.
