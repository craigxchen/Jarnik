# Reciprocal chords: exact denominator content and a positive integer

This note tests rational Cauchy and squared-Cauchy identities on eight
actual ordered Gaussian circle points. It gives an exact positive-integer
congruence that detects one balanced cut, including arbitrary prime powers.
It does not prove the desired inequality `D+W_4<=2W+B`.

The useful new quantity is a primitive positive permanent numerator `N`.
For the balanced conductor factor `X` separating two ordered groups of
four, it satisfies

```text
gcd(N,X)=1,       N>=X+24.                               (1)
```

The corresponding positive rational ratio `N/K` has reduced denominator
divisible by `X^6`. Its real size can nevertheless be subpower in the
critical full-profile regime. The exact identities below specify the
remaining arithmetic question; a bounded real value is not a height bound.

## 1. Actual chords and a Gaussian least common denominator

Use the common-unit normalization of
[ordered_residue_growth.md](ordered_residue_growth.md):

```text
z_i=epsilon product_p pi_p^(a_i(p)) bar(pi_p)^(e_p-a_i(p)),
R^2=product_p p^e_p,       W=log R^2.
```

Label the eight distinct points in their cyclic order, choose a cut of
that order, and put `I={1,2,3,4}`, `J={5,6,7,8}`. For an oriented pair,
write

```text
ell_ij=(z_i-z_j)/(2i epsilon)=g_ij tau_ij,
g_ij=product_p pi_p^min(a_i,a_j) bar(pi_p)^(e_p-max(a_i,a_j)),
tau_ij in Z\{0},       t_ij=|tau_ij|.
```

The `g_ij` are the actual Gaussian gcd factors, not independent formal
weights. The normalization removes the common factor `2i epsilon` from
every chord. Let `C` be the `4x4` matrix with entries `1/ell_ij` for
`i in I,j in J`. For a permutation matching `sigma`, put

```text
M_sigma=product_(i in I) ell_(i,sigma(i)),
G_sigma=product_(i in I) g_(i,sigma(i)).
```

Let `L` and `L_0` be Gaussian least common multiples of the `M_sigma`
and `G_sigma`, respectively. All `M_sigma` have the same complex phase:
each matching uses all eight point angles once, and all four cross chords
have the same sign of their half-angle sine. Choose the unit of `L` so

```text
A_sigma=L/M_sigma in Z_(>0).                            (2)
```

Their integer gcd is one, by the defining least-common-multiple property.
Define

```text
N=sum_sigma A_sigma=L per(C)>0,
K=|sum_sigma sign(sigma) A_sigma|=|L det(C)|>0.           (3)
```

Nonvanishing of the determinant follows from the Cauchy formula and
distinctness of all eight points. These are actual positive integers.

## 2. The exact conductor part of `L`

For a threshold layer at `p`, let `a` and `b` count its marked vertices
in `I` and `J`. At the first Gaussian orientation, the maximum valuation
of a core matching product is

```text
sum_layers min(a,b).
```

This maximum is simultaneous over every layer at that prime: sort the
four allocations in each group and match corresponding ranks. For every
threshold, the terminal marked ranks then overlap in exactly `min(a,b)`
positions. At the conjugate orientation the same matching attains
`sum_layers min(4-a,4-b)`. Consequently

```text
log|L_0|=sum_layers (2-|a-b|/2) log p.                   (4)
```

This proves the formula for arbitrary allocation depths, not just endpoint
allocations. Adding a rational residue to an edge raises its two Gaussian
valuations equally. Moreover the difference between the two orientation
orders of a matching product is independent of the matching. It follows
that, up to a unit,

```text
L=lambda L_0,       lambda in Z_(>0),
lambda | T_cross,   T_cross=product_(i in I,j in J)t_ij. (5)
```

One can also see divisibility directly: every core matching divides `L_0`,
so every actual matching divides `L_0 T_cross`, while the coordinatewise
maxima for `L` are at least those for `L_0`. Inert and ramified primes of
the residues are retained in `lambda`; no assumption on their support is
used.

## 3. The Cauchy determinant has an exact residue quotient

Put

```text
T_within=product_(i<j in I)t_ij product_(i<j in J)t_ij,
G=product_(p,layers) p^binom(|a-b|,2).
```

The ordinary Cauchy determinant formula, unchanged by the common chord
normalization, is

```text
det(C)=+/- [product_(within I)ell_ij product_(within J)ell_ij]
            /product_(cross I,J)ell_ij.
```

At one layer its first-orientation core order is
`binom(a,2)+binom(b,2)-ab=((a-b)^2-a-b)/2`.
Adding the least-common-denominator order `min(a,b)` gives
`binom(|a-b|,2)`. The conjugate orientation gives the same answer.
Therefore the exact integer identity is

```text
K=G lambda T_within/T_cross
 =G T_within/B,       B=T_cross/lambda in Z_(>0).        (6)
```

In particular it would be incorrect to assert unconditionally that `G`
divides `K`: the denominator `B` must be retained. Its entire height is
charged to actual primitive residues. Equation (6) gives the precise
denominator-content reduction, including the integrality constraint
`B/gcd(B,T_within) | G`.

After this reduction the logarithmic-modulus coefficient
`log|L_0|-log G` obtained from the two Gaussian orientations is exactly

```text
2-(a-b)^2/2.                                           (7)
```

For a fixed eight-row cut of size `r`, averaging over all `4+4` partitions
gives `E(a-b)^2=r(8-r)/7`, so the average of (7) is
`2-r(8-r)/14=((r-4)^2+12)/14`.
It is affine in the usual quadratic defect; no term supported on balanced cuts
survives this denominator cancellation. This is an exact calculation for
this rational identity, rather than an appeal to a general polynomial
barrier.

## 4. The permanent detects the separated balanced factor

For each conductor prime define

```text
m_p=max(0, min_(j in J)a_j-max_(i in I)a_i,
           min_(i in I)a_i-max_(j in J)a_j),
X=product_p p^m_p.                                     (8)
```

Thus `log X` is exactly the total weight of the `4|4` threshold layers
that separate these two groups. Suppose, for example, that the first
minimum in (8) gives a positive gap `m_p`. At `pi`, each cross chord
is its `I` endpoint times a unit congruent to one modulo `pi^m_p`.
At `bar(pi)`, each cross chord is minus its `J` endpoint times such a
unit. Matching products use every endpoint once. Their ratios therefore
satisfy

```text
M_sigma/M_rho == 1 mod pi^m_p and mod bar(pi)^m_p.
```

These ratios are rational `p`-adic units, so the congruence is modulo
`p^m_p`. Every cross primitive residue is a `p`-unit, since the two
allocations are unequal. Thus `lambda` and all the integers `A_sigma`
are also `p`-units. The reverse allocation order has the identical proof.
For every choice of reference matching `rho`,

```text
A_sigma == A_rho mod X,       gcd(A_sigma,X)=1,
N == 24 A_rho mod X.                                   (9)
```

All conductor primes are split odd primes, hence at least five. None
divides `24`, proving `gcd(N,X)=1`.

The twenty-four matching products cannot all coincide. Swapping two
partners changes their product by a nonzero product of two within-group
chords, by the two-by-two Cauchy numerator identity. Hence two distinct
positive integers in (2) differ by at least `X`. Since all twenty-four
are at least one, this proves `N>=X+24`, including the case `X=1`.

At a prime dividing `X`, every cross residue is a unit, so (6) gives
`v_p(K)>=v_p(G)>=6m_p`. Combining with (9) proves

```text
X^6 divides the reduced denominator of N/K.             (10)
```

This is a concrete rational invariant with a controlled sign and a large
mandatory denominator. Its positive real value is not constrained to
have small rational height.

## 5. Squared reciprocals and the exact remaining height budget

Borchardt's identity says

```text
det(C entrywise squared)=det(C) per(C).
```

It applies to these scaled differences by homogeneity; see the primary
proof in [Dan Singer, A Bijective Proof of Borchardt's Identity, Theorem 1.1](https://www.kurims.kyoto-u.ac.jp/EMIS/journals/EJC/Volume_11/PDF/v11i1r48.pdf).
Clearing denominators yields exactly the product of the signed integer
in (3) and `N`. Thus this identity introduces the permanent numerator;
it does not identify that numerator with a product of old primitive
residues.

The scale is explicit in a critical fair eight-row profile. Here `W`
is inherited conductor weight, including any constant layers. Suppose
each of the 256 oriented patterns has weight `W/256+o(W)`, and all
primitive chord residues have logarithm `o(W)`. Every chord `ell_ij`
then has logarithmic modulus `W/4+o(W)`. For independent
`a,b~Bin(4,1/2)`, one has

```text
E|a-b|=35/32,       E(a-b)^2=2,
E(2-|a-b|/2)=93/64,
E binom(|a-b|,2)=29/64.
```

Equations (3)--(6), and positivity of the permanent terms, consequently
give

```text
log N=(29/64)W+o(W),
log K=log G+o(W)=(29/64)W+o(W),
log X=W/128+o(W),       log(N/K)=o(W).                  (11)
```

Thus the new lower bound `N>=X+24` is compatible with the available
height budget by a factor of 58 in its leading logarithmic coefficient.
The forced sixth-power denominator in (10) is also compatible with the
subpower real size in (11). Taking products over all 35 partitions does
not fix this mismatch: the corresponding 35 numerator heights are still
charged, whereas each balanced cut is counted once in the product of
the `X` factors.

An additional arithmetic estimate on these positive rational quantities,
for example a height bound in terms of the old residue product, would
be new input. Neither the Cauchy identity, Borchardt's identity, positivity,
nor the least-common-denominator calculation supplies such an estimate.
The conclusion is scoped to this route: it is not an arithmetic
counterexample to the eight-point target or a general impossibility proof.

## Exact verification

[check_reciprocal_chord_content.py](check_reciprocal_chord_content.py)
checks all identities on 80 actual equal-norm Gaussian eight-tuples,
including 60 separating-prime cases and 32 gaps larger than one. It checks
the Gaussian LCD, the exact residue quotient (6), the cleared Borchardt
identity, positivity and primitive gcd of the summands, congruence (9),
coprimality, inequality (1), and the sixth-power denominator (10).
Those examples are not asserted to lie in endpoint arcs. The arguments
above establish the identities and inequalities for arbitrary prime
powers; the finite checks are independent arithmetic validation.
