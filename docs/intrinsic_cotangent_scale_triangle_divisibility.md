# The intrinsic cotangent scale divides triangle half-determinants

The minimal integer clearing scale for all half-angle cotangents is
intrinsic to the actual circle tuple. More strongly, each triangle
half-determinant is exactly the product of its three reduced edge
denominators, its Gaussian content norm, and a specified parity factor.
Multiplying gives a conductor-and-residue divisibility with a primewise
sharp residue exponent. This is exact Vandermonde accounting; its direct
endpoint consequence does not prove the desired exponent-four height bound.

## 1. Canonical clearing scale and edge divisibility

Let `z_0,...,z_(m-1)` be `m>=3` distinct Gaussian integers with common
positive norm `N`, and assume their full Gaussian gcd is a unit. Then `N`
is odd: if it were even, `1+i` would divide every row. Write `P_i` for
the corresponding integer coordinate vectors. For every pair define the
rational half-angle cotangent

```text
c_ij=(N+P_i dot P_j)/det(P_i,P_j)   if P_j != -P_i,
c_ij=0                            if P_j = -P_i.
```

The only zero central determinant for distinct rows occurs at antipodes,
so the separate second line covers all exceptional cases. Reversing the
edge changes the sign of `c_ij`. Let `s_ij` be its positive reduced
denominator, with the denominator of zero defined to be one. Put

```text
L_0=lcm_(i<j) s_ij.
```

This is the least positive integer for which all `L_0 c_ij` are integers.
It is unchanged by a common rotation, reflection, or Gaussian similarity
of the point configuration. Anchoring at any actual row and taking
`X_i=L_0 c_0i` gives the usual integer cotangent presentation. In particular
the finite pair quantities `(X_i X_j+L_0^2)/(X_i-X_j)` are integers;
they equal one of the signed values `L_0 c_ij`. Any integral clearing
scale for all these cotangents is a multiple of `L_0`.

For three distinct indices set

```text
D_ijk=det(P_j-P_i,P_k-P_i),
B_ijk=|P_k-P_i|^2-(P_j-P_i) dot (P_k-P_i).
```

The inscribed-angle identity, including its orientation, is

```text
D_ijk c_ij = B_ijk.                                  (1)
```

One direct verification uses `a=P_j-P_i`, `b=P_k-P_i`. Equal norms
give `2 P_i dot b=-|b|^2`. The vector `P_i+P_j` is perpendicular
to `a`; its coefficient relative to a quarter-turn of `a` is the
negative of `c_ij`. Taking its dot product with `b` yields (1).
At antipodes `P_i+P_j=0` and `B_ijk=0`, so (1) still holds.

Since `N` is odd, every row has coordinate sum odd. Thus both coordinates
of every difference have the same parity. It follows that `D_ijk` is even,
every squared difference norm is even, and the dot product of two
differences is even. Hence `B_ijk` is also even. All three distinct points
are noncollinear, because a line meets a circle in at most two points.
Dividing (1) by two and reducing the fraction proves

```text
s_ij | |D_ijk|/2             for every k distinct from i,j.  (2)
```

Equivalently, `s_ij` divides the gcd of all triangle half-determinants
containing the edge `ij`. This conclusion uses actual odd-norm integer
coordinates, rather than an arbitrary choice of cotangent clearing scale.

## 2. Exact triangle content and all three edge denominators

For an unordered triple `T={i,j,k}`, put

```text
a_T=|D_ijk|/2,       A_T=Norm(gcd_G(z_i,z_j,z_k)),
eta_T=2 if all three x-coordinates have the same parity, and 1 otherwise.
```

Then the exact identity is

```text
a_T = eta_T A_T s_ij s_ik s_jk.                        (3)
```

Here is a proof retaining both Gaussian orientations and the prime two.
Write each reduced edge cotangent as `a/s`, and put `epsilon=2` when
`a,s` are both odd, otherwise `epsilon=1`. Its primitive pair norm is

```text
n_ij=(a^2+s^2)/epsilon=N/Norm(gcd_G(z_i,z_j)).
```

The first equality follows by reducing `(a+is)/(a-is)` in `Z[i]`;
the only common factor of its coprime integer coordinates is the possible
ramified factor of norm two. The second is independent primitive Gaussian
normalization of the actual pair. These statements include antipodes,
where `a=0,s=1,n_ij=1`.

For the three pairs of a triple one has exactly

```text
n_ij n_ik n_jk=(N/A_T)^2.
```

For verification, first remove the triple Gaussian gcd. The remaining
primitive norm has only odd split prime factors. At any such prime its
three allocations have minimum zero and maximum `e`; the pair norm
exponents are their absolute differences, whose sum is `2e`. This proves
the identity with full prime-power multiplicities, and restores `A_T`
without assuming pairwise coprimality of the original norms.

The chord formula is

```text
|z_i-z_j|^2=4N s_ij^2/(epsilon_ij n_ij).
```

Since `N,n_ij` are odd, its two-adic valuation is one exactly when
`epsilon_ij=2`. This is also exactly when the endpoint x-parities differ:
the chord coordinates are both odd in that case and both even otherwise.
Among three vertices, either zero or two edges cross the two parity
classes. Thus the product of their `epsilon_ij` is one if all three
classes agree, and four otherwise. Finally the circumradius formula gives

```text
a_T^2=4 A_T^2 (s_ij s_ik s_jk)^2
                  /(epsilon_ij epsilon_ik epsilon_jk).
```

Taking positive square roots proves (3). In particular the product of
all three edge denominators, not merely their lcm, divides `a_T`.

The link to the existing primitive chord residues is also exact. Put
`c_ij=(z_i-z_j)/gcd_G(z_i,z_j)`. The chord formula gives

```text
|c_ij|^2=4 s_ij^2/epsilon_ij.
```

With `S=product_(i<j) s_ij`, if the two coordinate-parity classes
have sizes `m_0,m_1`, then

```text
product_(i<j) |c_ij| = 2^(binom(m,2)-m_0 m_1/2) S.
```

Thus the reduced cotangent denominators retain precisely the existing
primitive chord-residue sizes, with the explicit dyadic correction.
They are not a second independent source of residue savings; compare
the [inert-prime cofactor bound](inert_prime_cofactor_bound.md).

## 3. Global conductor and residue divisibility

Put

```text
S=product_(i<j) s_ij,
b_m=binom(floor(m/2),3)+binom(ceil(m/2),3),
b_T=lcm(s_ij,s_ik,s_jk) for T={i,j,k}.
```

A binomial coefficient with top entry less than three is zero.
Multiplying (3) gives the exact identity

```text
product_T a_T = (product_T eta_T A_T) S^(m-2).
```

Full tuple primitivity supplies the additional source conductor factor:

```text
(2N)^b_m S^(m-2) | product_T a_T.                        (4)
```

At an odd split source prime `p`, let `e=v_p(N)` and let `t_i` be
its allocations, with minimum zero and maximum `e`. The exponent of
`A_T` is `e-range(t_i:i in T)`. Expressing the allocations as the
`e` threshold cuts, this counts exactly the layers constant on the
triple. A cut with `r` rows on one side has
`binom(r,3)+binom(m-r,3)` monochromatic triples. This convex symmetric
expression is minimized by the balanced sizes and is at least `b_m`.
Summing over layers proves `N^b_m | product_T A_T`. There are no inert
or ramified factors in the full primitive norm. Independently,
`product_T eta_T` is two to the number of monochromatic triples in the
actual two parity classes, so it is divisible by `2^b_m` by the same
count. Combining these factors with the exact product identity proves
(4), without separating shared primes between `N` and `S`.

There is also the intrinsic-scale divisibility chain

```text
L_0^(m-2) | product_T b_T | product_T a_T.             (5)
```

For the first divisibility, at each rational prime choose an edge whose
denominator attains `v_p(L_0)`. Each of its `m-2` containing triples
has `v_p(b_T)>=v_p(L_0)`. The second follows from `b_T|a_T`.
Thus (5) retains the local lcms as well as the stronger residue-product
information in (4).

## 4. The exponent is primewise sharp for actual primitive circles

Fix `m>=3`, an odd inert prime `p>m` (`p=3 mod 4`), an integer `h>=1`,
and a positive integer `H` coprime to `p`. Take half-angle parameters

```text
t_0=0, t_1=p^h, t_2=1, ..., t_(m-1)=m-2,
q_j=(H+i t_j)/(H-i t_j).
```

These are distinct rational circle phases. Clear their reduced Gaussian
denominators and remove the full tuple gcd to obtain actual primitive
Gaussian integers `z_j`, of common norm `N`. Since `p` is inert and
`H` is a unit modulo `p`, every `H+i t_j` and `H-i t_j` is a unit at
the Gaussian prime `p`. Thus `p` does not divide `N`.

The edge cotangents are exactly

```text
c_ij=(H^2+t_i t_j)/(H(t_j-t_i)).
```

Only the difference `t_1-t_0` is divisible by `p`; all other differences
are units modulo `p`. The numerator on that exceptional edge is `H^2`,
also a unit. Consequently

```text
v_p(s_01)=h,      v_p(s_ij)=0 for every other edge,
v_p(L_0)=h.                                           (6)
```

The corresponding phase difference is
`q_j-q_i=2iH(t_j-t_i)/((H-i t_i)(H-i t_j))`. All denominator factors
and the physical anchor `z_0` are units at `p`. Hence the physical chord
valuations are `h` on edge `01` and zero on every other edge. The exact
circle triangle identity

```text
2i D_ijk = +/- N (z_j-z_i)(z_k-z_i)(z_j-z_k)/(z_i z_j z_k)
```

then gives

```text
v_p(a_T)=v_p(b_T)=h if {0,1} is contained in T,
v_p(a_T)=v_p(b_T)=0 otherwise.
```

Thus both products in (5) have valuation exactly `(m-2)h` at this prime.
The exponent `m-2` cannot be increased in a universal divisibility of
this form, even for actual primitive circle tuples with all cotangents
defined. This is primewise sharpness, not a claim that equality holds
at the other primes.

The phases lie on the minor arc of width `2 arctan(p^h/H)`. Taking `H`
arbitrarily large while preserving `p`-coprimality makes this angular
width small. There is **no assertion** that their normalized endpoint
constant stays bounded: their least radius changes with `H`. The examples
therefore do not rule out an endpoint-restricted improvement using extra
geometry, and are not counterexamples to the desired height bound.

## 5. Endpoint consequence and its limitation

Suppose the actual tuple lies on a minor arc of angular width
`Delta<=C N^(-1/4)`. For any three ordered rows with successive angular
gaps `a,b`, the chord-product formula and `sin(t)<=t` give

```text
|D_ijk|=4N sin(a/2)sin(b/2)sin((a+b)/2)
       <= N a b (a+b)/2 <= N Delta^3/8.
```

Hence every `a_T<=C^3 N^(1/4)/16`. Equation (4) yields the stronger
product bound

```text
S <= 2^(-b_m/(m-2)) (C^3/16)^(m(m-1)/6)
       * N^(floor(m/2)/4).             (7)
```

Indeed the radius exponent is
`[binom(m,3)/4-b_m]/(m-2)=floor(m/2)/4`, checked separately for even
and odd `m`. Since `L_0|S`, (7) also bounds the intrinsic scale, rather
than an artificially enlarged presentation scale. For five points its
radius exponent is `1/2`; the weaker chain (5) alone would give `5/6`.
Neither removes the genuine radius dependence in the existing local
cotangent capacities. Equation (4) is precisely a triangle-product
version of the existing source-cut/Vandermonde accounting; it does not
supply a new endpoint exponent-four height inequality.

## Verification

The [exact checker](check_intrinsic_cotangent_scale_triangle_divisibility.py)
checks the oriented identity, exact triangle content, parity, edge gcds,
and the global divisibilities (4)--(5) on actual circles, including
antipodal pairs. It separately
constructs and fully normalizes the sharpness family, compares its radius
with an ordinary all-edge denominator lcm, and checks the stated prime
valuations. The all-edge comparison never substitutes an anchor-only
ordinary lcm for the true primitive radius. Finite checks supplement the
proofs above. Astra derived the intrinsic-scale divisibility, root
strengthened it to the exact content identity and supplied the sharpness
construction, and Luna independently audited the triangle identity,
parity, source-cut product and endpoint constant. Root also checked the
exact identity on 6,064 additional literal circle triples.
