# Exact changes of anchor in the centrally truncated system

This note uses the notation and quantitative extraction in
[endpoint_central_truncation.md](endpoint_central_truncation.md).
It proves that the small imaginary residues at all anchors follow from
the one-anchor estimates once the extracted pair distances are retained.
It also isolates an additional condition that is not redundant: exact
divisibility by every pair's proposed core product. A concrete family
shows that this condition does not follow even from primitive, small-value
one-anchor equations with asymptotically uniform independent blocks.
However, the repair lemma in Section 7 shows that trimming a negligible
amount of core weight always restores this exact divisibility at the
quantitative coefficient scale supplied by central truncation. Thus this
distinction does not supply a new leading-height obstruction.
No uniform endpoint bound is proved.

## 1. The exact transition formula

Keep the original common Gaussian unit class and anchor `0`. Write

```text
h_i=product_p pi_p^((a_i-a_0)_+) bar(pi_p)^((a_0-a_i)_+),
h_ij=product_p pi_p^((a_j-a_i)_+) bar(pi_p)^((a_i-a_j)_+).
```

Thus `h_0=1`, `h_ji=bar(h_ij)`, and
`z_j/z_i=h_ij/bar(h_ij)`. Let `g_ij=gcd_G(h_i,h_j)`, choosing
its generator to be the literal product of the common prime powers.
Then

```text
h_j bar(h_i)=N(g_ij) h_ij.                       (1)
```

To verify this, put `x=a_i-a_0` and `y=a_j-a_0` at one rational
prime. The two conjugate exponents in `h_j bar(h_i)` are
`y_++(-x)_+` and `(-y)_++x_+`. Their minimum is the common
exponent of `h_i,h_j` when `x,y` have the same sign, and is zero
otherwise. Removing that equal pair of conjugate exponents leaves
`(y-x)_+` and `(x-y)_+`, as required. There are no Gaussian unit
ambiguities with the displayed generators.

In particular, writing `h_i=x_i+i t_i`, we get the exact integer
identity

```text
N(g_ij) t_ij=x_i t_j-x_j t_i,
t_ij=Im h_ij.                                    (2)
```

Taking norms in (1) gives the useful height identity

```text
log N(g_ij)=(d_0i+d_0j-d_ij)/2.                   (3)
```

Thus the pair-distance bounds in the central-truncation extraction give

```text
(1-zeta)W/4 <= log N(g_ij) <= (1+2zeta)W/4.        (4)
```

These are bounds on the actual Gaussian gcds, not merely on gcds of
the core products.

## 2. All pair residue bounds follow at the same exponential scale

Let `r_i=|h_i|` and `r_ij=|h_ij|`. Equation (1) implies

```text
|t_ij| <= r_ij (|t_i|/r_i+|t_j|/r_j).             (5)
```

The actual anchored endpoint estimates are

```text
|t_i|/r_i <= (C/2) exp(-W/4).
```

Combining them with `r_ij=exp(d_ij/2)` and the retained pair upper
bound yields

```text
|t_ij| <= C exp(zeta W/4).                        (6)
```

The direct original arc-diameter estimate has `C/2` instead of `C`;
reconstructing a diameter from two anchor distances loses only this
factor two. Consequently these all-anchor small-value statements are
not additional asymptotic restrictions beyond the actual one-anchor
estimates and pair-distance control.

Even if one remembers only the coarse inequalities from Section 6 of
the central-truncation note,

```text
r_i>=exp(W/4),
|t_i|<=(C/2)exp(zeta W/4),
```

equation (5) still gives

```text
|t_ij| <= C exp(zeta W/2).                        (7)
```

Both (6) and (7) have logarithms negligible relative to every
individual core block when `k=floor(c log_2 M)`, `c<1/2`, since
`zeta=O(k/sqrt(M))` and the block log norm is comparable to
`W/2^(k-1)`.

One must retain the strong pair moments for this conclusion at the
individual-block scale. If only the full-pattern error `eta_core`
is retained, set `G'_ij=product_(S contains i,j) H'_S`. Then

```text
log N(G'_ij)>=(1-eta_core)W_core/4,
log|h_i|<=s+(1+eta_core)W_core/4.
```

For an arbitrary one-anchor system with `|t_i|<=exp(v)`, these
estimates and `G'_ij|g_ij` yield only

```text
|t_ij| <= 2 exp(v+s+eta_core W_core/2).            (8)
```

With the available pattern error, dividing its logarithm by an
individual block scale gives an error of order `k 4^k/sqrt(M)`.
Thus this weaker bookkeeping does not preserve the full range
`c<1/2`. Equation (8) is a limitation of that bound, not a claimed
counterexample to the stronger extracted statement.

## 3. Pair core products and their small correcting factors

For the actual central truncation define

```text
G'_ij=gcd(A'_i,A'_j)=product_(S contains i,j) H'_S,
A'_{ij}=A'_j bar(A'_i)/N(G'_ij).
```

The pair core numerator is explicitly

```text
A'_{ij}=
 product_(j in S, i notin S) H'_S
 product_(i in S, j notin S) bar(H'_S).
```

Central truncation gives `A'_{ij}|h_ij` for every pair by the same
primewise separation estimate used for anchor zero. Write

```text
h_i=K_i A'_i,
h_ij=K_ij A'_{ij},
g_ij=G'_ij L_ij.
```

All displayed correcting quantities are Gaussian integers, and

```text
log|K_i|, log|K_ij|, log|L_ij| <= s.              (9)
```

The bound for `L_ij` is slightly sharper than the generic estimate
from multiplying `K_i,K_j`. At a prime, after removal of the common
core exponent, the gcd has at most `2r_p` additional exponents in
one Gaussian orientation. If both rows contain that core prime,
this follows from their correction exponents. If at least one
does not, its entire exponent in the gcd is at most `2r_p`.
Conjugate coprimality rules out two orientations simultaneously.
Hence `log|L_ij| <= sum_p r_p log p=s`.

Substituting these factorizations into (1) gives the further exact
compatibility identity

```text
N(L_ij) K_ij=K_j bar(K_i).                        (10)
```

Equivalently, in a proposed one-anchor model, pair-core divisibility
holds precisely when `N(L_ij)` divides `K_j bar(K_i)` in `Z[i]`.
In the actual central truncation this divisibility is guaranteed.

## 4. The additional condition is an ordering of residual allocations

Fix a prime with a nonzero core exponent `e'`, and orient it as
`gamma` so that the anchor lies outside its core cut `S`. Write
`b_i` for the signed exponent difference defining `h_i`, in this
orientation, and set

```text
delta_i=b_i-e' indicator(i in S),
delta_0=0.
```

The one-anchor factorization asserts `delta_i>=0` for `i in S`;
for outside rows the signed correction can have either sign.
For a pair on opposite sides of the cut, the signed exponent in
`h_ij` is

```text
e'+delta_j-delta_i               (i outside, j inside).
```

Consequently the proposed core numerator divides every pair numerator
if and only if

```text
max_(i outside S) delta_i <= min_(j inside S) delta_j.   (11)
```

Pairs on the same side impose no core-divisibility condition at this
prime. This proves both necessity and sufficiency.

For the actual central truncation, in either anchor orientation,

```text
delta_i=u_i-u_0                 for i outside S,
delta_j=2r-u_j-u_0              for j inside S.
```

Since every `u_i<=r`, the number `r-u_0` separates the two groups:

```text
delta_i<=r-u_0<=delta_j.
```

Thus the extra all-anchor divisibility is the exact ordering (11),
not an extra independent complex-valued equation of the circle.

## 5. A small-value family where the ordering fails

The following example satisfies the primitive one-anchor equations,
independent conjugate-coprime blocks, asymptotically equal full-pattern
weights, and coefficient and residue heights negligible relative to
each block. It nevertheless fails exact pair-core divisibility.

Take `n>=0`, and put

```text
T=1360+1690n,         gamma=3+2i,
H_empty=T+6+i,
H_1=T+2-i,
H_2=gamma(T+4-i),
H_12=(T+i)/gamma.
```

Here `T` is even, `T=0 mod 5`, and `T=8 mod 169`. The last
condition gives `v_gamma(T+i)=1` and `v_bar(gamma)(T+i)=0`.
Thus `H_12` is integral, and moving this one factor from the
common block into `H_2` preserves all required coprimalities.

For completeness, before that move consider the four linear factors
`T+6+i`, `T+2-i`, `T+4-i`, `T+i` and their conjugates. Every
individual norm is odd. The nonconstant differences relevant to
possible common divisors have norms supported on `2` and `5`;
the additional difference `6` causes no common divisor at the inert
prime `3` because all imaginary parts are units modulo `3`.
At `T=0 mod 5`, direct substitution excludes every possible shared
prime above `5`. Thus these eight factors are pairwise coprime.
Moving the unique factor `gamma` preserves that property.

Define `A_1=H_1 H_12`, `A_2=H_2 H_12`, and
`K_1=gamma`, `K_2=1`. Then

```text
h_1=K_1 A_1=(T+2-i)(T+i)=T^2+2T+1+2i,
h_2=K_2 A_2=(T+4-i)(T+i)=T^2+4T+1+4i.
```

Both `h_i` are coprime to their conjugates: their real parts are
odd, and their imaginary parts have only the rational prime factor
`2`. Their imaginary residues are the constants two and four.
Every block has log norm `2 log T+O(1)`, and the correcting
coefficients are constant.

But the actual common divisor is `g_12=T+i`, while the proposed
core common divisor is `G'_12=(T+i)/gamma`. Equation (1) gives

```text
h_12=(T+4-i)(T+2+i)=T^2+6T+9+2i,
A'_{12}=gamma h_12.
```

Thus `A'_{12}` does not divide `h_12`. Equivalently,
`L_12=gamma` and `N(L_12)=13` does not divide
`K_2 bar(K_1)=bar(gamma)`. At this prime the residual allocation
of the outside row exceeds that of the inside row, violating (11).

The actual pair residue is still the constant two. The example
therefore separates the additional exact divisibility from pair
residue smallness, which remains consistent with Section 2.

Exact Gaussian-integer computations checked the eight-factor
coprimality, primitivity, transition formula, and failed divisibility
for the first 100 values of `n`. The preceding calculations provide
the proof for the entire progression.

## 6. Resulting cocycles and the remaining issue

Eliminating the real coordinates from (2) gives, for distinct
nonreference indices `i,j,l`,

```text
t_i N(g_jl)t_jl-t_j N(g_il)t_il+t_l N(g_ij)t_ij=0.   (12)
```

Inserting the core factors expresses this as an exact relation with
small coefficients among the products
`product_(S contains i,j) N(H'_S)`. This is the familiar rank-two
determinant compatibility, now with the exact correcting factors
`N(L_ij)` retained. It follows algebraically from the one-anchor
coordinates and is not an additional independent relation.

What central truncation adds to an arbitrary small-value one-anchor
system is the coherent primewise ordering (11), or equivalently the
integrality in (10), simultaneously for every pair and cut. Whether
this is a genuinely stronger asymptotic requirement is settled by the
following repair lemma: at the available quantitative coefficient scale,
it is not.

## 7. Negligible trimming repairs every pair simultaneously

The primary agent proposed the following repair, which has been checked
prime by prime. It applies to an arbitrary proposed one-anchor system,
without assuming that its blocks came from central truncation.

Suppose independent conjugate-coprime Gaussian blocks `H_S` are given,
and put

```text
A_i=product_(S contains i) H_S,
h_i=K_i A_i,
gcd(h_i,bar(h_i))=1,
kappa_i=log|K_i|,       Kappa=sum_i kappa_i.
```

Then there are divisors `H''_S|H_S`, keeping the same cut labels and
orientations, such that all pair core numerators formed from the new
blocks divide the actual primitive numerators `h_ij`. Moreover, if
`L` is the total removed block log norm, then

```text
L=sum_S log N(H_S/H''_S) <= 2 Kappa.             (13)
```

The empty block can be retained unchanged. The new anchor and pair
correcting factors satisfy

```text
log|K''_i| <= kappa_i+L/2,
log|K''_ij| <= kappa_i+kappa_j+L/2 <= 2 Kappa.    (14)
```

### Construction and proof

Fix an oriented prime `gamma^e` in a nonempty block `H_S`. Write

```text
r_i=v_gamma(h_i)-v_bar(gamma)(h_i),      r_0=0,
a=min_(i in S) r_i,
b=max_(j outside S, including 0) r_j.
```

The primitive one-anchor factorization gives `a>=e`: each inside
row contains `gamma^e` and contains no conjugate factor. For an
outside row, the signed exponent comes entirely from `K_j`.
Consequently `b>=0` and

```text
b <= sum_j v_gamma(K_j).
```

Keep the exponent

```text
e''=min(e,max(0,a-b)).
```

If `e''>0`, every inside signed allocation exceeds every outside
allocation by at least `e''`. This is exactly the local criterion
for every new pair core to divide its primitive numerator. If
`e''=0`, there is no pair core condition at this prime. Also

```text
0<=e-e''<=b<=sum_j v_gamma(K_j),
```

because `a>=e`. Summing after multiplication by `log N(gamma)`
proves (13): no oriented prime is charged through two core blocks,
and a core prime and its conjugate do not both occur among the blocks.

For the anchor bound in (14), the new coefficient is the old one
times the factors removed from its incident blocks, so its logarithmic
modulus increases by at most `L/2`.

For the pair bound, at this prime write

```text
r_i=e indicator(i in S)+delta_i.
```

Primitivity gives
`|delta_i|=v_gamma(K_i)+v_bar(gamma)(K_i)`. Once the new pair core
has been divided out, its remaining norm exponent is

```text
|r_j-r_i|-e'' |indicator(j in S)-indicator(i in S)|
 <= |delta_i|+|delta_j|+(e-e'').
```

At primes outside all core supports, the full pair exponent is at
most the sum of the two coefficient exponents. Summing over all
primes proves the second bound in (14). The same bound is valid
when a prime disappears completely from the core. Thus the proof
does not assume disjoint support between core blocks and correcting
coefficients, nor bounded prime exponents.

### Consequence at the extracted growing-row scale

For central truncation, `kappa_i<=s<=kW/sqrt(M)`, so
`Kappa=O(k^2 W/sqrt(M))`. An individual block has log norm
comparable to `W/2^(k-1)`. Hence

```text
L/min_S log N(H_S) = O(k^2 2^k/sqrt(M))=o(1)
```

throughout `k=floor(c log_2 M)`, `c<1/2`. Every block therefore
retains its full relative weight, and every new correcting factor
still has height negligible relative to each block. The original
primitive numerators and all their imaginary residues remain exactly
the same.

Accordingly, at this scale a full-pattern primitive one-anchor model
can always be converted to a model with all-pair core divisibility
at negligible loss. Together with Section 2, this shows that the
all-anchor equations do not strengthen the available leading-height
reduction, provided the original strong pair-distance bounds are kept.
The exact counterexample in Section 5 remains valid; its deliberately
moved prime is precisely a negligible factor that this construction
removes from the misassigned core block.
