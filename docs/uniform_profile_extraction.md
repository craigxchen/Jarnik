# A backwards reduction for arbitrarily slowly growing clusters

## Status and purpose

This is a prose proof of a fixed-size profile extraction theorem. Its
Hilbert-space Bessel ingredient is also checked in
`GaussianChain/ObtuseBessel.lean`; the complete extraction is not formalized. It
applies to every sufficiently large endpoint cluster, with no assumption
that its cardinality approaches a logarithmic or sublogarithmic bound.
It does not prove the requested uniform bound. Its value is to replace
the general mixed-profile obstruction by one specified limiting profile.

## 1. Endpoint geometry gives an obtuse family of sign vectors

Work with distinct points in one common Gaussian-unit class, remove their
common Gaussian divisor, and write the remaining squared radius as
`N=R^2>1`. Removing a common divisor can only decrease the endpoint arc
constant. All inert and ramified common factors can be removed this way.
Use the split-prime representation

```
z_i = epsilon product_p pi_p^(a_i(p)) bar(pi_p)^(e_p-a_i(p)),
W = log N = sum_p e_p log p.
```

On the finite set of prime threshold layers `(p,t)`, `1<=t<=e_p`,
put probability weight `log p/W` and the sign function

```
f_i(p,t) = 2 indicator(a_i(p)>=t)-1.
```

These are unit vectors in the resulting real weighted `L^2` space.
Their inner products are

```
<f_i,f_j> = 1-2 d_ij/W,
d_ij = sum_p |a_i(p)-a_j(p)| log p.
```

The exact primitive-chord identity from `balanced_bonus_phase_audit.md`
is `|z_i-z_j|^2=4 t_ij^2 exp(W-d_ij)`, where `t_ij` is a positive
integer. If the arc constant is `C<=2`, it follows that

```
d_ij >= W/2+2 log(2t_ij/C) >= W/2,
<f_i,f_j> <= 0                  (i != j).            (1)
```

Only the weak final inequality will be used. The construction works
with arbitrary prime exponents and has no lower bound on `M/log R`.

## 2. An elementary Bessel inequality for obtuse unit vectors

If unit vectors `f_1,...,f_M` in a real Hilbert space satisfy (1), then
for every vector `h`,

```
sum_i <f_i,h>^2 <= 2 ||h||^2.                       (2)
```

Here is a proof avoiding a spectral estimate. Put `a_i=<f_i,h>`.
For the indices with `a_i>0`, let `v=sum_i a_i f_i` and
`S=sum_i a_i^2`. Nonpositive cross inner products give `||v||^2<=S`,
whereas `<v,h>=S`. Thus

```
0 <= ||h-v||^2 = ||h||^2-2S+||v||^2 <= ||h||^2-S.
```

So `S<=||h||^2`. Apply the same argument to `-h` and the remaining
strictly negative coefficients. Adding gives (2).

## 3. Extracting a fixed number of almost independent rows

Fix `1<=k<=M`. Select an ordered tuple of `k` distinct rows uniformly.
For a nonempty subset `S` of its positions, write

```
mu_S = E_layers product_(j in S) f_(i_j).
```

If `|S|=r`, condition on `r-1` of these row indices. Their product is
a sign function `h`, hence `||h||=1`. The last row is uniform among
`M-r+1` choices. Equation (2) therefore gives

```
E_rows mu_S^2 <= 2/(M-r+1).
```

Summing over subsets proves the existence of a tuple satisfying

```
sum_(empty != S subset [k]) mu_S^2
 <= sum_(r=1..k) binom(k,r) 2/(M-r+1)
 <= 2(2^k-1)/(M-k+1).                               (3)
```

Let `P(epsilon)` be the fraction of the inherited conductor weight
whose signs on this tuple equal `epsilon in {+1,-1}^k`. Fourier
inversion on the binary cube gives

```
P(epsilon) = 2^(-k) sum_(S subset [k]) epsilon_S mu_S,
mu_empty = 1.
```

By (3) and Cauchy--Schwarz, every pattern simultaneously obeys

```
|2^k P(epsilon)-1|
 <= (2^k-1) sqrt(2/(M-k+1)).                        (4)
```

In particular, given fixed `k` and `eta>0`, it suffices that

```
M >= k-1+2(2^k-1)^2/eta^2                          (5)
```

to extract a tuple for which every inherited pattern weight is within
relative error `eta` of `2^(-k)`. These constants are deliberately
unoptimized; their independence of the radius is the essential point.

## 4. Intrinsic normalization and the eight-point limit

Now divide out the selected tuple's own common Gaussian divisor.
Exactly the two constant sign patterns are removed from its conductor.
All other layer patterns are preserved. If their combined removed
weight is `aW`, the new total weight is `W'=(1-a)W`.

As `M` tends to infinity for fixed `k`, (4) gives `a -> 2^(1-k)`.
Every nonconstant oriented pattern therefore has intrinsic weight

```
(1/(2^k-2)+o(1)) W'.
```

Equivalently, after identifying a cut with its complement, every one
of the `2^(k-1)-1` nontrivial cuts has weight
`(1/(2^(k-1)-1)+o(1)) W'`. This identification concerns cut weights;
it does not equate the Gaussian factors at the two orientations.

Quantitatively, for `k=8`, if (4) is at most `eta<=1`, each intrinsic
unoriented cut has relative error at most `128 eta/(127-eta)` from
`W'/127`. This follows by dividing its inherited weight
`(1+error)/128` by `1-a`, with `|128a-1|<=eta`.

For `k=8`, the limiting intrinsic profile is uniform on 127 cuts.
In the notation of the intrinsic eight-point target,

```
D/W'   -> 240/127,
W_4/W' -> 35/127,
(D+W_4-2W')/W' -> 21/127 > 0.                     (6)
```

The arithmetic is

```
sum_(r=1..7) binom(8,r)(r-4)^2 = 480,
sum_(r=1..7) binom(8,r) = 254,
binom(8,4) = 70.
```

For example, in the neighborhood where every intrinsic cut weight
has relative error at most `1/20` from `W'/127`, one has the finite
inequality

```
D+W_4-2W' >= (217/2540) W'.                         (7)
```

To check it, the coefficients of `D+W_4-2W'` on cut sizes
`1,2,3,4` (up to complement) are `7,2,-1,-1`, with respective
multiplicities `8,28,56,35`. Their sum is 21 and their absolute
sum is 203. Thus the lower bound is
`(21-203/20)W'/127=217W'/2540`.

Thus proving the existing sufficient inequality
`D+W_4<=2W'+B` merely in a fixed neighborhood of this one profile
would already suffice, provided `B` is independent of the radius.
There is no need to establish it for every mixed eight-point profile.

The normalized radii of the extracted tuples must tend to infinity
in a hypothetical sequence with `M` tending to infinity:
the original `M` sign vectors satisfy (1) with the strict gap
`<f_i,f_j><=-4 log(2/C)/W` when `C<2`, and hence `W` grows with `M`.
Since `W'/W -> 1-2^(1-k)>0` by (4), `W'` also tends to infinity.
Indeed `0<=||sum_i f_i||^2` gives
`W>=4(M-1)log(2/C)`. We may fix `C=1/2`, so the strict-gap argument
always applies to the uniformity reduction.

## 5. Exact missing statement

A sufficient next theorem is:

> There exist fixed `eta_0>0` and a height threshold such that no
> primitive, common-unit eight-point arc of length at most
> `(1/2)sqrt(R)` above that threshold has all 127 intrinsic cut
> weights within relative error `eta_0` of `W/127`.

This statement is **unproved**. The extraction proves that excluding
this fixed neighborhood would rule out any unbounded sequence, including
sequences growing more slowly than every proposed unbounded rate.
Eight is chosen to match the existing local arithmetic target; excluding
this particular eight-point neighborhood is sufficient, not necessary.
If it can occur, the extraction remains available for any larger fixed `k`.
Bounded radii contribute only finitely many points and cause no problem.
One common unit costs a factor of at most four in selecting the original
large cluster; arbitrary fixed `C` is handled by subdividing into arcs
with constant at most `1/2`.

There is also an explicit conditional cardinality bound. Assume the
boxed assertion holds whenever `W>=W_0`, and reduce `eta_0` if needed
so `0<eta_0<=1`. A common-unit cluster of size at least

```
max(8,
    7+520200/eta_0^2,
    1+16 max(W_0,0)/(63 log 4))
```

would be contradictory: use inherited error `eta=eta_0/2` in (5),
then the intrinsic error is at most `eta_0`. Also `W'/W>=63/64`
and the strict-gap estimate give `W'>=(63/16)(M-1)log 4>=W_0`.
Integer ceilings, a factor four for units, and arc subdivision give
an unconditional-in-`R` count *conditional on the missing assertion*.

The nongenuine Gaussian relaxations from `integer_combination_separation.md`
and `finite_local_test_attack.md` still explain why the individual
divisor, determinant, and residue-size inequalities do not finish this
step. The missing assertion must use the actual Gaussian arguments and
integrality. Profile convergence by itself is not a contradiction.
In particular, a relative error `eta` in a logarithmic weight allows
an error of order `eta W` in that logarithm. Even when `eta` tends
to zero along a sequence, the absolute error need not be bounded.
No continuity argument passing from exact equal weights to a fixed
neighborhood is justified without controlling this issue.

## 6. Larger fixed tuples allow an arbitrarily small positive bonus

There is no need to insist on the size or coefficient of the eight-point
target. For even `k=2s`, define

```
D = sum_layers (r-s)^2 log p,
W_s = sum_(layers with r=s) log p.
```

The same exact product of primitive chord identities gives

```
D+2 log T_res+2 S_arc
 = (k/4)W+k(k-1)log(C/2),
T_res = product_(i<j) t_ij.
```

At the uniform intrinsic profile on the `2^k-2` nonconstant oriented
patterns,

```
D/W = k/4 - k(k-1)/(2(2^k-2)),
W_s/W = binom(k,k/2)/(2^k-2).
```

Consequently, a proof near that profile of

```
D+delta W_s <= (k/4)W+B                            (8)
```

with radius-independent `B` would give uniformity whenever

```
delta binom(k,k/2) > k(k-1)/2.                     (9)
```

Every fixed `delta>0` satisfies (9) for some sufficiently large even
`k`. For example, `k=32` and `delta=1/1000000` work, since
`binom(32,16)=601080390` and `32*31/2=496`.
The positive limiting excess persists in a sufficiently small fixed
profile neighborhood; Sections 3--5 then give a finite bound independent
of `R`. Its size is immaterial to the requested uniformity.

Thus a much weaker local arithmetic gain than the eight-point
coefficient one could suffice if it is available on a larger fixed
tuple. No such new positive gain has been proved here. In particular,
one cannot assume a family of gains `delta_k` works: it must satisfy
the explicit comparison (9), rather than decay too rapidly with `k`.

## 7. A checked limitation of passing to a function-field model

A possible next idea was to interpret a divergent sequence as a family
over a function field and apply a uniform theorem about squares of
quadratic polynomials. A relevant primary result is Pasten,
[Representation of squares by monic second degree polynomials](https://arxiv.org/pdf/1003.1969),
Theorem 2.7: over a characteristic-zero function field of a curve of
genus `g`, a monic quadratic taking square values at at least
`max(8,4(g+1))` distinct constants is constant-coefficient or a square.
His number-field Theorem 2.3 assumes Bombieri's conjecture and permits
a finite exceptional set depending on the chosen evaluation points.

Neither theorem is a uniform result for the integer configurations here.
The extraction above gives no curve of controlled genus, no square-value
parametrization at distinct constants, and no control of a number-field
exceptional set. Passing to a subsequence or observing convergence of
normalized weights does not supply these missing hypotheses. This route
has therefore not been used to assert a uniform bound.

## Verification boundary

The project build passes with the new `ObtuseBessel` module imported.
Its three lemmas use only `propext`, `Classical.choice`, and `Quot.sound`.
Exact rational checks confirm the eight-point coefficients and the
averaged moment bounds on finite orthogonal character families.
These checks support the stated intermediate results; they do not verify
the unproved exclusion in Section 5 or the full endpoint theorem.
