# Angular moment cancellation and a square-root product-length criterion

This strengthens the earlier
[polynomial product-rigidity bounds](short_arc_product_rigidity.md).
It applies to arbitrary Gaussian lattice points, but does not prove a
uniform bound on their number.

The later [linear-allocation argument](linear_allocation_affine_rigidity.md)
gives rigidity for products of every length whenever `C<=sqrt(2)`,
with arbitrary Gaussian units. Its real affine-independence theorem
also gives `M<=r+1` for r varying split primes. The present moment
theorem remains valid, including its larger-C cases.

## Theorem

Let `k` be a positive integer and `C>0`, with

```text
k C^2 <= 8.
```

On every centered-circle arc of length at most `C sqrt(R)`, equality
of products of `k` nonzero Gaussian lattice points implies equality
of the two multisets, including multiplicities. The same holds for
every smaller product length. In particular, every `C=1/2` arc is a
multiplicative `B_32` set, without a radius threshold.

The assertion for `k=1` is immediate. Empty or singleton arcs are
also immediate. In the proof below `k>=2` and `R>=1`.

## 1. The angular lift or a singleton arc

Let `Delta` be the angular width of the arc. If `k Delta < 2 pi`,
choose all arguments relative to the midpoint in
`[a,b]=[-Delta/2,Delta/2]`. Equality of the two products makes the
two argument sums equal modulo `2 pi`. Their difference has absolute
value at most `k Delta`, so the sums agree in the real numbers.

The other case is harmless. If `k Delta >= 2 pi`, the endpoint bound
`Delta <= C/sqrt(R)` gives

```text
sqrt(R) <= k C/(2 pi),
arc length <= C sqrt(R) <= k C^2/(2 pi) <= 4/pi < sqrt(2).
```

Distinct equal-norm Gaussian integers have squared distance a positive
even integer, so their distance is at least `sqrt(2)`. The arc then
contains at most one point. Thus every nontrivial case has the required
real equality of argument sums. This also disposes of any possible
wraparound of the arguments.

## 2. Equal argument sums force equal Gaussian sums

We first prove a sharper replacement for the Taylor bound in the
earlier note. Suppose two multisets of `k` arguments in `[a,b]` have
the same mean `mu`. Let `A,B` be their uniform probability measures,
and let `E` put masses at the endpoints `a,b` with this same mean.
For `s` in `[a,b]`, define

```text
K_A(s) = integral (x-s)_+ dA(x) - (mu-s)_+,
```

and define `K_B,K_E` likewise. Convexity of `x -> (x-s)_+` gives

```text
0 <= K_A(s), K_B(s) <= K_E(s).
```

The lower bound is Jensen's inequality. For the upper bound, the
linear interpolation of this convex function between `a,b` lies above
it, and its expectation depends only on the mean.

Twice integrating the second derivative of any complex-valued `C^2`
function `f` gives

```text
integral f dA - integral f dB
  = integral_a^b f''(s) (K_A(s)-K_B(s)) ds.
```

This follows as well from the exact formula
`f(x)=f(a)+f'(a)(x-a)+integral_a^b f''(s)(x-s)_+ ds`.
For `f(s)=exp(i s)`, it follows that

```text
|integral exp(i x) dA - integral exp(i x) dB|
 <= integral_a^b K_E(s) ds
 = (mu-a)(b-mu)/2
 <= (b-a)^2/8.
```

The middle equality is half the variance of the endpoint measure.
For `b>a`, the final bound on the complex difference is in fact
**strict**. If `K_A-K_B` vanishes identically, the difference is zero.
Otherwise this real continuous function is nonzero on an interval.
On such an interval the argument of `f''(s)=-exp(i s)` varies, so
equality in the complex integral triangle inequality is impossible.

Returning to the circle, multiplication by its midpoint phase and
radius therefore gives

```text
|sum z_j - sum w_j| < k R Delta^2/8 <= k C^2/8 <= 1.
```

The difference is a Gaussian integer and must be zero. The case
`Delta=0` already has identical multisets. Consequently the first
elementary symmetric coefficients `e_1(z),e_1(w)` agree, including
the boundary `k C^2=8`.

## 3. An exact cancellation function for each higher moment

For a positive integer `r`, define the complex-valued function

```text
E_r(theta) = integral_0^theta
    i r exp(-i(r-1)s) (exp(i s)-1)^(2r-1) ds.
```

Expanding the binomial inside the integral shows that

```text
E_r(theta) = exp(i r theta) + A_r theta + B_r
              + sum_(0<|j|<r) c_(r,j) exp(i j theta).
```

Indeed, the derivative has frequencies from `-(r-1)` through `r`.
The coefficient at frequency `r` is `i r`, which integrates to
`exp(i r theta)`. Frequency zero integrates to a linear term; all
other frequencies are strictly between `-r` and `r`.

The integral definition also gives the useful bound, independent of
the size of the interpolation coefficients,

```text
|E_r(theta)|
 <= r integral_0^|theta| |s|^(2r-1) ds
 = |theta|^(2r)/2.                                  (1)
```

Here `|exp(i s)-1| <= |s|` was used. The bound holds for both signs
of `theta`.

## 4. Induction through the integer coefficients

Suppose `e_j(z)=e_j(w)` for every `j<r`. Newton's identities imply
equality of the power sums through degree `r-1`. After dividing each
point by the circle midpoint `q`, these are the equalities of the
positive frequency sums `sum exp(i j theta)`. Taking complex
conjugates gives the corresponding negative frequency equalities.
The constant sums agree, and Section 1 supplied equality of the
argument sums.

Every term in the expansion of `E_r`, except its frequency `r`, thus
cancels between the two multisets. Bound (1) yields

```text
|sum z_j^r - sum w_j^r|
 <= k R^r (Delta/2)^(2r)
 <= k (C^2/4)^r.
```

Newton's identity at degree `r`, with all earlier coefficients equal,
now gives

```text
|e_r(z)-e_r(w)| <= (k/r) (C^2/4)^r.                 (2)
```

For `k>=4` and `2 <= r <= floor(k/2)`, use `C^2/4<=2/k`.
At `r=2` the bound is at most `2/k<=1/2`. Successive bounds have
ratio at most `(2/k) r/(r+1)<1`. Hence all these coefficient
differences have modulus strictly less than one and vanish in `Z[i]`.
For `k=2,3` no higher induction step is needed.

Finally, if the common product is `p`, equal moduli imply

```text
e_(k-r)(z) = p conjugate(e_r(z))/R^(2r),
e_(k-r)(w) = p conjugate(e_r(w))/R^(2r).
```

The reflected coefficients therefore agree as well. Both monic root
polynomials have all coefficients equal, so their root multisets are
identical. Padding shorter products by the same arc point and then
canceling proves the statement for all lengths at most `k`.

## Consequence and unresolved step

An `M`-point `C=1/2` cluster has exactly `binom(M+31,32)` distinct
unordered thirty-two-factor products. Thus no nontrivial relation with
at most thirty-two factors on each side can occur, even when indices
repeat. The theorem needs no conductor pattern, limiting template,
coprimality, or common Gaussian-unit class.

The product points, however, have radius `R^32` and angular width at
most `32 Delta`; their containing arc has length at most
`32 C R^(63/2)`, or exponent `63/64` relative to the new radius.
This distinct-product count does not give the requested uniform
cardinality bound on the original arc. The missing step is still an
arithmetic bound on configurations with no such short relations.

The [balanced-pattern rank audit](relation_free_balanced_patterns.md)
shows that the complete half-cut obstruction has no affine exponent
relations of any length. Increasing `k` alone therefore cannot exclude
that pattern or repair its existing short-phase relaxation.

This is a prose proof. The finite-difference identities underlying the
cancellation functions were checked by exact integer arithmetic for
`1<=r<=32`; the displayed binomial expansion proves them for all `r`.
The equal-norm `sqrt(2)` separation used in Section 1 is now checked
in `GaussianChain/DependentEndpoint.lean`. The angular integral,
moment induction, and complete rigidity theorem remain prose proofs.
