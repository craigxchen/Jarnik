# Higher-row aggregate parity and the six-row uniform-profile gap

The eight-row aggregate-parity argument has a useful general form.  It
also exposes a sharp limitation: its six-row version would exclude a
near-uniform extracted profile, but there is no fixed-size parity-selection
lemma when the set of split primes is unbounded.

## General balanced-character formula

Take an even number `k>=2` of distinct Gaussian points in one common unit
class, on an arc of radius `R>1` and length at most `C sqrt(R)`, and write
`W=log R^2`. Remove their complete
common Gaussian content, whose weight is `D`.  For each split-prime
threshold layer let `W_r` be the total layer weight whose minority size is
`r`, `1<=r<=k/2`; complementary orientations have the same coefficient.
Assume that at every split prime `p` the total allocation exponent across
the selected `k` rows is even.  For every balanced sign vector `lambda`
(`k/2` plus signs and `k/2` minus signs), the signed exponent is then even.
With common units, [affine allocation rigidity](linear_allocation_affine_rigidity.md) gives the nonvanishing needed
for the square-root character identity; this requires `C<=2`.  The balanced
signed argument sum is at most `k Delta/2`, so the square-root character is
within `k Delta/8` of its target axis.  With `Delta<=C/sqrt(R)`, the height
bound is

```
V_lambda >= W/2 - 2 log(k C/8).
```

There are

```
L_k = (1/2) binom(k,k/2)
```

characters modulo global sign.  For a layer with minority size `r`, put

```
c_(k,r) = (1/2) E |r-2A|,
A ~ Hypergeometric(k, r, k/2).
```

This is the average of `|lambda dot S|/4` over the `L_k` characters.  If
`K_k` denotes the prime-by-prime triangle-inequality cancellation loss,

```
K_k = sum_p log(p) * [
  sum_lambda sum_thresholds |lambda dot S_(p,t)|/4
  - (1/2) sum_lambda |sum_i lambda_i a_i(p)|
],
```

then `K_k>=0` and exact summation gives

```
sum_(r=1..k/2) c_(k,r) W_r
  >= (D + sum_r W_r)/2 + K_k/L_k - 2 log(k C/8).     (1)
```

Equivalently,

```
sum_r (c_(k,r)-1/2) W_r
  >= D/2 + K_k/L_k - 2 log(k C/8).                   (2)
```

This is valid for arbitrary nested prime-power layers; disjoint prime
support is unnecessary.  With arbitrary row units, use `C<=sqrt(2)` for
the affine-rigidity input, assume even total row-unit parity, and replace
`k C/8` by `k sqrt(2) C/8` in the height bound.  The displayed formulas
below use the common-unit normalization; the arbitrary-unit version has
that replacement in each logarithmic constant.

The first coefficients are

```
k       c_(k,1)  c_(k,2)  c_(k,3)  c_(k,4)
4          1/2      1/3       -        -
6          1/2      2/5      3/5       -
8          1/2      3/7      9/14    18/35
10          1/2      4/9      2/3      4/7
```

The checker evaluates these coefficients directly from all balanced sign
vectors and checks (1)'s coefficient identity.

## Six rows would contradict a uniform extracted profile

For `k=6`, the nonconstant intrinsic cuts are the `31` unoriented cuts,
with multiplicities `6,15,10` at minority sizes `1,2,3`.  If every cut has
weight `w=W'/31`, then the left side of (2), with `D=0`, is

```
6(0)w + 15(-1/10)w + 10(1/10)w = -w/2 = -W'/62.    (3)
```

The right side is at least `-2 log(3C/4)`, since `K_6>=0`.  Hence a six-row
tuple whose intrinsic cut weights are sufficiently close to uniform and
whose radius is sufficiently large cannot satisfy primewise aggregate
parity.  For example, if every cut lies within relative error `1/20` of
`w`, then the left side is at most

```
(1/10)[10(1+1/20)-15(1-1/20)]w = -3w/8,
```

Consequently (2) implies the stronger bound on the total norm weight

```
D/2 + 3W'/248 <= 2 log(3C/4),
W=D+W' <= (496/3) log(3C/4).                       (4)
```

Here the second line uses `D>=0` and `1/2>=3/248`. In particular no such
tuple exists once `W>(496/3) log^+(3C/4)` in the stated range `0<C<=2`.
When `3C/4<=1`, (4) already contradicts `R>1`.

This is a genuine higher-row arithmetic limitation on a near-uniform
profile, and it keeps exact heights throughout.  It does not yet apply to
an arbitrary large endpoint cluster because the required six-row parity
subtuple need not exist.

## Why fixed-size parity selection fails with unbounded support

Reduce every row's allocation vector modulo two, obtaining vectors in an
F2-space whose coordinates are the split primes.  Primewise aggregate
parity for a six-row subset means that its six vectors sum to zero.  There
is no such guarantee for an arbitrarily long sequence when the ambient
dimension is unbounded.  For `M` rows and distinct split primes
`p_1,...,p_M`, take row `i` to have allocation parity `1` only at `p_i`
(and zero elsewhere).  Every six-row subset then has six distinct nonzero
coordinates in its sum, so none has zero parity.

This parity pattern is compatible with equal Gaussian norms at the formal
allocation level: choose exponent `e_{p_i}=1` at each prime, put allocation
one in row `i` and zero in the other rows, and use the conjugate exponent
for the complementary allocation.  To make the associated cut profile
arbitrarily close to uniform while retaining this parity obstruction, add
for every nonconstant unoriented cut `S` a dominant even column with
allocation `2Q` on `S` and zero off `S`.  Its `2Q` threshold copies have
even exponent and contribute weight `2Q` to that cut.  The private odd
columns add only weight one to the singleton cuts, so every cut has relative
error at most `1/(2Q)` from the uniform profile. The total formal weight is
`2Q(2^(M-1)-1)+M`, and every pair distance is `2Q*2^(M-2)+2`.
Their difference from half the total is `Q+2-M/2`, which is strictly
positive when `Q>(M-4)/2`. The private columns give full affine rank. This
remains a formal equal-norm allocation fixture, not a claim about a short
endpoint arc.  It proves that affine independence, equal norms, an
arbitrarily near-uniform fixed-cut profile, and unbounded split-prime
support alone cannot supply the missing fixed six-row parity subtuple.

The unit-weight description can also be implemented with actual distinct
split primes, without equating their logarithms. Fix one prime `p_S` for
each full cut and separate private primes for the rows. At cut `S` use
the even exponent `e_S=2 floor(H/(2 log p_S))`. Its log-norm weight is
`H+O(log p_S)`; the private exponents remain one. For fixed `M`, as
`H` tends to infinity every full-cut weight is `H+O_M(1)`, the parity
vectors are unchanged, and the pair-distance excess is `H/2+O_M(1)>0`.
Any fixed selected row set has asymptotically equal intrinsic cut weights
as well, by counting the full cuts restricting to each selected cut.
Taking products of fixed Gaussian primes with these exponents gives
literal equal-norm Gaussian points and full affine rank. Their angles
are uncontrolled; no endpoint realization is asserted.

More generally, a zero-sum theorem in the F2 allocation space only gives a
zero-sum subset whose size depends on the dimension.  Passing to such a
subset does not preserve a fixed `k`; already at `k=8` the uniform
coefficient margin has changed sign (the six-row margin is negative but the
eight-row margin is positive), so this six-row contradiction cannot simply
be extended by taking a larger zero-sum subset.

If the selected six rows have odd total parity at even one split prime,
then every balanced character has an odd signed exponent at that prime:
the signs `lambda_i` are all odd modulo two, so the parity is independent
of `lambda`. The quarter-character exponents needed for `Beta` are then
half-integral. The Gaussian half-character `A_lambda` has
`Norm(A_lambda)=t h_lambda^2`, and its normalization
`A_lambda/(h_lambda sqrt(t))` lies in the common field `Q(i,sqrt(t))`.
As [the quadratic-twist note](odd_layer_character_quadratic_twist_obstruction.md)
explains, this does not supply a Gaussian square root `Beta` or even place
all quarter roots in that same quadratic field. The axis/diagonal
perpendicular-coordinate gap for `Beta` is unavailable. This is an exact
failure of the parity step, rather than a loss that can be absorbed into
`O_C(1)`.

The conclusion is therefore precise: higher-row aggregate parity yields the
conditional six-row obstruction (3), but no parity extraction of the
required fixed size follows without a bound on split-prime support or a new
geometric restriction on parity vectors.

[The exact coefficient and parity audit](check_higher_row_aggregate_parity_uniform_limit.py)
checks the coefficient tables, the parity independence, and the even
full-cut/private-odd formal family.  It contains no endpoint search and
makes no claim about realizability of that family on a short arc.
