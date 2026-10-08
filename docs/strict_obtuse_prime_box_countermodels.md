# Prime-weight integer boxes at linear support and logarithmic radius cost

The weighted obtuse-box conditions do not force a superlinear number of
varying prime coordinates or a radius cost larger than `O(M log M)`.
This remains true with genuine distinct split primes, all allocation
widths equal to one, the stronger common-unit pair-norm separation, and
enough separation to accommodate all pairwise lower chord bounds on an
abstract short arc. A second construction also makes every unordered
pair norm distinct.

These are countermodels to deductions from that package of necessary
conditions. They are literal primitive Gaussian circle configurations,
but their actual arguments are uncontrolled. They are not asserted to
lie on short arcs, and do not contradict the fixed-rank phase/Roth
obstruction or prove anything negative about the original uniformity
conjecture.

A later [exact phase audit](linear_support_character_height_obstruction.md)
excludes the canonical five-copy flip assignment used by the checker
from fixed-constant endpoint arcs once `M` is large. Its normalized
containing-arc constant is at least a constant times `M^3` in the prime
family below. This supplies information beyond the necessary norm
conditions tested here; their countermodel conclusions remain valid.

The [exact checker](check_strict_obtuse_prime_box_countermodels.py) tests
both constructions on actual split primes for `M=4,8,16`.

## 1. Many nearby prime weights from a dyadic pigeonhole argument

Fix `C>0`, let `M=2^t -> infinity`, and put

```text
X=(8M/C)^4.
```

For either `r=M-1` or `r=5(M-1)`, there are, for all sufficiently large
`M`, distinct primes `p_1,...,p_r`, all congruent to one modulo four,
in one interval

```text
[A,A+X/r),       X<=A<=2X.
```

Indeed, the fixed-modulus prime number theorem gives
`pi(2X;4,1)-pi(X;4,1) ~ X/(2 log X)`. Partition `[X,2X]` into `r`
equal intervals. Their total prime count exceeds `r(r-1)` eventually,
so one contains at least `r` primes. This uses no short-interval prime
number theorem. The needed fixed-modulus result is the classical PNT
for arithmetic progressions; see
[Vaughan, Chapter 11](https://personal.science.psu.edu/rcv4/568s20/568chapter11.pdf).

Set `ell=log A` and `log p_j=ell+epsilon_j`. Then

```text
0<=epsilon_j<1/r,
W=sum_j log p_j=r ell+O(1)=4r log M+O_C(r).            (1)
```

Thus these are distinct integer prime norms with almost equal logarithmic
weights, not formally assigned equal real weights.

## 2. Minimal support and an almost scalar positive spectrum

Index rows by `x in F_2^t` and columns by nonzero `a in F_2^t`. Set

```text
b_x(a)=x.a mod 2,       e_a=1.
```

Every column varies, and every pair of distinct rows differs in exactly
`M/2` columns. Hence, using the primes from Section 1 with `r=M-1`,

```text
d_xy-W/2 > ell/2-1/2 >= (log X)/2-1/2.                (2)
```

For `f_x(a)=sqrt(log p_a)(2b_x(a)-1)`, binary endpoints make the
linear embedding inequality an equality:

```text
<f_x,f_y>=W-2d_xy < -log X+1 < 0.
```

The truncated Walsh matrix has orthogonal columns of squared norm `M`.
Consequently the nonzero eigenvalues of the row Gram matrix are exactly

```text
M log p_a,        a!=0.                              (3)
```

There is one zero eigenvalue, with positive kernel vector consisting of
all ones. All other eigenvalues lie in an interval of length less than
`M/r`, around `M ell`. In particular this is an exceptionally stable
simplex, with affine rank `M-1`, smallest possible widths, and smallest
possible support `r=M-1`. Its radius satisfies

```text
log R=W/2=(2+o(1)) M log M.                           (4)
```

This rules out improving dimension, width, or spectral stability merely
from strict obtuseness, prime-weight integrality, or the lower chord
packing constraints checked in Section 4. Pair norms repeat in this
minimal-support model. The next construction removes that repetition.

## 3. Distinct pair norms with only five times as many coordinates

Take five copies of every nonzero Walsh column, so `r=5(M-1)`.
Choose an injection `j:F_2^t -> {1,...,r}` and flip exactly coordinate
`j(x)` in row `x`. Denote the modified binary rows by `c_x`.

Each coordinate has originally `M/2` zeros and `M/2` ones, and is flipped
in at most one row. Thus every coordinate still varies for `M>=4`.
For distinct rows, their two flips change Hamming distance by at most
two. Therefore

```text
dist_H(c_x,c_y)>=5M/2-2 = r/2+1/2.                   (5)
```

With `r` nearby distinct primes as above, (2) and strict negativity
hold unchanged. There are only `5(M-1)` varying primes, all widths are
one, and

```text
W=(20+o(1)) M log M,
log R=(10+o(1)) M log M.                              (6)
```

Every unordered pair has a different difference support. To prove this,
suppose the supports for `{x,y}` and `{u,v}` agree. Before the flips,
the difference words are the fivefold Walsh words of `x+y` and `u+v`.
If these two indices differ, the words differ in `5M/2>4` coordinates,
which four flips cannot remove. Thus `x+y=u+v`. Cancellation of the
Walsh words now gives

```text
e_(j(x))+e_(j(y))=e_(j(u))+e_(j(v))    over F_2.
```

Injectivity of `j` implies equality of the unordered pairs. Since the
prime labels are distinct, the integers

```text
P_xy=product_(columns where c_x!=c_y) p_j
```

are therefore distinct for all `binom(M,2)` pairs. In particular the
pair-norm injectivity required by actual sufficiently short arcs does
not repair the superlinear-support or super-`M log M` inference.

## 4. The lower chord bounds can all fit in an abstract short arc

For both constructions, (2) implies

```text
P_xy^2/N > X/e = (8M/C)^4/e
              > 256(M-1)^4/C^4,       N=exp(W).        (7)
```

The common-unit integer chord lower bound is `2R/sqrt(P_xy)`.
Since `R=sqrt(N)`, (7) gives

```text
2R/sqrt(P_xy) < C sqrt(R)/(2(M-1)).                   (8)
```

Put `M` abstract equally spaced directions on an arc of angular length
`C/sqrt(R)`. For sufficiently large `M`, all half-angle gaps are at most
one. The elementary inequality `sin u>=u/2` on `[0,1]` shows that the
smallest chord of these abstract points is at least
`C sqrt(R)/(2(M-1))`. Thus every pair's integer chord lower bound fits
strictly below the corresponding chord, including every consecutive
pair. This is stronger than merely keeping the fixed additive slack
`log(4/C^2)` in the pair-norm inequality.

The abstract directions are not claimed to be the directions determined
by the selected Gaussian primes. Precisely that phase identification is
missing; satisfying lower chord bounds is not an existence proof for a
short arc.

## 5. Literal Gaussian realization and scope

Choose a Gaussian prime `pi_j` of norm `p_j`, with positive coordinates,
and use the same literal unit in every row:

```text
z_x=product_j pi_j^c_x(j) bar(pi_j)^(1-c_x(j)).
```

Every `z_x` is a Gaussian integer of norm `N`. Since each allocation
column takes both endpoint values, their complete common Gaussian gcd
is a unit. Unique factorization gives exactly the stated pair norm
`P_xy`. The same construction using `b_x` realizes the minimal-support
model. Thus prime integrality, Gaussian integrality, common units,
primitive normalization and all norm calculations are genuine.

The minimal-support model has repeated pair norms and is already
incompatible with actual short-arc norm injectivity. The modified model
removes that particular defect. Neither model verifies small imaginary
chord numerators or simultaneous Gaussian phase concentration. Those
conditions can provide information that weighted boxes, their spectra,
distinct pair norms and lower chord packing inequalities do not.
