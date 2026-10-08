# Squareclass packing baselines and the correct higher-power group

This note separates elementary uniform bounds in restricted classes
from the new six-wise restriction in
[norm_squareclass_separation.md](norm_squareclass_separation.md).
It also specifies the group in which a higher-power relation must be
stated. The six-wise theorem is a real additional restriction, but
bounded squareclass rank already gives uniformity by elementary
arguments; it should not be credited with first proving that special
case. No bound uniform in the unrestricted radius is proved here.

## 1. Units must be included in the squareclass label

Remove the common Gaussian divisor of a cluster. This can only
decrease its normalized arc constant. Its squared radius is then
an odd product of split primes, and its points can be written

```text
z_i=epsilon_i product_p pi_p^a_ip bar(pi_p)^(e_p-a_ip),
epsilon_i/epsilon_0=i^(u_i),
W=log R^2=sum_p e_p log p.
```

Use the augmented parity label

```text
lambda_i=(u_i mod 2, (a_ip-a_0p mod 2)_p).         (1)
```

The anchor has label zero. This is equivalently the squareclass
of the norm of a coordinate-primitive Gaussian numerator `H_i`
representing the **full** ratio `z_i/z_0=H_i/bar(H_i)`.
For an odd unit ratio, `H_i` includes one factor `1+i` or `1-i`,
and its norm acquires the prime-two squareclass. Thus (1) is the
ordinary canonical allocation norm squareclass augmented by one
unit-parity bit.

That bit cannot be omitted. For primitive `gamma=a+i b` solving
`a^2+2ab-b^2=+/-1`, the points

```text
z_0=bar(gamma)^2,       z_1=i gamma^2
```

have common radius `a^2+b^2` and distance `sqrt(2)`. Their
canonical allocation norm labels agree, but their unit parities
differ. The full half-angle numerator is `(1+i)gamma^2`, up to
a rational content and a unit, so its norm has squareclass two.
This Pell family rules out any all-unit claim based only on the
old allocation norm squareclasses. Alternatively, one may partition
the points into the two unit classes modulo `+/-1`.

Let `d` be the binary dimension of the augmented labels. If the old
allocation norm squareclasses have dimension `d_0`, then only
`d<=d_0+1` is automatic.

## 2. The elementary pairwise bound

Suppose an arc has length at most `C sqrt(R)` with `C<2`.
If two labels in (1) agree, their unit ratio is `+/-1`, and their
primitive Gaussian half-angle numerator `g` has square norm.
Furthermore

```text
log N(g)=sum_p |a_ip-a_jp| log p<=W,
|g|<=R.
```

Since the points are distinct, `Im g!=0`. The elementary
all-unit square-norm separation estimate gives

```text
|Im g|/|g| >= |g|^(-1/2) >= R^(-1/2).
```

The arc gives the opposite bound `C/(2 sqrt(R))`, a contradiction.
Thus the labels are distinct and

```text
M<=2^d.                                          (2)
```

This bound requires no near-uniform profile or central truncation.
For the old canonical labels alone the safe all-unit conclusion
is `M<=2^(d_0+1)`, or `M<=2^d_0` inside a fixed unit-parity class.

## 3. A stronger elementary baseline: four-wise separation

For `C<1`, an arbitrary endpoint cluster already has **no nonempty
relation supported on at most four nonanchor augmented labels**.
This again requires no profile extraction. The zero-word case below
is essential: an averaged short multiplicative word may be exactly
one, in which case its imaginary part gives no contradiction.

First consider four distinct points whose four augmented labels sum
to zero. Average the three balanced sign partitions

```text
(+,+,-,-),       (+,-,+,-),       (+,-,-,+).
```

At one prime, scale its allocations to `x_1,...,x_4 in [0,1]`.
The average of the three absolute signed sums is convex in each
coordinate. Its maximum is therefore achieved at a binary vertex.
For vertex weights zero through four, the averages are

```text
0, 1, 2/3, 1, 0.
```

Consequently the average logarithmic primitive norm load of the
three words is at most `W`. If all three words have nonunit
primitive numerators, some word therefore has norm load at most
`W` and is suitable for the square-root separation argument.

### Repair when the averaged word is a unit

At most one of the three balanced words can be a unit. Indeed every
word has even total unit parity by the assumed label relation. A
unit word therefore has phase `+1` or `-1`. Its phase angle has
absolute value at most twice the arc width, which is less than two
because `C<1` and `R>=1`; hence its phase must be `+1`. Two such
balanced multiplicative identities would force two of the four
points to agree, since all their angles lie in an interval of
width less than one.

Suppose, after relabeling, the first word is a unit. Then at every
prime

```text
a_1+a_2=a_3+a_4.
```

The other two signed loads are `2(a_1-a_4)` and `2(a_1-a_3)`.
For allocations in `[0,e]` satisfying the displayed relation,

```text
|a_1-a_4|+|a_1-a_3|<=e.                         (3)
```

If `a_1` lies between `a_3,a_4`, the left side is
`|a_3-a_4|`. Otherwise it equals `|a_1-a_2|`. Both are at most
`e`. Thus the two remaining logarithmic norm loads sum to at most
`2W`. Neither remaining word is a unit, so one has nonunit
primitive numerator `g` with `log N(g)<=W`.

For every chosen sign partition, the parity relation forces all
signed prime exponents to be even and the total Gaussian unit
ratio to be `+/-1`. Hence this nonunit `g` has square norm.
A balanced four-point phase word has argument of absolute value
at most twice the arc width. Therefore

```text
R^(-1/2) <= |Im g|/|g| <= C R^(-1/2),
```

contrary to `C<1`.

A relation supported on three nonanchor labels is handled by adding
the anchor as the fourth distinct point. Supports one and two are
excluded by the pairwise argument. This proves the four-wise claim.

An exact allocation check covered 6,048 four-row weighted instances
in dimensions one through three, including 642 with a zero balanced
word. It verified the existence of a nonzero word of load at most
`W` and the repair inequality (3). The proof above is unrestricted.

### The resulting packing bound

Products over subsets of at most two of the `M-1` nonanchor labels
are all different. Thus, for augmented label rank `d`,

```text
1+(M-1)+binom(M-1,2)=1+binom(M,2)<=2^d.          (4)
```

In particular `M=O(2^(d/2))`. This is the appropriate elementary
baseline for comparing the six-wise theorem. For an arbitrary fixed
arc constant, subdividing into a fixed number of arcs with constant
less than one gives the same uniform-in-radius conclusion with a
factor depending only on that constant.

## 4. What six-wise separation adds, and where it gives uniformity

For six balanced signed rows, the worst binary cut has three rows
on either side. Its average absolute load is `6/5`, so the preceding
unconditional averaging argument no longer gives load at most `W`.
The near-uniform profile instead averages over all cuts and produces
`15/16<1`. This is the genuine new input in the six-wise theorem.
It improves the packing bound for the extracted tuple to

```text
sum_(j=0..3) binom(m,j)<=2^d,                   (5)
```

where `m` is its number of nonanchor rows. Distinct nonzero labels
alone do not imply either (4) or (5); for example the three labels
`e_1,e_2,e_1+e_2` are distinct but dependent.

The following restricted classes have uniform bounds, but bounded
rank already supplies them through (2) or (4).

* **Bounded squareclass rank.** If the augmented rank is bounded by
  a fixed `d`, (4) bounds the cluster independently of `R`.
* **Bounded split-prime support.** If the primitive squared radius
  has at most `s` distinct split primes, then `d<=s+1`. Arbitrarily
  large prime exponents are allowed. This is an elementary parity
  consequence, not a new consequence of central truncation.
* **Even core with bounded correction support.** In an independent
  block system with every core exponent even, `N(A_i)` is a square.
  If the correcting norms use at most `s` rational primes, the
  conjugate-primitive ratios `h_i/bar(h_i)` have squareclass rank
  at most `s`. Whenever the reconstructed arc constant is below
  one, (4) already gives `1+binom(m+1,2)<=2^s`. The six-wise
  hypotheses strengthen this to (5), but are not what first makes
  the count uniform. Bounded correction support imposed on growing
  central extractions likewise rules out unbounded clusters, already
  by this elementary argument.

The new restriction does not exclude the squarefree full independent
block model: any prime in each singleton block isolates that row's
squareclass coordinate, making all row labels linearly independent.
In contrast, even core exponents push all these classes into the
correcting factors, where (5) gives a genuine lower bound on their
required support or rank.

A perfect-power squared radius by itself does not bound that rank.
The original exponents `e_p` may all be divisible by a large fixed
integer while the intermediate allocations `a_ip` still take odd
values independently at arbitrarily many primes. Nor does a small
logarithmic correction height bound the number of available primes
independently of the radius. These are the missing upper bounds if
one seeks unrestricted uniformity from the packing inequalities.

## 5. The exact group for higher powers

Let

```text
U={u in Q(i)^*: u bar(u)=1}.
```

Unique factorization gives

```text
U = mu_4 times direct_sum_(p=1 mod4) Z,
u = i^a product_p (pi_p/bar(pi_p))^r_p.
```

For an integer `ell>=2`, the correct group of power classes is

```text
U/U^ell = Z/gcd(4,ell)Z
          plus (direct_sum_(p=1 mod4) Z/ell Z).     (6)
```

A relation in this group asserts an actual `ell`-th power of a
rational norm-one number. The torsion coordinate is part of the
statement. On the common-unit tuples used in central truncation it
is zero, so the condition reduces to the signed allocation vectors

```text
v_i=(a_ip-a_0p mod ell)_p,
sum_i c_i v_i=0.
```

For those common-unit tuples, the canonical primitive numerator of
the corresponding phase word is literally `gamma^ell`: every
signed prime exponent is divisible by `ell`, and the chosen unit
is one. This is the clean setting for the higher-power phase theorem.

For `ell=2`, (6) is exactly the augmented norm squareclass map in
(1): the unit generator `i` maps to the prime-two class, and
`pi/bar(pi)` maps to the prime-`p` class. For higher powers,
ordinary primitive numerator norms do **not** provide this map.
Already inversion changes a signed exponent from `1` to `-1` but
leaves its absolute norm exponent equal to one. Their difference
is invisible modulo two, but not modulo four.

There is a concrete two-row failure even if one tries to choose
balanced signs afterward. For two distinct split primes `pi,rho`,
take

```text
h_1=pi rho,
h_2=pi^3 bar(rho)^3.
```

Then `N(h_1)N(h_2)=(pq)^4` is a fourth power, but their signed
vectors modulo four are `(1,1)` and `(3,-3)`. Neither their sum
`(4,-2)` nor their difference `(-2,4)` is zero modulo four.
Thus a norm fourth-power relation does not supply the required
signed Gaussian relation, even after either choice of relative
sign. These are exact arithmetic data, not an endpoint example.

Finally, groups modulo a composite `ell` are not vector spaces.
The safe packing statement counts the finite generated group's
order. For example a subgroup of `(Z/4)^s` has order at most
`4^s`, whereas its free rank alone need not determine its order.
Equal-cardinality subset-sum collisions supply balanced signed
relations, so a theorem excluding such relations yields
`|G_ell|>=binom(m,r)` for the applicable subset size `r`.
It does not automatically assert independence of every support
up to twice that size with arbitrary signs.
