# Kurlberg--Wigman measures and the shrinking-arc gap

## Scope

This note audits Kurlberg--Wigman, [*On probability measures arising from
lattice points on circles*](https://arxiv.org/abs/1501.01995), specifically
for a possible bound on lattice points in an angular interval of length
`N^(-1/4)`.  The paper supplies an exact squarefree convolution model, but
it does **not** supply a quantitative small-ball or inverse-concentration
estimate at a scale that shrinks with `N`.

The useful output below is a finite divisor reduction. An unbounded
endpoint cluster, with arbitrary prime powers, would force quadratically
many distinct Gaussian divisors of norm just above `sqrt(N)` to have
subpower imaginary part. The reduction is proved in Section 4; the
additional divisor upper bound needed for uniformity remains unproved.

## 1. What the paper actually proves

For a split prime `p`, Section 4.1 writes the de-symmetrized prime measure
as

```text
nu_p = (delta_(theta_p) + delta_(-theta_p))/2.
```

For coprime factors these measures convolve. For odd squarefree

```text
N = product_(j=1)^k p_j,
```

then

```text
nu_N = convolution_j (delta_(theta_j)+delta_(-theta_j))/2,
hat(nu_N)(h) = product_j cos(h theta_j),
r_2(N) = 4*2^k.                                      (1)
```

Equivalently, `nu_N` is the pushforward of the uniform measure on the
Rademacher cube `{+1,-1}^k` by

```text
epsilon |-> sum_j epsilon_j theta_j  modulo 2 pi.      (2)
```

Theorem 1.4 classifies the possible pairs `hat(mu)(4),hat(mu)(8)`
for squarefree weak limits. The convolution closure argument controls
fixed Fourier modes. Lemma 4.1 uses Hecke equidistribution to locate
prime angles in fixed neighborhoods. These statements give no
shrinking-scale error estimate. The recalled Cilleruelo sequences
converge to four axis atoms while their representation counts grow;
their concentration rate is not compared with `N^(-1/4)`.

## 2. The exact weak-convergence gap

Let `I_N` be an angular interval of length `Delta_N=N^(-1/4)` containing
`M_N` lattice points.  Its normalized mass is

```text
mu_N(I_N) = M_N/r_2(N).                                (3)
```

A uniform endpoint bound asks for the much stronger estimate

```text
mu_N(I_N) = O(1/r_2(N)).                               (4)
```

Weak convergence controls fixed continuity sets.  It neither controls the
moving sets `I_N` nor resolves masses of size `1/r_2(N)`.  Even after the
centers of `I_N` are made convergent, a cluster with

```text
M_N -> infinity,             M_N=o(r_2(N))
```

has mass tending to zero and is invisible in every weak limit.  If the
limit has an atom, even an `o(1)` small-ball conclusion is unavailable.
Proposition 1.7, which permits `r_2(N)->infinity` along sequences attaining
any attainable measure, does not close this gap.

There is the same obstruction on the Fourier side.  A trigonometric
majorant for an interval of length `Delta_N` needs frequencies through

```text
K asymp Delta_N^(-1) = N^(1/4).                        (5)
```

By (1), this would require simultaneous information on
`product_j cos(h theta_j)` for `h<=K`.  Fixed Fourier coefficients, including
the two coefficients classified in Theorem 1.4, give no such information.
Even convergence of every fixed Fourier coefficient is compatible with an
unbounded cluster of vanishing relative mass.

## 3. A direct sign-cube inverse statement

It is clearer here to undo the paper's fourfold angular
de-symmetrization.  Choose Gaussian primes

```text
pi_j = sqrt(p_j) exp(i phi_j).
```

In one fixed Gaussian-unit class, the `2^k` points are

```text
z_epsilon = product_j pi_j^((1+epsilon_j)/2)
                       conjugate(pi_j)^((1-epsilon_j)/2),
epsilon in {+1,-1}^k.                                  (6)
```

Their physical angles are `sum_j epsilon_j phi_j`.  A general cluster can
be partitioned into four such unit classes, costing only a factor four.

For `epsilon != epsilon'`, put

```text
S(epsilon,epsilon') = {j: epsilon_j != epsilon'_j},
P(epsilon,epsilon') = product_(j in S) p_j,
A(epsilon,epsilon') = product_(j in S, epsilon_j=+1) pi_j
                      product_(j in S, epsilon_j=-1) conjugate(pi_j).
                                                               (7)
```

Then `Norm(A)=P` and

```text
arg(z_epsilon/z_epsilon') = 2 arg A  modulo 2 pi.       (8)
```

The Gaussian integer `A` has nonzero imaginary part.  Otherwise unique
factorization would make a nonempty product selecting exactly one prime
above each `p_j` associate to its conjugate, which is impossible.
Consequently, if the two points have angular distance at most `Delta<pi`,

```text
1 <= |Im A| <= sqrt(P) sin(Delta/2),
P >= csc(Delta/2)^2 >= 4/Delta^2.                       (9)
```

At the endpoint scale `Delta=C N^(-1/4)`, every pair in the cluster thus
differs on a prime support satisfying

```text
P >= (4/C^2) sqrt(N).                                  (10)
```

This is the prime-support form of the weighted half-distance condition.
It is stronger than a statement about the limiting measure, but by itself
it is only the critical Plotkin threshold and does not bound the code size.

There is also an exact local inverse-energy identity.  For a fixed cube
point `epsilon`, its neighbors within angular distance `Delta` are in
bijection with the nonempty subsets `S` for which

```text
dist(2 arg A_(epsilon,S), 2 pi Z) <= Delta,             (11)
```

where `A_(epsilon,S)` is the oriented product in (7).  Hence bounding the
number of short-arc neighbors is exactly a near-real Gaussian-divisor
problem, not a weak-limit problem.

Equivalently, summing over all ordered pairs gives the ternary relation
energy

```text
E_Delta = sum_(eta in {-1,0,1}^k, eta != 0,
                 dist(2 sum_j eta_j phi_j,2 pi Z)<=Delta)
                  2^(#{j:eta_j=0}).                    (12)
```

An interval containing `M` sign-cube points satisfies

```text
M(M-1) <= E_Delta.                                     (13)
```

The multiplicity in (12) is exact: a fixed ternary difference `eta`
comes from two forced signs on its support and two common choices at every
zero coordinate.  A useful inverse theorem must retain this multiplicity;
merely finding one short signed sum is insufficient.

## 4. A general concentration lemma for near-critical supports

The following consequence no longer requires squarefreeness or a binary
sign cube.  Start with a fixed-unit configuration and divide its common
Gaussian divisor.  This decreases its normalized arc constant, so it is
harmless for an upper-bound argument.  For the resulting primitive
configuration, inert factors and the ramified factor have disappeared;
write

```text
N = product_p p^(s_p),
z_a = product_p pi_p^(e_(a,p)) conjugate(pi_p)^(s_p-e_(a,p)),
0 <= e_(a,p) <= s_p.                                   (14)
```

Primitivity says that, for every remaining `p`, the exponents
`e_(a,p)` attain both zero and `s_p` somewhere in the configuration.
For an ordered pair define

```text
A_ab = product_(p:e_(a,p)>e_(b,p)) pi_p^(e_(a,p)-e_(b,p))
       product_(p:e_(a,p)<e_(b,p))
                    conjugate(pi_p)^(e_(b,p)-e_(a,p)),
P_ab = Norm(A_ab) = product_p p^|e_(a,p)-e_(b,p)|.      (15)
```

Then `P_ab` divides `N`, and unique factorization gives

```text
z_a/z_b = A_ab/conjugate(A_ab).                         (16)
```

The integer `A_ab` is nonreal: it is nontrivial and at each rational split
prime uses only one of the two conjugate Gaussian primes.  Thus the
integrality estimate (9) applies unchanged.

The result below combines this exponent model with the angular
integrality and short-product rigidity already proved for circle points.

**Lemma.**  Assume `0<C<=2`.  If one Gaussian-unit class contains `M>=3`
points in an interval of angular length at most `C N^(-1/4)`, then there
are at least

```text
ceil(M(M-2)/4)                                          (17)
```

distinct Gaussian associate-and-conjugacy classes represented by
integers `A_ab` of the form (15) satisfying

```text
(4/C^2) N^(1/2) <= Norm(A_ab) <= N^(1/2+1/M),           (18)
1 <= |Im A_ab| <= (C/2) N^(1/(2M)).                     (19)
```

**Proof.**  Put `W=log N` and, for the `M` exponent rows in the interval,
let

```text
d_ab = sum_p |e_(a,p)-e_(b,p)| log p = log P_ab.
```

For a fixed prime, the layer-cake identity is

```text
sum_(a<b) |e_(a,p)-e_(b,p)|
  = sum_(t=1)^s_p n_(p,t)(M-n_(p,t)),                   (20)
```

where `n_(p,t)=#{a:e_(a,p)>=t}`.  Each summand is at most `M^2/4`.
Multiplying by `log p` and summing over the `s_p` layers gives

```text
sum_(a<b) d_ab <= W M^2/4.                              (21)
```

The points of each pair have angular distance at most
`Delta=C N^(-1/4)`.  The three-point chord argument in (26) below gives
`Delta<pi/2`, so applying (9) to (15)--(16) gives

```text
P_ab >= (4/C^2) sqrt(N),
d_ab >= W/2+log(4/C^2) >= W/2                           (22)
```

when `C<=2`.  Subtracting `W/2` for every pair from (21) gives

```text
sum_(a<b) (d_ab-W/2) <= WM/4.                           (23)
```

Fewer than `M^2/4` pairs can therefore have
`d_ab>W/2+W/M`.  At least

```text
binom(M,2)-M^2/4 = M(M-2)/4                            (24)
```

pairs obey the reverse inequality.  For each such pair, (15) and (22) give
(18), while

```text
|Im A_ab| <= (C/2) N^(-1/4) sqrt(Norm(A_ab))
            <= (C/2) N^(1/(2M)),                       (25)
```

which is (19).

It remains to show that two pairs counted in (24) cannot yield the same
associate-and-conjugacy class.  First note that the cluster's angular
width is less than `pi/2`.  Indeed, order any three of its points along
the arc.  The chord between consecutive distinct equal-norm Gaussian
integers has length at least `sqrt(2)`, so the arc length is at least
`2 sqrt(2)`.  Since it is also at most `C sqrt(R)`, where `N=R^2`,

```text
sqrt(R) >= 2 sqrt(2)/C,
Delta <= C/sqrt(R) <= C^2/(2 sqrt(2)) <= sqrt(2) < pi/2. (26)
```

For an ordered pair `(a,b)`, (16) is the exact quotient.  Conjugating
`A_ab` reverses the ordered pair, while changing its sign does not change
(16).  If two pair classes agree under these operations, orient
the pairs so that their quotients agree.  Then

```text
z_a z_d = z_c z_b.                                      (27)
```

The `B_2` product-rigidity theorem in
`angular_moment_product_rigidity.md` applies because `2C^2<=8`.  It says
that the two product multisets in (27) are equal.  Since the points in
each pair are distinct, this forces the same unordered pair.

Multiplication of `A` by `i` negates the quotient in (16).  It cannot
identify two classes arising here: both point quotients have unwrapped
arguments in `[-Delta,Delta]`, and (26) excludes a difference of `pi`.
Equivalently, both `A` representatives lie within `Delta/2<pi/4` of the
real axis modulo `pi`, whereas multiplication by `i` moves one to within
`Delta/2` of the imaginary axis.  Thus every pair in (24) supplies a
different full associate-and-conjugacy class, proving (17).  `square`

Thus a hypothetical sequence with `M->infinity` forces an unbounded
collection of oriented Gaussian divisors with norm
`N^(1/2+o(1))` and imaginary part `N^o(1)`.  This conclusion uses both
prime support and concentration of many code pairs; it is not produced by
the convolution-limit theorem alone.

## 5. The remaining arithmetic statement

A sufficient next lemma would uniformly bound the number of
associate-and-conjugacy classes of the conjugate-coprime Gaussian divisors
`A_ab` in (15).  At each rational prime they select one of the two Gaussian
primes to an exponent at most `s_p`, and their rational norms divide `N`.
The relevant region is

```text
Norm(A) in [(4/C^2) sqrt(N), N^(1/2+o(1))],
0 < |Im A| <= N^o(1).                                  (28)
```

The `o(1)` rates in (28) are linked: the lemma above gives exponents
`1/M` and `1/(2M)` for an `M`-point cluster.  Proving such a bound for
arbitrarily slow `M->infinity` would rule out unbounded primitive clusters,
up to the four unit classes.  For a cluster containing several unit
labels, first retain a largest unit class and then divide that subcluster's
common Gaussian divisor; this costs a factor at most four in cardinality
and only improves the arc constant.  Section 4 already includes arbitrary
nested prime exponents `0<=e_(a,p)<=s_p`; no binary or two-level assumption
is used in its layer-cake proof.

No result in Kurlberg--Wigman implies this divisor statement.  Prime-angle
equidistribution concerns individual primes at fixed angular resolution;
(28) concerns correlated products of exponent differences at a resolution
determined by the full product `N`.  Standard Littlewood--Offord bounds
also need a nonconcentration or dissociation hypothesis on the angles.  Here angles
may be arbitrarily close to the axes, and modulo-circle signed sums can
resonate.  The Gaussian integrality in (9) is the available replacement,
but at present it yields the critical half-support threshold and the
near-real-divisor reduction rather than a uniform count.

Accordingly, the paper offers a useful exact coordinate system for a new
attack, not an endpoint estimate: one must prove a quantitative inverse
theorem for the oriented divisor products in (28), or an equivalent bound
for the local degrees in (11).  Weak convergence and its classified fixed
Fourier projections are insufficient.

The [exact arithmetic checker](check_gaussian_near_real_divisor_reduction.py)
verifies 25,368 layer-cake cases, quotient identities for an actual tuple
with three and four valuation levels, and three members of the four-point
Pell family. For the Pell examples the endpoint certificate `C<2` uses
only rational half-angle comparisons. The admissible divisor-class counts
are five, six, and six. These tests support the general prose proof;
they do not test or assume the missing divisor upper bound.

For M>=3, the [ordinary-norm reduction](gaussian_strip_divisor_route_audit.md)
shows more: all `binom(M,2)` pairs have distinct rational norms d, because
`|Im A|^2<=C^2d/(4sqrt(N))<=sqrt(d)` forces
`|Re A|=floor(sqrt(d))`. Thus each d is a divisor of N whose difference
from the square immediately below it is a small positive square.
At least `ceil(M(M-2)/4)` lie in the near-critical interval (18).
The same note explains why the
known Gaussian residue-class algorithm cannot directly group these
divisors into a bounded number of classes.
