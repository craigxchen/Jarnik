# Prime incidences of the primitive chord residues

This note keeps actual Gaussian integrality while examining residue
primes outside the conductor. It gives an exact congruence hierarchy,
finite collision counts, and a genuine bounded-residue obstruction
family. It does not prove a bound uniform in the radius.

The associated global inert-prime estimate, including its incidental
improvement of the leading sublogarithmic constant, is recorded
separately in `inert_prime_cofactor_bound.md`. That improvement does
not change the growth rate sought in the main problem.

## 1. Outside-prime residues are exact congruence depths

Let `z_1,...,z_M` be distinct Gaussian integers of common norm `N`.
For each pair choose a Gaussian gcd and write

```text
g_ij = gcd_G(z_i,z_j),
c_ij = (z_i-z_j)/g_ij in Z[i].
```

No common-unit restriction is needed. If an odd rational prime `p`
does not divide `N`, then `g_ij` is a unit modulo every power of `p`.
Consequently

```text
p^a divides c_ij in Z[i]
    if and only if z_i is congruent to z_j modulo p^a.        (1)
```

In the common-unit primitive-chord convention, `|c_ij|=2|t_ij|`,
so for odd `p` the valuation in (1) is precisely `v_p(t_ij)`.
It is safer to use the Gaussian cofactors `c_ij` when different unit
classes occur; no factor-four selection is then introduced.

At every fixed depth `a`, the divisible pairs form disjoint complete
graphs, one for each occupied residue class. As `a` increases these
classes refine. Equivalently the pairwise congruence depths satisfy
the ultrametric inequality

```text
v_p(c_ij) >= min(v_p(c_ik),v_p(c_kj)).             (2)
```

Here `v_p` of a Gaussian number means the largest rational power
`p^a` dividing both coordinates. The minimum of the three depths in
a triangle is attained at least twice. Hence arbitrary independent
edge assignments of residue primes are not permissible: a prime
dividing two edges of a triangle must divide the third.

The parallel global-content investigation strengthens the analogous
statement for the integer primitive residues to conductor primes as
well, using primitive half-angle numerators. The congruence identity
(1) itself is asserted here only outside `N`.

## 2. The exact number of available conic classes

For odd `p` not dividing `N`, the number of solutions of
`x^2+y^2=N` modulo `p^a` is

```text
q_(p,a) = (p-chi_p(-1)) p^(a-1).                  (3)
```

For split `p`, choose a square root `i` of `-1` in the finite field.
The equation becomes `(x+i y)(x-i y)=N`, with `p-1` choices of the
first nonzero factor. For inert `p`, the equation is a nonzero norm
fiber in `F_(p^2)`, with `p+1` elements. Every solution is smooth
because the gradient `(2x,2y)` is nonzero, and each solution lifts in
exactly `p` ways at each further level. This proves (3).

Put

```text
E(M,q) = q binom(floor(M/q),2)
         +(M mod q) floor(M/q).
```

Distributing `M` objects among at most `q` classes creates at least
`E(M,q)` coincident pairs: the sum of binomial class counts is
minimized when the occupancies differ by at most one. Summing (1)
over levels gives the exact finite inequality

```text
sum_(i<j) v_p(c_ij)
 >= sum_(a>=1) E(M,(p-chi_p(-1))p^(a-1)).         (4)
```

The sum is finite because `E(M,q)=0` when `q>=M`.

For a primitive Gaussian cluster, every inert prime is outside the
conductor. Thus all primes `p=3 mod 4` can be used simultaneously,
without an assumption about missing conductor primes. At inert primes,
each divisibility by `p` costs modulus `p`, rather than `sqrt(p)`.
The exact finite product bound is therefore

```text
product_(i<j) |c_ij|
 >= product_(p=3 mod 4) p^[sum_(a>=1) E(M,(p+1)p^(a-1))].     (5)
```

This is genuine information about the collective residues, stronger
than applying only their individual nonzero lower bounds. The right
side nevertheless depends on the number of points, not on the radius.
For each fixed tuple size it is a fixed constant, compatible with
residues that grow subpolynomially with the radius.

The asymptotic conversion of (5) uses the arithmetic-progression
Mertens theorem; a primary reference is Kenneth S. Williams,
[Mertens' Theorem for Arithmetic Progressions](https://people.math.carleton.ca/~williams/papers/pdf/057.pdf),
Theorem 1 and Section 5. The full derivation and the all-unit Gaussian
gcd upper bound are in `inert_prime_cofactor_bound.md`.

## 3. A genuine arbitrarily large-radius family with bounded residues

Small residues and bounded common content alone do not bound the
radius, even if all Gaussian reconstruction and quadruple identities
are imposed. For any fixed `M>=2`, let the positive integer `n` grow
and define

```text
H_j = 2(n+j)+i,               0<=j<M,
z_j = H_j product_(l!=j) conjugate(H_l).
```

These are distinct Gaussian integers on the circle of radius

```text
R_n = product_j |H_j| ~ (2n)^M.
```

Every `H_j` is coprime to its conjugate: a common divisor divides
`2i`, while `N(H_j)=4(n+j)^2+1` is odd. For a pair, put
`g=gcd_G(H_i,H_j)`. The exact reduced numerator for the ratio
`z_i/z_j` is

```text
h_ij = H_i conjugate(H_j)/N(g),
z_i/z_j = h_ij/conjugate(h_ij).
```

Indeed the common Gaussian divisor of `H_i conjugate(H_j)` and
its conjugate is the rational integer `N(g)`. The remaining factors
are conjugate-coprime. Taking the imaginary part gives the exact
primitive residue

```text
|t_ij| = 2|i-j|/N(g) <= 2(M-1).                  (6)
```

There is no divisibility assumption hidden in this quotient. The
Gaussian gcd `g` divides the rational integer `2(i-j)` and is
conjugate-coprime, so its rational norm divides `2(i-j)`.

For adjacent indices, `N(g)` is odd and divides `2`, so it equals
one. Thus adjacent residues are exactly two. All residues in (6)
are even, and their global gcd is exactly two. Equivalently the
Gaussian cofactors have moduli at most `4(M-1)`, with adjacent
cofactor modulus four and global rational content four.

In particular all primes appearing in these primitive residues are
bounded in terms of `M`, even as the circle radius diverges. Their
outside-prime congruence partitions satisfy (1)--(4) automatically.
Every Gaussian metric identity and every rational quadruple identity
also holds because these are actual Gaussian points.

### Removing a common Gaussian factor does not remove the example

The common Gaussian gcd of all the `z_j` is bounded in terms of `M`,
independently of `n`. A prime common to every `z_j` must occur in
at least two of the factors among the `H_j` and their conjugates.
Any such repeated divisor divides one of the fixed differences

```text
H_i-H_j = 2(i-j),
H_i-conjugate(H_j) = 2(i-j)+2i.
```

Their valuations bound all pairwise overlaps. At one Gaussian prime,
the common exponent in the `z_j` is at most the sum of all but the
largest exponent among the `H_j` and their orientations; that sum
is bounded by the sum of their pairwise overlap exponents. Thus the
common factor is bounded by a fixed product of the displayed
differences and their conjugates. The primitive radius still has
order `n^M`, and the pair residue ratios are unchanged.

### The missing hypotheses are visible in this family

The angular span is

```text
Delta_n = 2[arctan(1/(2n))-arctan(1/(2(n+M-1)))]
        ~ (M-1)/n^2.
```

Hence the normalized endpoint length has order

```text
Delta_n sqrt(R_n) ~ (M-1) 2^(M/2) n^(M/2-2).
```

For `M>=5` it diverges. This family therefore does not give an
endpoint counterexample. Its nonconstant conductor mass is
asymptotically supported on singleton cuts: apart from bounded
pairwise overlaps, a prime from `H_j` distinguishes only row `j`.
It is far from the extracted uniform profile.

The example rules out a radius or growing-common-content assertion
based only on bounded primitive residues and Gaussian reconstruction.
It does not rule out a theorem using the near-uniform conductor
profile together with the endpoint span. Its rational projective shape
is fixed, so it is also consistent with the new fixed-template theorem.

## 4. What would still produce a uniform bound

The extracted profile permits all primitive residues to be
`R^o(1)`. The congruence hierarchy controls how their prime powers
can be distributed, and (5) gives a nontrivial cost in the tuple size.
Neither result gives a positive power of `R` for a fixed tuple.

A sufficient genuinely new input would be a statement that, for some
fixed tuple size and a fixed neighborhood of the uniform conductor
profile, the joint residue hierarchy and Gaussian metric equations
force

```text
max_(i<j) |t_ij| >= R^delta
```

with a fixed `delta>0`, or force a projective or multiplicative
degeneracy already excluded by the endpoint geometry. No such
statement is proved here. The bounded-residue family shows why the
profile hypothesis cannot be omitted from this next step.

## Verification

Across 1,035 actual-circle configurations with squared radius at most
1,000, including mixed Gaussian-unit classes, exact arithmetic checked
the pair-gcd product bound

```text
product_(i<j) N(g_ij) >= N^[binom(M,2)-floor(M^2/4)].
```

It also checked the equivalent cofactor upper bound with the largest
actual squared chord, and 7,859 outside-prime collision bounds from
(4). These are finite sanity checks of the stated normalization;
the proofs do not depend on their range. None of these new residue
results is asserted to be Lean-formalized.

For the bounded-residue family, 28,600 exact pair reductions in 1,100
configurations (`2<=M<=12`, `1<=n<=100`) checked formula (6),
conjugate-coprimality of the reduced numerators, adjacent residues two,
global residue gcd two, and global Gaussian-cofactor content four.

## 5. Exact feasible residue partitions at a conductor prime

There is a sharper classification than ultrametricity alone at a
split conductor prime. Let `p^e` divide `N` exactly, choose an
orientation `pi` above `p`, and let

```text
a_i = v_pi(z_i),       v_conjugate(pi)(z_i)=e-a_i,
I_a = {i:a_i=a},       n_a=|I_a|.
```

Use the Gaussian cofactors `c_ij`; the conclusion permits arbitrary
unit classes. The following are exact.

1. If `a_i!=a_j`, then `v_p(c_ij)=0`.
2. Within one allocation group `I_a`, put
   `w_i=z_i/[pi^a conjugate(pi)^(e-a)]`. All `w_i` are Gaussian
   integers of the same norm `nu=N/p^e`, which is a `p`-unit. Then

   ```text
   v_p(c_ij)=v_p(w_i-w_j)       (i,j in I_a).       (7)
   ```

For the first statement, divide the difference by the pair gcd.
At `pi`, exactly one of the two reduced terms is a unit and the
other is divisible by `pi`; the conjugate place gives the analogous
conclusion. For the second, the displayed common factor removes the
entire `p`-part of the pair gcd, and any remaining pair gcd is a
unit at `p`.

Consequently each allocation group has at most `p-1` residue classes
at depth one, and every class has at most `p` children at each later
depth. Residue edges cannot join two different allocation groups.
The exact finite lower bound is

```text
sum_(i<j) v_p(c_ij)
 >= sum_(a=0..e) sum_(h>=1) E(n_a,(p-1)p^(h-1)).  (8)
```

The multiplicity `e` enters through the allocation groups, not as a
multiplier of the collision cost. In particular, raising the exponent
of a prime while keeping the rows at allocations `0` and `e` does
not increase the right side of (8).

### The branching conditions are also sufficient locally

For a fixed `p`, `e`, and unit `nu`, these conditions describe all
possible finite residue-depth partitions of a local Gaussian norm
configuration. This is a statement over `Z_p[i]`, not an assertion
of a simultaneous global integer realization.

Because `p` splits,

```text
Z_p[i] = Z_p x Z_p,
conjugation swaps the factors,
norm(u,v)=uv.
```

The unit norm-`nu` circle is parametrized exactly by

```text
w(u)=(u,nu/u),       u in Z_p^*.
```

Its congruence depth satisfies

```text
v_p(w(u)-w(v))=v_p(u-v),
```

since `nu/u-nu/v=nu(v-u)/(uv)` and all the extra factors are units.
The units have `p-1` first digits and `p` choices at every subsequent
digit. Thus any finite hierarchy with precisely the stated branching
limits can be realized by selecting `p`-adic digits of the parameters
`u_i`. Multiplying each allocation group by
`pi^a conjugate(pi)^(e-a)` realizes its prescribed allocation.

All these local points have norm `p^e nu`. Their metric matrices are
actual rank-one Hermitian matrices and their cocycles and quadruple
identities hold exactly. Hence those local Gaussian identities do not
further restrict the finite partition types classified above. They
can still impose global restrictions when combined with size bounds
and simultaneous integrality at different places.

## 6. The combined weighted bound and its fixed-profile limitation

Combining (4) and (8) gives, for a primitive cluster, the explicit
finite lower bound

```text
sum_(i<j) log |c_ij|
 >= sum_(odd p not dividing N) log(p)
      sum_(h>=1) E(M,(p-chi_p(-1))p^(h-1))
  + sum_(split p^e exactly dividing N) log(p)
      sum_(a=0..e) sum_(h>=1) E(n_a,(p-1)p^(h-1)).   (9)
```

Only finitely many terms on the right are nonzero. This is valid
because a rational prime power dividing a Gaussian cofactor costs
its full rational modulus, and distinct rational prime powers can be
multiplied. The prime two is simply omitted, preserving a valid
lower bound.

The extra conductor information in (9) does not exclude any fixed
neighborhood of the uniform conductor profile by itself.

* If `p>M`, each allocation group can occupy distinct first-level
  unit classes. Its entire forced residue contribution in (8) is zero.
  All conductor weight may come from such primes without violating
  any of these local partition conditions.
* If a fixed prime is raised to an arbitrarily large exponent and
  its two endpoint allocation groups stay fixed, the conductor weight
  is `e log p`, while the forced residue contribution remains unchanged.
  Thus no lower bound proportional to that conductor weight follows
  from the partition count.
* For every fixed tuple size, the forced outside-prime cost is a fixed
  number independent of `N`. Allowing `R^o(1)` residues accommodates
  such a cost; a fixed tuple cannot turn it into `R^delta`.

For example, in the squarefree uniform-pattern model, choose all
conductor primes larger than the fixed tuple size. At each of these
places choose different unit parameters within each allocation group.
Every conductor-prime residue depth is then zero. At the finitely
many smaller outside primes, choose any allowed finite hierarchy,
with depths bounded in terms of the tuple size. Each local norm
configuration exists by the explicit parametrization above, for the
same prescribed integer norm `N`.

This produces compatible **local partition data and local Gaussian
metric identities**, not a lattice-point counterexample. Gluing the
local points to actual Gaussian integers on the specified short arc
is precisely the unresolved global size problem. The classification
shows that a new fixed-profile exclusion cannot rest only on the
number, branching, and nesting of residue classes at individual primes.

There is a related bookkeeping caution. If
`T=product_(i<j)|t_ij|` has subpower height, the conductor primes
dividing `T` have radical log weight at most `log T=o(log R)`.
For unrestricted prime powers this does **not**, by itself, bound
their full conductor weight `sum e_p log p`: a prime can divide a
same-allocation residue to low order while its exponent in `N` is
large. Discarding those primes as negligible conductor mass requires
an additional argument, or a squarefree-conductor hypothesis.

For this extension, exact arithmetic checked 45,872 pair classifications
and 32 grouped prime-power occupancy bounds on the sixteen full,
mixed-unit circles of norms `5^e 13^f`, with `1<=e,f<=4`.
