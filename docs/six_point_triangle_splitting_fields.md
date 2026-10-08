# Exact splitting fields of the ten triangle quadratics

Fix nonzero rational weights `u_0,...,u_5` with total sum zero, and assume
`beta=diag(u)` is hyperbolic over `Q`. Fix a rational hyperbolic frame
containing `e=1`, and either of the two rational isotropic pencils from
[`six_point_isotropic_circle_cover.md`](six_point_isotropic_circle_cover.md).
For a triple `I`, put `J=I^c` and assume

```text
S=sum_(i in I)u_i != 0,
D_I=-S product_(i in I)u_i.
```

Let `M_I(s,t)` be the actual determinant of the three pencil columns in
`I`, with their common clearing retained. Then `M_I` is a nonzero
squarefree binary quadratic on either ruling and

```text
disc(M_I) = D_I times a nonzero rational square.     (1)
```

This assertion for a selected partition needs no condition on the other
proper subsums. Full nonresonance ensures it applies to all ten partitions
and, separately, that roots from different partitions are disjoint.

## 1. The two binary spaces supported on the triples

Let `E_I` consist of vectors supported on `I` satisfying
`sum_I u_i x_i=0`; define `E_J` likewise. They are orthogonal binary
spaces. Their restrictions are nondegenerate because the weights and `S`
are nonzero. Their direct sum gives an isometric model of

```text
e^perp/<e> = E_I orthogonal_sum E_J.                 (2)
```

Indeed it has dimension four and lies in `e^perp`; its intersection with
`<e>` is zero since a nonzero constant vector has weighted sum `S` on
`I`. Thus projection gives the asserted isomorphism.

For weights `a,b,c` on `I`, use the basis `(b,-a,0),(c,0,-a)` of `E_I`.
The Gram determinant is `a^3 b c(a+b+c)`. The ordinary discriminant of
the binary quadratic has squareclass the negative of this determinant,
namely `-abc(a+b+c)=D_I`. Its two geometric isotropic lines are distinct
and are defined over `Q(sqrt(D_I))`; if `D_I` is a square they are both
rational. Similarly

```text
D_J=S product_(j in J)u_j,
D_I D_J=-S^2 product_(i=0)^5 u_i.                   (3)
```

Hyperbolicity of the six-dimensional form implies
`-product_i u_i` is a nonzero rational square, by taking the determinant
in a hyperbolic basis. Hence `D_I,D_J` have the same squareclass.

## 2. Four planes, two on each ruling

Each choice of an isotropic line `ell` in `E_I` and an isotropic line `m`
in `E_J` gives a maximal isotropic plane

```text
W=<e,ell,m>.                                       (4)
```

Its triple-`I` determinant vanishes because the row supported on `J`
restricts to zero on `I`. Conversely, if that determinant vanishes, the
restriction `W->Q^I` has a nonzero kernel. This is an isotropic line in
`E_J`: the binary space `E_J` has no two-dimensional isotropic subspace.
In the quotient (2), a maximal isotropic two-space containing `m` lies
in `E_I direct_sum m`, since `m^perp` within `E_J` is exactly `m`.
It therefore equals `ell direct_sum m` for an isotropic line `ell` of
`E_I`. This proves that the four choices in (4) are exactly all roots
across the two rulings.

Label the lines `ell_+,ell_-` and `m_+,m_-`. One ruling contains
`ell_+ direct_sum m_+` and `ell_- direct_sum m_-`; the other contains
`ell_+ direct_sum m_-` and `ell_- direct_sum m_+`. For a direct check,
choose generators `e1,e2` of `ell_+,m_+`, and rescale generators `f1,f2`
of `ell_-,m_-` so their cross-pairings are one. The first and second
pencils in (1) of the cover note have respectively these pairs as their
parameters zero and infinity. Each pencil's quadratic therefore has two
distinct geometric roots and cannot vanish identically.

If the common squareclass is nonsquare, its nontrivial Galois automorphism
switches both signs at once. It exchanges the two roots on each ruling
and does not interchange the rulings. The residue field is exactly the
quadratic field: from `W` one recovers `ell=W intersection Q^I`, so a
rational `W` could not conceal a nonrational isotropic line. If the
squareclass is square, all four planes and all four parameter roots are
rational. Since each chosen pencil is an isomorphism over `Q` with its
ruling, these field statements prove (1), including a rational root at
parameter infinity. In that case the homogeneous quadratic discriminant
is still nonzero, and is a square as required.

## 3. Exact good-prime consequence and its exceptions

Now take integral weights and a positive integer `h` clearing a normalized
hyperbolic frame. Use the raw pencil obtained by multiplying both moving
rows by `h`, so its minor forms are integral. For the selected partition
it is enough to exclude primes dividing

```text
E_I=2 h S product_i u_i.                            (5)
```

At any remaining prime the normalized frame reduces to a hyperbolic
basis: its exact Gram identities still hold modulo the prime. The two
binary spaces in (2) remain nondegenerate, and the same proof works over
`F_p`. Consequently

```text
M_I(s,t)=0 for some (s:t) in P^1(F_p)
  iff (D_I/p)=+1.                                  (6)
```

There are exactly two projective roots in the split case and none in the
nonsplit case. In particular, if an integral pair `(s,t)` is primitive
and a good prime divides `M_I(s,t)`, that prime splits in the indicated
quadratic field. Conversely every good split prime has two simple roots,
which lift by the elementary simple-root Hensel argument; there is no
upper bound on its possible valuation in primitive values of `M_I`.

The raw integer scaling is relevant. At a prime dividing `h`, a minor
may acquire content or the chosen pencil may lose rank modulo that prime.
At a prime dividing a weight or `S`, the binary restriction degenerates;
ramified primes of the quadratic field are among these exceptions (with
the prime two also excluded). Such primes can have roots, repeated roots,
or valuation phenomena not covered by (6). Formula (1) remains an exact
characteristic-zero squareclass identity, but it does not justify using
a nonzero Legendre symbol or bounding bad-prime valuations there.

For all ten partitions one may use the single exception integer
`2h product_i u_i product_I S_I`. With the controlled frame it has size
`U^O(1)`, but finite prime support alone gives no bound on valuations in
parameter values. Nor does the split-prime restriction itself bound
parameter height, radius, or the number of points in an endpoint arc.

## 4. Joint reciprocity does not obstruct the ten local cuts

Here is a precise simultaneous compatibility statement. Fix a fully
nonresonant `u` and a ruling with a nonempty real circle interval. Index
the ten partitions by the triples `I` contained in `{0,...,4}`. Then

```text
product_I D_I = (product_(i=0)^4 u_i)^6 product_I S_I. (7)
```

Thus the product squareclass is `product_I S_I`; it need not be trivial.
Any further relation among these squareclasses is automatically respected
by primes completely split in the multiquadratic field
`Q(i,sqrt(D_I) for all I)`.

Such primes can be chosen ten at a time, distinct and comparable. To avoid
needing a field-density theorem, take a fixed positive integer modulus
divisible by `8 product_I |D_I|`. A prime equal to one modulo this modulus
has `(D_I/p)=1` for every `I`, by quadratic reciprocity, and also has
`(-1/p)=1`. Enlarge the modulus by all nonempty proper subsums and the
finitely many frame, discriminant, and resultant exceptions so that all
raw minor roots are simple, distinct
across partitions, and disjoint from `Delta`. The fixed-modulus prime
number theorem gives at least ten such primes in `[Y,2Y]` for every
sufficiently large `Y`; see
[Sutherland's MIT notes, Lecture 18](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_full_lec.pdf).
This statement has `u` fixed and supplies no useful threshold at a
prescribed power of `U`.

Assign one chosen prime `p_I` to each partition and prescribe any positive
integer valuation `e_I`. Its quadratic `M_I` has two simple projective
roots modulo `p_I`, at least one in the chart `(s:t)=(1:r)`. At either
root every other minor is a unit. The degenerate conic consists of two
nonparallel lines defined over `F_p`, so its quadratic-part discriminant
satisfies

```text
Delta = -(nonzero square) mod p_I.
```

Since `p_I=1 mod 4`, the circle cover has two nonsingular local lifts
there. Hensel lifting for the simple minor root allows one to prescribe
`v_p(M_I)=e_I` exactly by a residue modulo `p_I^(e_I+1)`; throughout this
residue neighborhood `Delta` remains a nonzero square in `Q_p`.

CRT now gives a single rational base parameter with these assigned
valuations at all ten primes and with all nonassigned minors units at
each of them. It can also lie in any chosen real circle interval. For an
elementary construction, let `P` be the product of the prescribed prime
powers, choose a large denominator `q=1 mod P`, and choose a numerator in
the required residue class modulo `P` within distance `P` of `q` times a
point of that interval. Reduction to lowest terms does not change the
specified valuations, because the denominator is a unit at every chosen
prime. Homogeneity of degree twelve preserves the local squareclass of
`Delta` when passing from `r` to the primitive integer pair.

If the actual central vector comes with a rational seed point on its
circle cover, one can simultaneously require the parameter to remain
sufficiently close to that seed at any additional finite set of primes,
including every prime dividing the weights or subset sums. The same CRT
argument with a fixed extra denominator gives this finite weak
approximation. Choose each local square root near the seed square root;
the valuations of any specified finite set of nonzero node or chord
functions then stay unchanged. Thus existing local singleton and pair
profiles at the seed can also be retained. This still gives a collection
of local lifts, not one global rational lift.

This produces the actual local circle-cover lifts at every selected
prime, as well as a real lift. It does **not** produce a common rational
square root of `Delta`, assert solubility at every other prime, or control
the relative sizes of the weights and selected primes. In particular it
does not construct a full fair circle core. It proves the narrower useful
point: the ten splitting restrictions and their quadratic-reciprocity
relations are jointly compatible with arbitrarily prescribed local
balanced-cut valuations. Global rational lifting and the fair height
requirements remain additional conditions.
