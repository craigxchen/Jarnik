# Independent audit of the disjoint-support elliptic six-curve

The construction in
[disjoint_mixed_elliptic_six_curve.md](disjoint_mixed_elliptic_six_curve.md)
passes the independent geometric audit. Its normalization is a
geometrically connected genus-one curve over `Q` with rational points.
The six-point map is birational onto its image, and the twenty-five
stable boundary divisors pull back to twenty-five **distinct** rational
points, each with multiplicity one. The stronger assertion about
distinct support is valid here.

## Degree, connectedness, and the coordinate field

The homogeneous cross-product of two distinct pencil members is

```text
(mu-lambda)(x+z)(x-3z)(x-4z)(x-5z).
```

Each pencil member has numerator leading coefficient one. Its
denominator is nonzero at all four displayed roots. Consequently it
has no numerator-denominator cancellation, finite or infinite, and
defines a degree-two morphism of projective lines. The cross-product
also has nonzero value at `[1:0]`; there is no fifth common argument
at infinity. Explicitly, the three values there are
`infinity,-5329/240,257/8`, all different.

The supplied exact critical-point calculation verifies the three
branch pairs `{u,v}`, `{u,s}`, `{v,t}` with four distinct values.
There are no omitted critical points at poles: the denominators have
simple roots, and the derivative numerator has degree two with the
two listed roots. Infinity is also unramified for each map.

The quadratic function-field classes are independent over `Qbar(T)`:
valuation at `s` detects the second class, valuation at `t` the third,
and valuation at `u` then detects the first. Hence their compositum
has degree eight and its smooth normalization is geometrically
connected. Each of the four inertia groups has order two, including
at the two shared branch values; normalization does not turn shared
quadratic ramification into ramification of order four. Thus

```text
2g(E)-2 = -2*8 + 4*(8/2) = 0.
```

Each common target has a rational simple inverse root in every
quadratic map, hence two distinct rational inverse roots in each.
The unramified fiber therefore consists of eight rational points.
There are thirty-two such points across the four distinct targets.
Choosing one as origin makes the source an elliptic curve over `Q`.

By construction `Q(E)=Q(a,b,c)` since `T=g_a(a)`. On the open moduli
chart the three fixed labels `infinity,0,1` determine the projective
normalization uniquely, so `a,b,c` are coordinates on that chart.
The source map is therefore birational onto its image. Properness of
the stable moduli space extends it across the finitely many collisions.

## Exhaustiveness, distinctness, and multiplicity of contacts

Here are the *other* inverse roots, after the source normalization
`X=(x+1)/(16-4x)`. All entries are exact.

| Common root | Other root of `g_a` | Other root of `g_b` | Other root of `g_c` |
|---|---|---|---|
| `0` | `3/28` | `-22173/297052` | `227/252` |
| `1` | `1/4` | `2209/42436` | `361/324` |
| `infinity` | `7/32` | `15463/334688` | `2527/2752` |
| `-3/2` | `9/44` | `19881/466796` | `361/396` |

Within every row the alternatives are pairwise distinct and avoid
the three anchors and the common root. Thus a binary inverse-root
choice produces exactly one collision cluster, or none. In each
anchored fiber the seven nonempty subsets of moving labels give
seven different source points. In the nonanchored fiber precisely
the three pairs and the triple give contacts. Fibers of different
targets are disjoint. This proves `7+7+7+4=25` distinct supports,
not just twenty-five incidences counted with repetition. Their
partitions exhaust the `15+10` stable boundary divisors.

Any equality between moving coordinates is a zero of the homogeneous
cross-product and hence has one of the four listed common arguments.
An anchor equality must be at its listed common target. Consequently
there are no extra collisions at denominator poles, at the pole of
the normalizing coordinate change, at target infinity, or in a
ramified fiber. In particular the branch fibers cannot contribute:
their targets differ from every possible collision target.

At each listed contact `T-T_0` is a source uniformizer. In a local
coordinate around the common argument, including `1/X` at infinity,
all inverse slopes are nonzero. The simple cross-product zero makes
them pairwise distinct. The anchor has slope zero when present.
One blowup therefore separates every marked point in the cluster;
there are no further nested clusters. Its node smoothing parameter
has order one. Every boundary pullback is exactly the corresponding
single reduced point.

The original checker passed. The additional independent checker
[check_disjoint_mixed_elliptic_six_curve_audit.py](check_disjoint_mixed_elliptic_six_curve_audit.py)
compares polynomial coefficients, checks degree and infinity behavior,
enumerates the actual thirty-two coordinate triples, verifies all
twenty-five distinct contacts, and checks the inverse slopes exactly.

## Arithmetic scope

The repeated-support obstruction of the first mixed example does not
apply. Distinct rational boundary points also remain distinct away
from finitely many primes after spreading out the fixed source and
map. Thus generic good-prime boundary contacts have pairwise disjoint
support. This observation does not prescribe their sizes or their
split-versus-inert prime content.

No arithmetic orbit is constructed by this geometric audit itself.
The separate
[integer contact-profile lemma](disjoint_mixed_elliptic_integer_contact_profile.md)
proves that compact local recurrence realizes twenty-five pairwise
coprime equal-scale unoriented integer cores with bounded pair
determinant remainders. The nontorsion point exists by the exact
positive-rank certificate in Section 4 of the construction. The lemma bounds
the bad-prime and archimedean terms rather than inferring finite contact
sizes from divisor degrees alone.

The remaining arithmetic conditions include negligible inert-prime
contributions, separation of each unoriented contact into the two
oriented Gaussian core blocks, private singleton factors, and
simultaneously small primitive residues, row corrections, and arc
width after a common lift. The global uniform `C sqrt(R)` lattice-arc
bound remains unproved.
