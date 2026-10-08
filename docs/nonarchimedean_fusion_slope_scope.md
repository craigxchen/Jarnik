# What the fusion filtration does not yet give the integer lattice

The coalition calculation fixes the order of the **global determinant
section** of the coherent conformal-block kernel.  It does not by
itself select an integer vector shorter than the determinant average.
For `(n,d)=(4,2)` the [tensor coefficient floor](higher_rank_tensor_multiplicity_floor.md)
reaches that average exactly, leaving no rank-one slope to separate.
For `(n,d)=(3,4)` a strict first-minimum gain would need a global
integral sublattice whose individual height is controlled across the
entire full core profile.  No such sublattice follows from the local
fusion decomposition alone.

## The local kernel has trivial inclusion Smith invariants

Let `R` be the discrete valuation ring of one Gaussian core prime,
and let `E:R^r -> R^N` be the coherent evaluation matrix at an actual
profile, or a polynomial one-parameter collision matrix.  Its kernel
`K` is saturated: `R^r/K` is isomorphic to the image of `E`, a
torsion-free finite module over a DVR, hence is free.  Therefore the
inclusion `K subset R^r` splits, and **all its Smith exponents are
zero**.  This remains true after dividing output rows by their forced
core powers, since a diagonal matrix of nonzero divisors does not
change the kernel over `R`.

The nonzero coalition order `nu_T` in
[the fusion-content note](terminal_kernel_fusion_content.md) has a
different meaning.  A primitive polynomial Pluecker tuple `W` on the
all-distinct parameter space is a global section of `det K`.  Along a
collision path it can equal `t^nu_T` times a **local** nonvanishing
frame of that line.  Primitivity removes a common polynomial factor
on the original parameter space; it does not prevent vanishing on the
higher-codimension collision center.  In the six-row rank-one case,
the matching-vector numerator already has positive majority orders
at large coalitions while the local kernel inclusion stays saturated.
Thus neither a full-rank row-normalized special fiber nor the Smith
form of `K subset R^r` contradicts positive `nu_T`.

## A precise obstruction to deducing a first minimum from local data

For every integer `N>=2`, put `D=N^2+1` and consider the index-`D`
lattices in `Z^2`

```text
L_axis = {(x,y): x=0 mod D},
L_tilt = {(x,y): x+N*y=0 mod D}.
```

Both have determinant `D` and quotient `Z/D`, hence the same Smith
exponents `(0,v_p(D))` at **every** rational prime `p`.  Yet their
Euclidean first minima are respectively `1` and `sqrt(D)`.
For the second assertion, `(-N,1)` has length `sqrt(D)`.  A nonzero
vector with smaller length has `|x|<=N`, `|y|<=N`.  If `|y|<N`, then
`|x+Ny|<=N^2<D`, so the congruence forces `x=-Ny` and its length is
at least `sqrt(D)` unless it is zero.  If `|y|=N`, the strict length
bound forces `x=0`, but `N^2` is not divisible by `D`.  Thus no
shorter vector exists.

The abstract local elementary divisors and determinant are identical;
the embedding of the local flag in the fixed coefficient lattice
changes the first minimum by a factor of order `sqrt(D)`.  The example
does not model a circle configuration.  It isolates the logical step
missing from any inference that a fusion channel with favorable local
order yields a short **actual integral** invariant relation.

## Twelve-row local-generator diagnostic

The deterministic [finite checker](check_nonarchimedean_fusion_slope_scope.py)
uses the existing 462 standard four-triangle source vectors and the
4620 localized-six generators over `F_101`.  At the six-label
collision `T={0,...,5}`, rows in `T` have slopes `t*(i+1)` and the
other rows have slopes `i+17`.  A generator supported on six labels
`I` has forced order `max(0,|I intersection T|-3)`.  The exact
finite-field ranks are

```text
raw specialization rank:                      311,
rank after including order-one leading terms: 341,
rank after all normalized leading terms:      341.
```

The row-normalized coherent matrix has rank `121`, so the generic
kernel rank is `462-121=341` in this finite-field sample.  Since
the raw generator matrix drops rank by `30` and 30 order-one leading
columns restore full rank, its selected presentation has local Smith
multiset `311` zeros and `30` ones over `F_101[[t]]`.  These are
**presentation slopes** for the localized-six generating family.
The global primitive kernel Pluecker order at a size-six coalition is
`nu_6=15`; it compares a different global determinant trivialization
with the local saturated frame.  The numbers `30` and `15` therefore
must not be identified.  The finite ranks are a specified
characteristic-101 fixture, not a proof of an all-residue Smith form.

## The missing gluing statement for a strict exponent gain

At `(n,d)=(4,2)`, the conformal-block kernel has rank one.  The
full-core determinant upper exponent is `B/h=116`, and the tensor
coefficient floor is also `116`.  Hence, under the same full-profile
and fixed coefficient-conversion hypotheses, the primitive relation
has leading height `116w+o(w)`.  There is no second individual slope
or shorter vector in that kernel.

At `(n,d)=(3,4)`, the kernel has rank `341`, determinant exponent
`B=171600`, and average `B/h=171600/341>503`; the all-coordinate
coefficient floor remains `90`.  One sufficient refinement would
construct a **globally defined rational subbundle** `F` of positive
rank `a` inside the coherent kernel, together with a saturated
integral lattice `F(P) intersection Z^462` for every actual profile,
such that its primitive covolume has a full-core upper exponent
`B_F<90*a`.  Then Minkowski would produce a nonzero relation below
the coefficient floor.  Local fusion channels at one coalition give
candidate fiber subspaces, but they vary with the coalition, and the
factorization theorem supplies neither a common rational subbundle
nor bounded integral transition matrices among those flags.  Without
that gluing and height control, their individual residues cannot be
summed across the disjoint Gaussian core primes for one lattice
vector.  This gives a sufficient additional lemma for this route;
it neither proves that lemma nor rules out other global constructions.

The [generic local sharpness construction](local_tensor_coefficient_floor_sharpness.md)
attains the tensor coefficient floor at each chosen coalition and
color occupancy.  It reinforces that any stronger exclusion would
need simultaneous information across the core primes or special
actual-residue constraints, rather than another pointwise occupancy
count.
