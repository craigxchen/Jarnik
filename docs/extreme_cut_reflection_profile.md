# Extreme-cut cost for pair reflections in a fair allocation profile

This note records a profile-level consequence of the exact prime-power
reflection denominator. It is separate from the radial-count obstruction for
full reflection unions in
[the radial-count note](reflection_union_radial_count_obstruction.md).

Let `Z=(z_1,...,z_k)`, with `k>=3`, be a primitive Gaussian tuple of
common modulus `R`. Put `N=R^2`. At a split prime `p=pi bar(pi)`, write

    e_i=v_pi(z_i),       min_i e_i=0,       max_i e_i=W_p.

For distinct anchors `a,b`, the reflection coefficient is
`rho_ab=z_a z_b/N`. Its reduced Gaussian denominator has norm exponent

    q_ab,p = |e_a+e_b-W_p|.                              (1)

Thus
`Q_ab=Norm(d_ab)=product_p p^q_ab,p` and
`log Q_ab=sum_p q_ab,p log p`. For a restricted subset `J`, first
divide its common Gaussian divisor; then (1) becomes

    q^J_ab,p = |e_a+e_b-u_J-v_J|,

where `u_J=min_(i in J)e_i` and `v_J=max_(i in J)e_i`.

## Extreme-cut lemma

For a threshold `1<=t<=W_p`, let `S_t={i:e_i>=t}`.

If `S_t={a,b}`, every other exponent is below `t`, so the maximum
exponent belongs to `a,b`. Hence

    q_ab,p=e_a+e_b-W_p=min(e_a,e_b)>=t.

If `S_t` is the complement of `{a,b}`, primitivity forces the minimum
exponent into `a,b`, and

    q_ab,p=W_p-e_a-e_b=W_p-max(e_a,e_b).

In either orientation, the number of threshold layers producing the
unoriented extreme cut `C_ab={{a,b}, [k]\{a,b}}` is at most `q_ab,p`.
Both orientations cannot occur at one prime: threshold sets are nested,
whereas these two nonempty complementary sets are not nested.
If `H_ab` denotes the total layer weight of this cut, counting `log p`
per threshold layer, then

    log Q_ab >= H_ab.                                    (2)

Write `W=sum_p W_p log p=log N`. If every nontrivial unoriented `k`-cut has
weight in

    [(1-eta)W/(2^(k-1)-1), (1+eta)W/(2^(k-1)-1)],

then every pair satisfies

    log Q_ab >= (1-eta)W/(2^(k-1)-1).                    (3)

For the eight-row profile this is `(1-eta)W/127`. The stronger roughly
`W/2` estimate holds for independent one-layer blocks, but cannot be deduced
from the aggregate cut profile because same-prime threshold layers can cancel
inside the absolute value in (1).

The coefficient `1/(2^(k-1)-1)` is sharp for this information. Fix
`a,b`, put `T=[k]\{a,b}`, and give every cut class formal weight
`h`. For each complementary pair of nonempty proper subsets `A,B` of `T`,
use one block with allocations `0` on `A`, `1` on `{a,b}`, and `2`
on `B`, assigning weight `h` to each layer. Its two layers produce the cut
classes represented by `A` and `B`, while `q_ab,p=0`. The remaining
agreeing class represented by `A=T` is supplied by a one-layer prime with
allocations `1` on `{a,b}` and `0` on `T`; it contributes exactly
`h` to `log Q_ab`. Every separating cut class gets its own one-layer
prime with `a,b` on opposite sides and contributes zero. Thus all
`2^(k-1)-1` unoriented cut classes have weight `h`. Writing
`D_ab=sum_p |e_a-e_b| log p` for the primitive pair conductor gives

    log Q_ab=h,    W=(2^(k-1)-1)h,    D_ab=2^(k-2)h,
    log Q_ab=2D_ab-W.

This is an exact formal allocation. Choose distinct split primes for
the finitely many blocks and multiply each block's displayed allocations
by `t_p=floor(h/log p)`. For fixed `k` and these fixed primes, every
cut weight is `h+O_k(1)` and `log Q_ab=h+O_k(1)` as `h` tends to infinity.
The resulting equal-norm Gaussian tuples are primitive. Their angles
are uncontrolled; no short-arc or endpoint realization is claimed.

For a retained subset `J` of size `ell>=3`, the same proof after subset
primitive normalization gives `log Q^J_ab>=H^J_ab`. Under a uniform
`k`-cut profile,

    W_J = 2^(k-ell)(2^(ell-1)-1) W/(2^(k-1)-1),
    H^J_ab = 2^(k-ell) W/(2^(k-1)-1).

Thus `H^J_ab=W_J/(2^(ell-1)-1)`. For fixed `k`, any retained subset of
at least three rows has a positive proportional denominator cost; that
proportion need not be bounded below when `k` grows. The case `ell=2`
is degenerate: the reflection exchanges the two retained rows and
`Q^J_ab=1`.

## Cancellation fixture

Counting agreeing threshold layers is false in general. At one prime take

    (e_0,e_1,e_2,e_3)=(1,1,0,2),   W_p=2.

For the anchor pair `(0,1)`, the two threshold cuts are
`{0,1,3}` and `{3}`. The anchors agree on both layers, but

    q_01,p=|1+1-2|=0.

This is a literal nested-prime allocation with minimum exponent `0` and
maximum exponent `2`; it is a cancellation fixture, not an assertion that
the formal allocation is realized by a short-arc endpoint tuple.

## Clearing denominators and radius bookkeeping

Write the reduced coefficient as `rho_ab=A/d`, with
`gcd_G(A,d)=1`. The least Gaussian integer multiplier making the
reflected tuple integral is `d`. If `d=x+iy`, the least positive integer
multiplier is `Norm(d)/gcd(|x|,|y|)`. The canonical pair denominator has
primitive coordinates, so this integer is `Q_ab=Norm(d)`.

Multiplying the reflected tuple alone gives `A bar(Z)`. Its common Gaussian
gcd is `A`, so primitive normalization returns (up to a unit) to
`bar(Z)`, with the original radius `R`. Therefore `R sqrt(Q_ab)` is
the radius cost for the cleared source-plus-reflection union

    d Z union A bar(Z),

whose gcd is one. Before primitive normalization, the reflected-only cleared
arc has radius `R sqrt(Q_ab)` and normalized constant
`C Q_ab^(1/4)`; after normalization its intrinsic radius and constant
return to `R` and `C`. A positive integer clearing multiplier has the
analogous factor `sqrt(m)` in the raw endpoint constant.

The lower bound `log Q_ab>=H_ab` alone does not force the union
endpoint constant to diverge when the source constant is allowed to shrink.
For distinct actual directions, the ordinary nonzero Gaussian pair
residue gives `Delta>=exp(-D_ab/2)`, with units retained. Indeed, write
`z_b/z_a=epsilon H/bar(H)` with `Norm(H)=exp(D_ab)`. Its normalized
chord is `|epsilon H-bar(H)|/|H|`; the numerator is a nonzero Gaussian
integer. Angular span dominates this chord. Combining with the union
radius gives only

    C_union>=exp((W+log Q_ab-2D_ab)/4).

At the sharp formal profile, `log Q_ab=h=2D_ab-W`, so this exponent
is exactly zero. The separate radial-count theorem gives a stronger
count-dependent restriction for actual symmetric unions. Neither
calculation proves a bound for the original unsymmetrized cluster.

[The exact checker](check_extreme_cut_reflection_profile.py) exhaustively
tests the local inequality for three through seven rows with prime
width four, enumerates restricted-cut multiplicities, and verifies the
sharpness construction and cancellation fixture.
